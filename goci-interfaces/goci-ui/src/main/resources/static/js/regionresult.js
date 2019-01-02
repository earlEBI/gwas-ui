/** DRY. From Xin original code. We must refactor all these 'action'result.js in a common way! */

var global_fl;
var global_raw;
global_fl = 'pubmedId,title,author_s,orcid_s,publication,publicationDate,catalogPublishDate,authorsList,' +
    'initialSampleDescription,replicateSampleDescription,ancestralGroups,countriesOfRecruitment,' +
    'ancestryLinks,' + 'fullPvalueSet,' + 'genotypingTechnologies,' + 'authorAscii_s,' +
    'association_rsId,' + //size per study
    'traitName,mappedLabel,mappedUri,traitUri,shortForm,' +
    'label,' + 'efoLink,parent,id,resourcename,';
global_fl = global_fl + 'riskFrequency,qualifier,pValueMantissa,pValueExponent,snpInteraction,multiSnpHaplotype,rsId,strongestAllele,context,region,entrezMappedGenes,reportedGene,merged,currentSnp,studyId,chromosomeName,chromosomePosition,chromLocation,positionLinks,author_s,publication,publicationDate,catalogPublishDate,publicationLink,accessionId,initialSampleDescription,replicateSampleDescription,ancestralGroups,countriesOfRecruitment,numberOfIndividuals,traitName_s,mappedLabel,mappedUri,traitUri,shortForm,labelda,synonym,efoLink,id,resourcename,range,orPerCopyNum,betaNum,betaUnit,betaDirection'
global_raw = 'fq:resourcename:association or resourcename:study';

// Gene page specific constans:

/**
 * Other global setting
 */
var pageRowLimit=5;

$(document).ready(() => {

    $('html,body').scrollTop(0);
    var searchTerm = getTextToSearch('#query');

    // console.log("Query:" + searchTerm);
    // console.log("Loading search module!");
    if (searchTerm != '') {

        // console.log("Start search for the text " + searchTerm);
        var elements = {};
        searchTerm.split(',').forEach((term) => {
            elements[term] = term;
        });

        // console.log(elements);
        executeQuery(elements, true);
    }
});

/**
 * The elem to search is defined by the url, as a main entry of the page. It is stored in the div id
 * in the date attribute of the <global_elem_info_tag_id>`
 * @return Eg. String efoID - 'EFO_0000400'
 * @example getElemToSearch()
 */

getTextToSearch = function(divId){
    return $(divId).text();
}

executeQuery = function(data={}, initLoad=false) {
    // console.log("executeQuery");
    updatePage(initLoad);
}



updatePage = function(initLoad=false) {

    //start spinner. The spinner will be stopped whenever the data is ready, thus closed by the corresponding data loading function.
    if(initLoad){
        showLoadingOverLay('#summary-panel-loading');
    }
    showLoadingOverLay('#study-table-loading');
    showLoadingOverLay('#association-table-loading');
    
    var main = getTextToSearch('#query');
    
    //******************************
    // when solr data ready, process solr data and update badges in the selection cart
    //******************************
    var solrPromise = getDataSolr(main, initLoad);
    
}


/*
Extracting data from the fat solr based on the query term
 */
function getDataSolr(main, initLoad=false) {
    // initLoad will be pass to processEfotraitData, controlling whether to upload the triat information(initload)
    // or just reload the tables(adding another efo term)
    
    var searchQuery = main;
    var solrQuery = '';
    var regionTest = /([XY0-9]{1,2}):(\d+)-(\d+)/gi; // matches regions 6:234511-23500
    var cytobandTest = /([XY0-9]{1,2})([PQ][0-9]+\.[0-9]+)/gi; // matches cytobands eg 6p33.1

    // console.log("Solr research request received for " + searchQuery);

    // Testing if the query was a cytoband or region:
    if ( searchQuery.match(cytobandTest) ){
        solrQuery = searchQuery;
    }
    else if (searchQuery.match(regionTest) ){
        var coordinates = regionTest.exec(searchQuery);
        solrQuery =  "chromosomeName: "+coordinates[1]+" AND chromosomePosition:[ "+coordinates[2]+" TO "+coordinates[3]+" ]"
    }
    else {
        $('#lower_container').html("<h2>The provided query term <em>"+searchQuery+"</em> cannot be interpret as a region.</h2>");
    }

    // console.log("** solr Query: "+ solrQuery)

    return promisePost( gwasProperties.contextPath + 'api/search/advancefilter',
        {
            'q': solrQuery,
            'max': 99999,
            'group.limit': 99999,
            'group.field': 'resourcename',
            'facet.field': 'resourcename',
            'hl.fl': 'shortForm,efoLink',
            'hl.snippets': 100,
            'fl' : global_fl == undefined ? '*':global_fl,
            'raw' : global_raw == undefined ? '' : global_raw,
        },'application/x-www-form-urlencoded').then(JSON.parse).then(function(data) {

        // Check if Solr returns some results
        if ( data.grouped.resourcename.groups.length == 0 ) {

            // console.log(data)
            $('#lower_container').html("<h2>No associaitons could be found in the region: <em>"+searchQuery+"</em> in the GWAS Catalog database</h2>");
        }
        else {
            processSolrData(data, initLoad, searchQuery);
            setDownloadLink(searchQuery);
        }

        // console.log("Solr research done for " + searchQuery);
        return data;
    }).catch(function(err) {

        console.error('Error when seaching solr for' + searchQuery + '. ' + err);
        throw(err);
    })
    
}


function processSolrData(data, initLoad=false, searchTerm) {
    var isInCatalog=true;
    
    data_association = [];
    data_study = [];
    data_association.docs = [];
    data_study.docs = [];
    
    if (data.grouped.resourcename.matches == 0) {
        isInCatalog = false;
    }

    //split the solr search by groups
    //data_study, data_association
    data_facet = data.facet_counts.facet_fields.resourcename;
    data_highlighting = data.highlighting;

    $.each(data.grouped.resourcename.groups, (index, group) => {
            switch (group.groupValue) {
        case "efotrait":
            data_efo = group.doclist;
            break;
        case "study":
            data_study = group.doclist;
            break;
        case "association":
            data_association = group.doclist;
            break;
            //not sure we need this!
        case "diseasetrait":
            data_diseasetrait = group.doclist;
            break;
        default:
        }
    });
    
    //remove association that annotated with efos which are not in the list
    var remove = Promise.resolve();

    remove.then(()=>{
        // If no solr return, greate a fake empty array so tables/plot are empty
        if( !isInCatalog ){
            data_association.docs = [];
            data_study.docs = [];
        }

        var PAGE_TYPE = "region";

        displayDatatableAssociations(data_association.docs);
        console.log("[Info] displayDatatableAssociations - OK");

        generateGeneInformationTable(searchTerm, data_study);
        console.log("[Info] generateGeneInformationTable - OK");

        // when chr:pos is queried, there's no returned study. We have to specifically fetch those.
        if ( data_study.docs.length == 0 ){
            console.log("** There's no sudy!!!");
            data_study.docs = fetchStudies(data_association.docs);
            displayDatatableStudies(data_study.docs, PAGE_TYPE);
        }
        else {
            displayDatatableStudies(data_study.docs, PAGE_TYPE);
        }
        console.log("[Info] displayDatatableStudies - OK");

        checkSummaryStatsDatabase(data_study.docs);
        console.log("[Info] checkSummaryStatsDatabase - OK");

    })

}

// Helper function to retrieve Ensembl data through the REST API
// SYNC!!
function getEnsemblREST( URL ){
    var result = null;
    $.ajax({
        url: URL,
        type: 'get',
        dataType: 'json',
        async: false,
        success: function(data) {
            result = data;
        }
    });
    return result;
}

/**
 * This function fills up the gene table.
 * Input:
 *    Gene name      -  searchTerm
 *    Study document -  data_study.docs
 *
 * Behaviour:
 *    Fills up gene information table
 *    1) Extracts data from Ensembl based on gene name.
 *    2) Extracts cross reference data from Ensembl.
 *    3) Extracts reported traits from study documents.
 */
function generateGeneInformationTable(geneName, studies) {
    // Extracting gene data from Ensembl:
    // var geneQueryURL = gwasProperties.EnsemblRestBaseURL + "/lookup/symbol/homo_sapiens/" + geneName + "?content-type=application/json"
    // var geneData = getEnsemblREST(geneQueryURL);

    // adding gene data to html:
    // $("#geneSymbol").html(`${geneData.display_name}`);
    // var description = geneData.description.split(" [S")[0];
    // $("#description").html(`${description}`)
    // $("#genomicCoordinates").html(`${geneData.seq_region_name}:${geneData.start}-${geneData.end}`);
    // $("#biotype").html(`${geneData.biotype.replace("_", " ")}`);
    //
    // console.log(studies.length)
    //
    // // Looping through all studies and parse out repoted genes:
    // var reportedTraits = {};
    // for ( var study of studies.docs) {
    //     reportedTraits[study.traitName_s] = 1;
    // }
    //
    // // joining reported traits & sort:
    // reportedTraits = Object.keys(reportedTraits).sort()
    // var joinedTraits = reportedTraits.join("</li>\n\t<li>")
    // $("#reportedTraits").html(`<ul>\n\t<li>${joinedTraits}</li></ul>`);
    // console.log(joinedTraits)
    //
    // // Extracting cross-references:
    // var xrefQueryURL = gwasProperties.EnsemblRestBaseURL + '/xrefs/id/' + geneData.id + '?content-type=application/json'
    // var xrefData = getEnsemblREST(xrefQueryURL);
    // var entrezID = "NA";
    // var OMIMID = "NA";
    // for ( xref of xrefData ){
    //     if ( xref.dbname == "EntrezGene" ){
    //         entrezID = xref.primary_id
    //     }
    //     if ( xref.dbname == 'MIM_GENE' ){
    //         OMIMID = xref.primary_id
    //     }
    // }
    //
    // // Adding automatic cross references pointing to Ensembl:
    // $("#ensembl_button").attr('onclick', "window.open('"+gwasProperties.EnsemblURL+"Summary?db=core;g="+geneData.id+"',    '_blank')");
    // $("#ensembl_phenotype_button").attr('onclick', "window.open('"+gwasProperties.EnsemblURL+"Phenotype?db=core;g="+geneData.id+"',    '_blank')");
    // $("#ensembl_pathway_button").attr('onclick', "window.open('"+gwasProperties.EnsemblURL+"Pathway?db=core;g="+geneData.id+"',    '_blank')");
    // $("#ensembl_regulation_button").attr('onclick', "window.open('"+gwasProperties.EnsemblURL+"Regulation?db=core;g="+geneData.id+"',    '_blank')");
    // $("#ensembl_expression_button").attr('onclick', "window.open('"+gwasProperties.EnsemblURL+"ExpressionAtlas?db=core;g="+geneData.id+"',    '_blank')");
    //
    // // Adding automatic cross reference pointing to Open targets:
    // $("#opentargets_button").attr('onclick', "window.open('"+gwasProperties.OpenTargetsURL+ geneData.id+"',    '_blank')");
    //
    // // Looping through the cross references and extract entrez id:
    // if ( entrezID != "NA" ){
    //     $("#entrez_button").attr('onclick', "window.open('"+gwasProperties.EntrezURL+ entrezID + "',    '_blank')");
    // }
    // // Looping through the cross references and extract OMIM id:
    // if ( OMIMID != "NA" ){
    //     $("#OMIM_button").attr('onclick', "window.open('"+gwasProperties.OMIMURL+ OMIMID + "',    '_blank')");
    // }

    // Print out some info to make sure things are not messed up completely:
    // console.log("[Info] Number of reported traits:" + reportedTraits.length)
    // console.log("[Info] ID: " + geneData.id);
    // console.log("[Info] Biotype: " + geneData.biotype);
    // console.log("[Info] Description: " + geneData.description);
    // console.log("[Info] Genomic location: " + geneData.seq_region_name + ":" + geneData.start + "-" + geneData.end)

    // OK, loading is complete:
    hideLoadingOverLay('#summary-panel-loading');
}


// This function returns studies based on GCSC accession id extracted from association documents:
function fetchStudies(associationDocs){
    // Look through all the docs and get accessions:
    var accessionIDs = [];
    var stepSize = 20;
    var studyData = [];

    for (var assoc of associationDocs){
        accessionIDs.unshift(assoc.accessionId)
    }

    accessionIDs = Array.from(new Set(accessionIDs))

    // Generate unqiue list of accessions:
    // Loop through the accessions and download data by 50 studies at a time: (might need to be an other function for that)
    for (var i = 0; i < accessionIDs.length; i += stepSize) {
        temparray = accessionIDs.slice(i, i + stepSize);
        console.log("** Looping through the accession IDs. Chunk: ", i);
        studyData = studyData.concat(getStudyData(temparray))
    }
    return(studyData)
};

// Query slim solr to return rsIDs that are mapped to a given gene:
// WARNING: syncronous call!!
function getStudyData(studyIDs){
    var queryString = studyIDs.join(" OR ")
    var result = null;
    console.log("** queried accessions: " + queryString)
    $.ajax({
        url : gwasProperties.contextPath + 'api/search/advancefilter',
        data : {'q': "( " + queryString + ') AND resourcename:study'},
        type: 'get',
        dataType: 'json',
        async: false,
        success: function(data){
            // Parsing out response:
            result = data.grouped.resourcename.groups[0].doclist.docs
        }
    });
    console.log(result)
    return result;
}
