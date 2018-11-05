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
const EnsemblRestBaseURL = "https://rest.ensembl.org";
const EnsemblURL = "http://www.ensembl.org/Homo_sapiens/Gene/";
const OpenTargetsURL = "https://www.targetvalidation.org/target/";
const EntrezURL = "https://www.ncbi.nlm.nih.gov/gene/"
/**
 * Other global setting
 */
var pageRowLimit=5;

$(document).ready(() => {

//jump to the top of the page
    $('html,body').scrollTop(0);

var searchTerm = getTextToSearch('#query');
console.log("Query:" + searchTerm);
console.log("Loading search module!");
if (searchTerm != '') {
    // console.log("Start search for the text " + searchTerm);
    var elements = {};
    searchTerm.split(',').forEach((term) => {
        elements[term] = term;
    });
    //first load
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
    
    //start spinner. The spinner will be stoped whenever the data is ready, thus closed by the coresponding data loading function.
    if(initLoad){
        showLoadingOverLay('#summary-panel-loading');
//            showLoadingOverLay('#highlight-study-panel-loading');
    }
    showLoadingOverLay('#study-table-loading');
    showLoadingOverLay('#association-table-loading');
    
    var main = getTextToSearch('#query');
    
    //******************************
    // when solr data ready, process solr data and update badges in the selection cart
    //******************************
    var solrPromise = getDataSolr(main, initLoad);
    
}


/**
 * Make solr query.
 * @param {String} mainEFO
 * @param {[]String} additionalEFO
 * @param {[]String} descendants
 * @param {Boolean} initLoad
 * @returns {Promise}
 */
function getDataSolr(main, initLoad=false) {
    // initLoad will be pass to processEfotraitData, controlling whether to upload the triat information(initload)
    // or just reload the tables(adding another efo term)
    
    var searchQuery = main;
    
    console.log("Solr research request received for " + searchQuery);
    return promisePost( gwasProperties.contextPath + 'api/search/advancefilter',
        {
                'q': 'entrezMappedGenes : ' + searchQuery + ' OR association_entrezMappedGenes : ' + searchQuery,
            'max': 99999,
            'group.limit': 99999,
            'group.field': 'resourcename',
            'facet.field': 'resourcename',
            'hl.fl': 'shortForm,efoLink',
            'hl.snippets': 100,
            'fl' : global_fl == undefined ? '*':global_fl,
            // 'fq' : global_fq == undefined ? '*:*':global_fq,
            'raw' : global_raw == undefined ? '' : global_raw,
        },'application/x-www-form-urlencoded').then(JSON.parse).then(function(data) {
        // Check if Solr returns some results
        if (data.grouped.resourcename.groups.length == 0) {
            $('#lower_container').html("<h2>The Gene name <em>"+searchQuery+"</em> cannot be found in the GWAS Catalog database</h2>");
        }
        else {
            processSolrData(data, initLoad, searchQuery); // gene name is now added to the process solr data function.
            //downloads link : utils-helper.js
            setDownloadLink(searchQuery);
        }
        // console.log("Solr research done for " + searchQuery);
        return data;
    }).catch(function(err) {
        // console.error('Error when seaching solr for' + searchQuery + '. ' + err);
        throw(err);
    })
    
}

/**
 * Parse the Solr results and display the data on the HTML page
 * @param {{}} data - solr result
 * @param {Boolean} initLoad
 */

/**
 * Parse the Solr results and display the data on the HTML page
 * @param {{}} data - solr result
 * @param {Boolean} initLoad
 */
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
    //TODO not repeat yourself!!!!
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
        //If no solr return,greate a fake empyt array so tables/plot are empty
        if(!isInCatalog) {
        data_association.docs = []
        data_study.docs = []
    }

    var PAGE_TYPE = "gene";
    
    //update association/study table
    displayDatatableAssociations(data_association.docs);
    console.log("[Info] displayDatatableAssociations - OK")
    displayDatatableStudies(data_study.docs, PAGE_TYPE);
    console.log("[Info] displayDatatableStudies - OK")
    checkSummaryStatsDatabase(data_study.docs);
    console.log("[Info] checkSummaryStatsDatabase - OK")
    generateGeneInformationTable(searchTerm, data_study)
    console.log("[Info] generateGeneInformationTable - OK")
    //displaySummaryPublication(data_study.docs);
    
})

}

// Helper function to retrieve Ensembl data through the REST API
// SYNC!!
function getEnsemblREST( URL )
{
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
    var geneQueryURL = EnsemblRestBaseURL + "/lookup/symbol/homo_sapiens/" + geneName + "?content-type=application/json"
    var geneData = getEnsemblREST(geneQueryURL);

    // adding gene data to html:
    $("#geneSymbol").html(`${geneData.display_name}`);
    var description = geneData.description.split(" [S")[0];
    $("#description").html(`${description}`)
    $("#genomicCoordinates").html(`${geneData.seq_region_name}:${geneData.start}-${geneData.end}`);
    $("#biotype").html(`${geneData.biotype.replace("_", " ")}`);

    console.log(studies.length)

    // Looping through all studies and parse out repoted genes:
    var reportedTraits = {};
    for ( var study of studies.docs) {
        reportedTraits[study.traitName_s] = 1;
    }

    // joining reported traits:
    var joinedTraits = Object.keys(reportedTraits).join("</li>\n\t<li>")
    $("#reportedTraits").html(`<ul>\n\t<li>${joinedTraits}</li></ul>`);
    console.log(joinedTraits)

    // Extracting cross-references:
    var xrefQueryURL = EnsemblRestBaseURL + '/xrefs/id/' + geneData.id + '?content-type=application/json'
    var xrefData = getEnsemblREST(xrefQueryURL);
    var entrezID = "NA"
    for ( xref of xrefData ){
        if ( xref.dbname == "EntrezGene" ){
            entrezID = xref.primary_id
        }
    }

    // Adding automatic cross references pointing to Ensembl:
    $("#ensembl_button").attr('onclick', "window.open('"+EnsemblURL+"Summary?db=core;g="+geneData.id+"',    '_blank')");
    $("#ensembl_phenotype_button").attr('onclick', "window.open('"+EnsemblURL+"Phenotype?db=core;g="+geneData.id+"',    '_blank')");
    $("#ensembl_pathway_button").attr('onclick', "window.open('"+EnsemblURL+"Pathway?db=core;g="+geneData.id+"',    '_blank')");
    $("#ensembl_regulation_button").attr('onclick', "window.open('"+EnsemblURL+"Regulation?db=core;g="+geneData.id+"',    '_blank')");
    $("#ensembl_expression_button").attr('onclick', "window.open('"+EnsemblURL+"ExpressionAtlas?db=core;g="+geneData.id+"',    '_blank')");

    // Adding automatic cross reference pointing to Open targets:
    $("#opentargets_button").attr('onclick', "window.open('"+OpenTargetsURL+ geneData.id+"',    '_blank')");

    // Looping through the cross references and extract entrez id:
    if ( entrezID != "NA" ){
        $("#entrez_button").attr('onclick', "window.open('"+EntrezURL+ entrezID + "',    '_blank')");
    }

    // Print out some info to make sure things are not messed up completely:
    console.log("[Info] Number of reported traits:" + reportedTraits.length)
    console.log("[Info] ID: " + geneData.id);
    console.log("[Info] Biotype: " + geneData.biotype);
    console.log("[Info] Description: " + geneData.description);
    console.log("[Info] Genomic location: " + geneData.seq_region_name + ":" + geneData.start + "-" + geneData.end)

    // OK, loading is complete:
    hideLoadingOverLay('#summary-panel-loading');
}


