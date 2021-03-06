function displayDatatableUnpublishedSummaryStats(data) {
    var filename = buildFileName('list_gwas_unpublished_summary_statistics_');

    $.each(data, (index, summary_stats) => {
        var tmp = {};
        var ftpPath = gwasProperties.FTP_PATH_PREFIX.concat(summary_stats.file);
        var ftplink = "<a href='" + ftpPath.concat("' target='_blank'>") + "Download</a>";
        summary_stats['path'] = ftplink;
        summary_stats['author'] = summary_stats['body_of_work'][0]['first_author']
        summary_stats['pubmed_id'] = summary_stats['body_of_work'][0]['pubmed_id']
        summary_stats['journal'] = summary_stats['body_of_work'][0]['journal']
        summary_stats['title'] = summary_stats['body_of_work'][0]['title']
        summary_stats['publication_date'] = summary_stats['body_of_work'][0]['publication_date']
    });
    $('#summary-stats-unpublished-table').bootstrapTable('removeAll');
    $('#summary-stats-unpublished-table').bootstrapTable({
        exportDataType: 'all',
        exportOptions: {
            fileName: filename
        },
        columns: [
            {
                field: 'author',
                title: 'First Author',
                sortable: true
            },
            {
                field: 'pubmed_id',
                title: 'PubMed ID',
                sortable: true
            },
            {
                field: 'accession',
                title: 'Study accession',
                sortable: true
            },
            {
                field: 'publication_date',
                title: 'Publication date',
                sortable: true
            },
            {
                field: 'journal',
                title: 'Journal',
                sortable: true
            },{
                field: 'title',
                title: 'Title',
                sortable: true,
                //          width:"120", //This works when the table is not nested into other tag, for example, in a simple Div
            },
            {
                field: 'trait',
                title: 'Reported trait',
                sortable: true
            },
            {
                field: 'path',
                title: 'FTP Path',
                visible: true
            }
        ],
        data: data,

    });
    $('#summary-stats-unpublished-table').bootstrapTable('load',data);
    if(data.length>5){
        $('#summary-stats-unpublished-table').bootstrapTable('refreshOptions',{showRefresh: true,pagination:true,pageSize: pageRowLimit, pageList: [5,10,25,50,100,'All']})
    }
    // Add custom tooltip text for button
    $('.keep-open').attr('title','Add/Remove Columns');
    hideLoadingOverLay('#summary-stats-table-loading');

}
function buildFileName(base) {
    // get current date
    var d = new Date();
    var curr_date = d.getDate();
    var monthNames = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
        "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"];
    var curr_month = monthNames[d.getMonth()];
    var curr_year = d.getFullYear();
    var date = (curr_date + "_" + curr_month + "_" + curr_year);

    var filename = base+date;
    return filename;
}
    /** First refactoring action: common js - DRY. From Xin original code
 *  Datatable action to rendering summary stats panel.
 *  SAB 2018. MVP.
 * */

/**
 * display summary_stats table
 * @param {Object} data - summary_stats solr docs
 * @param {Boolean} cleanBeforeInsert
 */
function displayDatatableSummaryStats(data, summaryStatsStudyAccessions) {
    //by default, we clean the table before inserting data
    var summary_stats_ids = [];
    $('#summary-stats-table').bootstrapTable('removeAll');
    
    var data_json = []
    $.each(data.response.docs, (index, summary_stats) => {
        var tmp={};
    var summary_stats_id = summary_stats.accessionId;
        summary_stats_ids.push(summary_stats_id);
        tmp['accessionId'] = '<a href="'+gwasProperties.contextPath+'studies/'+summary_stats_id+'">'+summary_stats_id+'</a>';
        tmp['author'] = summary_stats.author_s;
        tmp['pubmedId'] = summary_stats.pubmedId;
        tmp['title'] = summary_stats.title;
        tmp['journal'] = summary_stats.publication;
    // Publication date
    var p_date = summary_stats.publicationDate;
    var publi = p_date.split('T')[0];
    tmp['publication_date'] = publi;

    // var catalog_publish_date = summary_stats.catalogPublishDate;
    // var cp_date = catalog_publish_date.split('T')[0];
    // tmp['catalog_publish_date'] = cp_date;

    tmp['reported_trait'] = summary_stats.traitName_s;
    
    var mappedTraits = summary_stats.mappedLabel;
    if (mappedTraits) {
        $.each(mappedTraits, function (index, trait) {
            var link = gwasProperties.contextPath + 'efotraits/' + summary_stats.mappedUri[index].split('/').slice(-1)[0];
            mappedTraits[index] = setExternalLinkText(link, trait);
        });
        tmp['mappedTraits'] = mappedTraits.join(', ');
    } else {
        tmp['mappedTraits'] = '-';
    }
    // Number Associations
    var nr_association = 0;
    if ('association_rsId' in summary_stats) {
        var arraysize= summary_stats.association_rsId;
        nr_association = arraysize.length;
    }
    
    tmp['nr_associations'] = nr_association.toString();
    var a = (summary_stats.authorAscii_s).replace(/\s/g,"");
    var dir = a.concat("_").concat(summary_stats.pubmedId).concat("_").concat(summary_stats.accessionId);

    // Data Access
    var ftpPath = gwasProperties.FTP_PATH_PREFIX.concat(dir);

    var ftplink = "<a href='"+ftpPath.concat("' target='_blank'>");

    var linkFullPValue = ftplink.concat("FTP Download</a>");

    if ($.inArray(summary_stats_id, summaryStatsStudyAccessions) != -1) {
        var apiLink = "&nbsp;&nbsp;or&nbsp;&nbsp;<a href='http://www.ebi.ac.uk/gwas/summary-statistics/docs' target='_blank'>API access</a>";
        tmp['link'] = linkFullPValue+apiLink;
    } else {
        tmp['link'] = linkFullPValue;
    }

    // FTP Path
    tmp['ftpPath'] = ftpPath;

    data_json.push(tmp);
    
    });

    var filename = buildFileName('list_gwas_summary_statistics_');
    $('#summary-stats-table').bootstrapTable({
        exportDataType: 'all',
        exportOptions: {
            fileName: filename
        },
        columns: [
            {
                field: 'author',
                title: 'First Author',
                sortable: true
            },
            {
                field: 'pubmedId',
                title: 'PubMed ID',
                sortable: true
            },
            {
                field: 'accessionId',
                title: 'Study accession',
                sortable: true
            },
            {
                field: 'publication_date',
                title: 'Publication date',
                sortable: true
            },
            {
                field: 'journal',
                title: 'Journal',
                sortable: true
            },{
                field: 'title',
                title: 'Title',
                sortable: true,
      //          width:"120", //This works when the table is not nested into other tag, for example, in a simple Div
            },
            {
                field: 'mappedTraits',
                title: 'Trait(s)',
                sortable: true
            },
            {
                field: 'reported_trait',
                title: 'Reported trait',
                sortable: true
            },
            {
                field: 'nr_associations',
                title: 'Association count',
                sortable: true
            },
            {
                field: 'link',
                title: 'Data access',
                sortable: true
            },
            {
                field: 'ftpPath',
                title: 'FTP Path',
                visible: false
            }
        ],
        data: data_json,
        
    });
    
    $('#summary-stats-table').bootstrapTable('load',data_json);
    if(data_json.length>5){
        $('#summary-stats-table').bootstrapTable('refreshOptions',{showRefresh: true,pagination:true,pageSize: pageRowLimit, pageList: [5,10,25,50,100,'All']})
    }
    // Add custom tooltip text for button
    $('.keep-open').attr('title','Add/Remove Columns');
    hideLoadingOverLay('#summary-stats-table-loading');
}


