// Plot heatmap by celltype as a callback for the AJAX request
// Use global variables to store persistent data
function HeatmapDifferential(result, html_element_id, dataScale, order) {
    if (!result) {
        alert("Error: Nothing to plot or gene names invalid")
        return;
    }

    let x_axis, y_axis;
    if (order == "original") {
        x_axis = result['celltypes'];
        y_axis = result['genes'];
    } else {
        x_axis = result['celltypes_hierarchical'];
        y_axis = result['genes_hierarchical'];
    }
    var ngenes =  y_axis.length;
    var graph_width = 1300;
    var graph_height = 370 + 26 * ngenes;

    // Add hyperlinks to gene names
    let yticktext = [];
    for (let i = 0; i < y_axis.length; i++) {
        const gene = y_axis[i];
        const geneId = result['gene_ids'][gene];
        if (geneId === "") {
            yticktext.push(gene);
        } else {
            let geneUrl = gene;
            if (geneId.startsWith('MGI')) {
                geneUrl = 'http://www.informatics.jax.org/marker/'+geneId;
            } else {
                geneUrl = 'https://www.genecards.org/cgi-bin/carddisp.pl?gene='+geneId;
            }
            const tickText = '<a href="'+geneUrl+'">'+gene+'</a>'
            yticktext.push(tickText);
        }
    }

    // Fill heatmap data
    let data_content = [];
    for (let i = 0; i < y_axis.length; i++) {
        const gene = y_axis[i];
        data_content.push([]);
        for (let j = 0; j < x_axis.length; j++) {
            const ct = x_axis[j];
            let geneExp = result['data'][gene][ct];
            if (dataScale == "log10") {
                geneExp = Math.log10(geneExp + 0.5);
            } else if (dataScale == "log2FC") {
                // If the comparison is against a certain cell type, we subtract
                // that subtype from all others
                let geneExpBaseline;
                if (heatmapData['comparison'] == 'celltypes') {
                    let ct2 = result['celltype_baseline'];
                    geneExpBaseline = result['data_baseline'][gene][ct2];
                } else {
                    geneExpBaseline = result['data_baseline'][gene][ct];
                }
                geneExp = Math.log2(geneExp + 0.5) - Math.log2(geneExpBaseline + 0.5);
            }
            data_content[i].push(geneExp);
        }
    }

    let title;
    if (dataScale == "log2FC") {
        title = 'Differential expression ';
        if (heatmapData['comparison'] == 'hyperoxia/normal') {
            title += 'hyperoxia vs normal at '+heatmapData['timepoint'];
        } else if (heatmapData['comparison'] == 'timepoints') {
            title += heatmapData['timepoint']+' vs '+heatmapData['timepoint_baseline'];
        } else if (heatmapData['comparison'] == 'celltypes') {
            title += ' at '+heatmapData['timepoint']+' in '+heatmapData['celltype']+' vs '+heatmapData['celltype_baseline'];
        }
    }

    var data = {
            type: 'heatmap',
            hoverongaps: false,
    }
    if (dataScale === "log2FC") {
        data['colorscale'] = 'RdBu';
        data['zmid'] = 0;
    } else {
        data['colorscale'] = 'Reds';
        data['zmid'] = '';
    }

    // Make new plot if none is present
    if ($('#'+html_element_id).html() === "") {
        data['z'] = data_content;
        data['x'] = x_axis;
        data['y'] = y_axis;

        var layout = {
            autosize: true,
            width: graph_width,
            height: graph_height,
            title: title,
            xaxis: {
                //title: 'Cell types',
                automargin: true,
                tickangle: 70,
                scaleanchor: 'y',
                scaleratio: 1,
                type: 'category',
            },
            yaxis: {
                //title: 'Genes',
                automargin: true,
                autorange: "reversed",
                type: 'category',
                tickvals: y_axis,
                ticktext: yticktext,
            },
        };
            
        Plotly.newPlot(
            document.getElementById(html_element_id),
            [data],
            layout,
        ); 

    // Update existing plot if present
    } else {
        data['z'] = [data_content];
        data['x'] = [x_axis];
        data['y'] = [y_axis];
        Plotly.update(
            document.getElementById(html_element_id),
            data,
            {
                height: graph_height,
                title: title,
                yaxis: {
                    autorange: "reversed",
                },
            },
            [0],
        ); 
    }
} 


// NOTE: this is why react was invented...
function updatePlot() {
    let dataScale = "original";
    if ($("#log2FCTab").hasClass('is-active')) {
        dataScale = "log2FC";
    } else if ($("#logTab").hasClass('is-active')) {
        dataScale = "log10";
    }
    let celltypeOrder = "original";
    if (!$("#originalOrderTab").hasClass('is-active')) {
        celltypeOrder = "hierarchical";
    }

    // NOTE: heatmapData is the global persistent object
    HeatmapDifferential(
        heatmapData, 
        "h5_data_plot",
        dataScale,
        celltypeOrder,
    );
}


// gene of interest: Car4,Vwf,Col1a1,Ptprc,Ms4a1
// Col1a1,Fsd1l
function AssembleAjaxRequest() {

    // Get the list of genes to plot from the search box
    let geneNames = $('#searchFeatures').val();
    // NOTE: you cannot cache the genes because the hierarchical clustering
    // will differ anyway

    let requestData = {
        comparison: heatmapData['comparison'],
        ct1: heatmapData['celltype'],
        ct2: heatmapData['celltype_baseline'],
        ds1: heatmapData['dataset'],
        ds2: heatmapData['dataset_baseline'],
        tp1: heatmapData['timepoint'],
        tp2: heatmapData['timepoint_baseline'],
        dis1: heatmapData['disease'],
        dis2: heatmapData['disease_baseline'],
        genestr': geneNames,
        species: species,
    }

    // sent conditions and gene names to the API
    $.ajax({
        type:'GET',
        url:'/data/differential',
        data: $.param(requestData),
        dataType:'json',
        success: function(result) {
            // Clear mobile DOM elements
            $("#h5_data_plot").html("");

            heatmapData = result;
            updatePlot();

        },
        error: function (e) {
            console.log(e);
            alert('Request data Failed');
        }
    });

};


// On search click, keep the same conditions but change the genes
$("#searchOnClick").click(function() {
  // action here when clicking the search button
  AssembleAjaxRequest();
});

// On load, the heatmap data are already embedded in the template
$(document).ready(updatePlot);


// normalization
$("#log2FCOnClick" ).click(function() {
    $("#log2FCTab").addClass('is-active');
    $("#cpmTab").removeClass('is-active');
    $("#logTab").removeClass('is-active');
    updatePlot()
});

$("#log10OnClick" ).click(function() {
    $("#log2FCTab").removeClass('is-active');
    $("#cpmTab").removeClass('is-active');
    $("#logTab").addClass('is-active');
    updatePlot()
});

$("#CPMOnClick" ).click(function() {
    $("#log2FCTab").removeClass('is-active');
    $("#cpmTab").addClass('is-active');
    $("#logTab").removeClass('is-active');
    updatePlot()
});


// order of cell types
$("#hClusterOnClick" ).click(function() {
    $("#hierachicalTab").addClass('is-active');
    $("#originalOrderTab").removeClass('is-active');
    updatePlot();
});


$("#originalOnClick" ).click(function() {
    $("#originalOrderTab").addClass('is-active');
    $("#hierachicalTab").removeClass('is-active');
    updatePlot();
});

