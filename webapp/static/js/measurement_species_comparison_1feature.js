// Plot heatmap by celltype as a callback for the AJAX request
// Use global variables to store persistent data
function HeatmapSpeciesComparison(result, html_element_id, dataScale, order) {
    if (!result) {
        alert("Error: Nothing to plot or gene names invalid")
        return;
    }

    let x_axis, y_axis;
    if (order == "original") {
        x_axis = result['celltypes'];
        if ((dataScale == "originalBaseline") | (dataScale == "log10Baseline")) {
            y_axis = result['genes_baseline'];
        } else {
            y_axis = result['genes'];
        }
    } else {
        x_axis = result['celltypes_hierarchical'];
        if ((dataScale == "originalBaseline") | (dataScale == "log10Baseline")) {
            y_axis = result['genes_hierarchical_baseline'];
        } else {
            y_axis = result['genes_hierarchical'];
        }
    }
    var ngenes =  y_axis.length;
    var graph_width = 1300;
    var graph_height = 370 + 26 * ngenes;

    // Add hyperlinks to gene names
    let yticktext = [];
    for (let i = 0; i < y_axis.length; i++) {
        const gene = result['genes'][i];
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
        const geneBaseline = result['genes_baseline'][i];
        data_content.push([]);
        for (let j = 0; j < x_axis.length; j++) {
            const ct = x_axis[j];
            let geneExp;
            if (dataScale == "original") {
                geneExp = result['data'][gene][ct];
            } else if (dataScale == "log10") {
                geneExp = Math.log10(result['data'][gene][ct] + 0.5);
            } else if (dataScale == "log2FC") {
                geneExp = result['data'][gene][ct];
                let geneExpBaseline = result['data_baseline'][geneBaseline][ct];
                geneExp = Math.log2(geneExp + 0.5) - Math.log2(geneExpBaseline + 0.5);
            } else if (dataScale == "originalBaseline") {
                geneExp = result['data_baseline'][gene][ct];
            } else if (dataScale == "log10Baseline") {
                geneExp = Math.log10(result['data_baseline'][gene][ct] + 0.5);
            }
            data_content[i].push(geneExp);
        }
    }

    let title;
    if (dataScale == "log2FC") {
        title = 'Differential expression '+heatmapData['species']+' vs '+heatmapData['species_baseline'];
    } else if ((dataScale == "originalBaseline") | (dataScale == "log10Baseline")) {
        title = 'Expression in '+heatmapData['species_baseline'];
    } else {
        title = 'Expression in '+heatmapData['species'];
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
                    tickvals: y_axis,
                    ticktext: yticktext,
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
    } else if ($("#cpmBaselineTab").hasClass('is-active')) {
        dataScale = "originalBaseline";
    } else if ($("#logBaselineTab").hasClass('is-active')) {
        dataScale = "log10Baseline";
    }
    let celltypeOrder = "original";
    if (!$("#originalOrderTab").hasClass('is-active')) {
        celltypeOrder = "hierarchical";
    }

    // NOTE: heatmapData is the global persistent object
    HeatmapSpeciesComparison(
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
    let geneNames = $('#searchGeneName').val();
    // NOTE: you cannot cache the genes because the hierarchical clustering
    // will differ anyway

    let requestData = {
        genes: geneNames,
        species: species,
        species_baseline: heatmapData['species_baseline'],
    }

    // sent conditions and gene names to the API
    $.ajax({
        type:'GET',
        url:'/data/speciescomparison',
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
$("body").keyup(function(event) {
    if (event.keyCode === 13) {
        $("#searchOnClick").click();
    }
});

// On load, the heatmap data are already embedded in the template
$(document).ready(updatePlot);

// normalization
$(".dataScaleButton" ).click(function() {
    $(".dataScaleButton").removeClass('is-active');
    $(this).addClass('is-active');
    updatePlot()
});
// order of cell types
$(".dataOrderButton" ).click(function() {
    $(".dataOrderButton").removeClass('is-active');
    $(this).addClass('is-active');
    updatePlot();
});
