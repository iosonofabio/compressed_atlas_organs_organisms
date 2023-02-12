// Plot heatmap by celltype as a callback for the AJAX request
// Use global variables to store persistent data
var heatmapData = {};

function HeatmapDisease(result, dataScale, order) {
    if (!result) {
        alert("Error: Nothing to plot or feature names invalid")
        return;
    }

    const title = result['dataset']+", "+result['timepoint'];
    const html_element_id = "measurementPlot";

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

    let data_content = [];
    for (let i = 0; i < y_axis.length; i++) {
        const gene = y_axis[i];
        data_content.push([]);
        for (let j = 0; j < x_axis.length; j++) {
            const ct = x_axis[j];
            let geneExp = result['data'][gene][ct];
            if (dataScale == "log2FC") {
                const geneExpBaseline = result['data_baseline'][gene][ct];
                geneExp = Math.log2(geneExp + 0.5) - Math.log2(geneExpBaseline + 0.5);
            } else if (dataScale == "log10") {
                geneExp = Math.log10(geneExp + 0.5);
            }
            data_content[i].push(geneExp);
        }
    }
    var data = {
        z: data_content,
        x: x_axis,
        y: y_axis,
        type: 'heatmap',
        hoverongaps: false,
    }
    if (dataScale == "log2FC") {
        data['colorscale'] = 'RdBu';
        data['zmid'] = 0;
    } else {
        data['colorscale'] = 'Reds';
        data['zmid'] = '';
    }

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
        },
    };
        
    Plotly.newPlot(
        document.getElementById(html_element_id),
        [data],
        layout,
    ); 
} 


// gene of interest: Car4,Vwf,Col1a1,Ptprc,Ms4a1
// Col1a1,Fsd1l
function AssembleAjaxRequest() {

    // Get the list of genes to plot from the search box
    var gene_names = $('#searchFeatures').val();
    let requestData = {
        gene_names: gene_names,
        species: species,
    }

    // sent gene names to the API
    $.ajax({
        type:'GET',
        url:'/data/disease',
        data: $.param(requestData),
        dataType:'json',
        success: function(result) {
            heatmapData = result;
            console.log(result);
            updatePlot();
        },
        error: function (e) {
            console.log(e);
            alert('Request data Failed');
        }
    });

};


function updatePlot() {
    // NOTE: heatmapData is the global persistent object
    let dataScale, order;
    if ($("#log2FCTab").hasClass('is-active')) {
        dataScale = 'log2FC';
    } else if ($("#cpmTab").hasClass('is-active')) {
        dataScale = "cpm";
    } else {
        dataScale = "log10";
    }
    
    if ($("#originalOrderTab").hasClass('is-active')) {
      order = "original";
    } else {
      order = "hierachical";
    }

    // Plot each condition
    for(let i = 0; i < heatmapData.length; i++) {
        item = heatmapData[i];
        let dsTp = item['dataset']+'_'+item['timepoint'];
        console.log(dsTp);
        if ($("#dropdownItem_"+dsTp).hasClass('is-active')) {
            HeatmapDisease(
                  item, 
                  dataScale,
                  order,
                );
        }
    }
}

// Both on click and load, plot the heatmap
$("#searchOnClick").click(AssembleAjaxRequest);
$(document).ready(function() {
    AssembleAjaxRequest();

    // add event listeners for dropdown menu
    $(".datasetTimepointDropdown").click(function() {
        $(".datasetTimepointDropdown").removeClass('is-active');
        $(this).addClass('is-active');
        updatePlot();
    });
});

// normalization
$("#log2FCOnClick" ).click(function() {
    $("#log2FCTab").addClass('is-active');
    $("#cpmTab").removeClass('is-active');
    $("#logTab").removeClass('is-active');
    updatePlot();
});

$("#log10OnClick" ).click(function() {
    $("#log2FCTab").removeClass('is-active');
    $("#cpmTab").removeClass('is-active');
    $("#logTab").addClass('is-active');
    updatePlot();
});

$("#CPMOnClick" ).click(function() {
    $("#log2FCTab").removeClass('is-active');
    $("#cpmTab").addClass('is-active');
    $("#logTab").removeClass('is-active');
    updatePlot();
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

