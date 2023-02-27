import { dotPlotSizeToFrac, dotPlotFracToSize, getDomains, getPseudocount, getTickTexts } from './plotUtils.js';

// Plot heatmap by celltype as a callback for the AJAX request
// Use global variables to store persistent data
var plotData = {};

function plotMeasurementSpeciesComparison1Feature(
    result,
    dataScale,
    tableOrder,
    heatDot) {

    let htmlElementId = 'plotDiv';
    let htmlElement = document.getElementById(htmlElementId);

    let x_axis, y_axis;
    if (tableOrder == "original") {
        y_axis = result['celltypes'];
        x_axis = result['speciess'];
    } else {
        y_axis = result['celltypes_hierarchical'];
        x_axis = result['speciess_hierarchical'];
    }

    let longestXlabel = 0, longestYlabel = 0;
    for (let i=0; i < y_axis.length; i++) {
        longestYlabel = Math.max(longestYlabel, result['celltypes'][i].length);
    }
    for (let i=0; i < x_axis.length; i++) {
        longestXlabel = Math.max(longestXlabel, result['speciess'][i].length);
    }

    let ncelltypes = y_axis.length;
    let nspecies = x_axis.length;
    let pxCell = 40, pxChar = 4.4, plotGap = 10;
    let ytickMargin = 85 + pxChar * longestYlabel;
    let xtickMargin = 15 + pxChar * longestXlabel;
    let graphWidth = ytickMargin + pxCell * ncelltypes + 60;
    let dendrographHeight = 50;
    let graphHeight = pxCell * ncelltypes + dendrographHeight + plotGap + xtickMargin;

    let yAxisDomains = [
        [0, 1.0 * (pxCell * ncelltypes) / graphHeight],
        [1.0 * (pxCell * ncelltypes + plotGap) / graphHeight, 1.0],
    ];

    // Fill trace data
    let zs = [];
    for (let i = 0; i < y_axis.length; i++) {
        const celltype = y_axis[i];
        zs.push([]);
        for (let j = 0; j < x_axis.length; j++) {
            const spec = x_axis[j];
            let measurement;
            if (dataScale == "original") {
                measurment = result['data'][celltype][spec];
            } else {
                let pseudoCount = getPseudocount(result['feature_type']);
                measurement = Math.log10(result['data'][celltype][spec] + pseudoCount);
            }
            zs[i].push(measurement);
        }
    }

    // Layout for plotly
    let layout = {
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

    // Config for plotly
    let config = {
      scrollZoom: false,
      editable: false,
      staticPlot: false,
      responsive: true,
      modeBarButtonsToRemove: ['toImage'],
      modeBarButtonsToAdd: [
        {
          name: 'Download plot as a PNG',
          icon: Plotly.Icons.camera,
          click: function(gd) {
            Plotly.downloadImage(gd, {format: 'png'})
          }
        },
        {
          name: 'Download plot as an SVG',
          icon: Plotly.Icons.camera,
          click: function(gd) {
            Plotly.downloadImage(gd, {format: 'svg'})
          }
        },
      ],
    }

    // Trace for plotly
    let trace = {
        type: 'heatmap',
        hoverongaps: false,
        z: zs,
        x: x_axis,
        y: y_axis,
    }
    if (dataScale === "log2FC") {
        trace['colorscale'] = 'RdBu';
        trace['zmid'] = 0;
    } else {
        trace['colorscale'] = 'Reds';
        trace['zmid'] = '';
    }

    if ($('#'+htmlElementId).html() === "") {
        Plotly.newPlot(htmlElement, [trace], layout, config); 
    } else {
        Plotly.react(htmlElement, [trace], layout, config); 
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
    let tableOrder = "original";
    if (!$("#originalOrderTab").hasClass('is-active')) {
        tableOrder = "hierarchical";
    }

    // NOTE: plotData is the global persistent object
    plotMeasurementSpeciesComparison1Feature(
        plotData, 
        dataScale,
        tableOrder,
    );
}


// gene of interest: Car4,Vwf,Col1a1,Ptprc,Ms4a1
// Col1a1,Fsd1l
function AssembleAjaxRequest() {
    // Get the list of genes to plot from the search box
    let feature = $('#searchGeneName').val();
    // NOTE: you cannot cache the genes because the hierarchical clustering
    // will differ anyway

    let requestData = {
        feature: feature,
        tissue: tissue,
        species: species,
    }

    // sent conditions and gene names to the API
    $.ajax({
        type:'GET',
        url:'/data/speciescomparison',
        data: $.param(requestData),
        dataType:'json',
        success: function(result) {
            plotData = result;
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
