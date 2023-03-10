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
    let graphWidth = ytickMargin + pxCell * nspecies + 60;
    let dendrographHeight = 50;
    let graphHeight = pxCell * ncelltypes + dendrographHeight + plotGap + xtickMargin;

    let yAxisDomains = [
        [0, 1.0 * (pxCell * ncelltypes) / graphHeight],
        [1.0 * (pxCell * ncelltypes + plotGap) / graphHeight, 1.0],
    ];

    // Fill trace data
    let zs = [];
    let measurement;
    for (let i = 0; i < y_axis.length; i++) {
        const celltype = y_axis[i];
        zs.push([]);
        for (let j = 0; j < x_axis.length; j++) {
            const spec = x_axis[j];
            if (dataScale == "original") {
                measurement = result['data'][spec][celltype];
            } else {
                let pseudoCount = getPseudocount(result['feature_type']);
                measurement = Math.log10(result['data'][spec][celltype] + pseudoCount);
            }
            zs[i].push(measurement);
        }
    }

    // Layout for plotly
    let layout = {
        autosize: true,
        width: graphWidth,
        height: graphHeight,
        xaxis: {
            automargin: true,
            tickangle: 90,
            scaleanchor: 'y',
            scaleratio: 1,
            type: 'category',
        },
        yaxis: {
            automargin: true,
            autorange: "reversed",
            type: 'category',
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
    if ($("#logTab").hasClass('is-active')) {
        dataScale = "log10";
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
function AssembleAjaxRequest(errorCallback) {
    // Get the list of genes to plot from the search box
    let feature = $('#searchFeatures').val();
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
        url:'/data/speciescomparison/1feature',
        data: $.param(requestData),
        dataType:'json',
        success: function(result) {
            plotData = result;

            updateSimilarFeatures();

            updatePlot();
        },
        error: function (e) {
            errorCallback(e);
            alert('Request data Failed');
        }
    });

};

function updateSimilarFeatures() {
    let similarFeatures = plotData['similar_features'];
    let featureTypes = ['gene_expression', 'chromatin_accessibility'];
    let divIds = ['gene', 'region'];

    for(let k=0; k < featureTypes.length; k++) {
        let htmlDiv = $('#'+divIds[k]+'SuggestionsDropdown');
        let suggestions = similarFeatures[featureTypes[k]];
        if (!suggestions) {
            continue;
        }
        // Empty current suggestions
        htmlDiv.html("");
        for(let i=0; i < suggestions.length; i++) {
            if (i != 0) {
                htmlDiv.append("<hr class=\"dropdown-divider\">");
            }
            htmlDiv.append("<a href=\"#\" class=\"dropdown-item featureSuggestion\" id=\"featureSuggestion_" + suggestions[i] + "\">" + suggestions[i] + "</a>");
        }
    }

    // Rebind the callback since the old elements are gone
    $(".featureSuggestion").click(onClickFeatureSuggestion);

}

// Check out another tissue
function onClickTissueSuggestions() {
    var newTissue = $(this).text().trim();
    tissue = newTissue;
    $("#tissueSuggestionActive").text(tissue);
    AssembleAjaxRequest();
}

// Check out a similar feature
function onClickFeatureSuggestion() {
    let oldFeature = $('#searchFeatures').val();
    let newFeature = $(this).text();
    $('#searchFeatures').val(newFeature);

    // Get new data and plot
    AssembleAjaxRequest(function(e) {
        $('#searchFeatures').val(oldFeature);
    });
}


////////////////////
// EVENTS
////////////////////
// Normalise the data with log10 and generate a new plot (when user click the button)
$("#log10OnClick" ).click(function() {
    // if User has input their features of interest, generate the heatmap with that value
    // otherwise use the default data
    $("#logTab").addClass('is-active');
    $("#cpmTab").removeClass('is-active');
    updatePlot();
});

$("#CPMOnClick" ).click(function() {
    $("#logTab").removeClass('is-active');
    $("#cpmTab").addClass('is-active');
    updatePlot();
});

// Second set of buttons
$("#hClusterOnClick" ).click(function() {
    // if User has input their features of interest, generate the heatmap with that value
    // otherwise use the default data
    $("#hierachicalTab").addClass('is-active');
    $("#originalOrderTab").removeClass('is-active');
    updatePlot();
});


$("#originalOnClick").click(function() {
    $("#originalOrderTab").addClass('is-active');
    $("#hierachicalTab").removeClass('is-active');
    updatePlot();
});

// Third set of buttons
$("#heatOnClick").click(function() {
    $("#heatTab").addClass('is-active');
    $("#dotTab").removeClass('is-active');
    updatePlot(true);
});

$("#dotOnClick").click(function() {
    $("#dotTab").addClass('is-active');
    $("#heatTab").removeClass('is-active');
    updatePlot(true);
});

// Suggestions
$(".tissueSuggestion").click(onClickTissueSuggestions);
$("#geneSimilar").click(function() {
    return onClickFeatureSimilarSuggestions("gene_expression");
});
$("#regionSimilar").click(function() {
    return onClickFeatureSimilarSuggestions("chromatin_accessibility");
});
$("#geneNearby").click(function() {
    return onClickFeatureNearbySuggestions("gene_expression");
});
$("#regionNearby").click(function() {
    return onClickFeatureNearbySuggestions("chromatin_accessibility");
});

// Search
$("#searchOnClick").click(function() { AssembleAjaxRequest() });
$("body").keyup(function(event) {
    if (event.keyCode === 13) {
        $("#searchOnClick").click();
    }
});

$(document).ready(function() {
    AssembleAjaxRequest();
});
