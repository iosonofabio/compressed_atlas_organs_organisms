import { dotPlotSizeToFrac, dotPlotFracToSize, getDomains, getPseudocount, getTickTexts } from './plotUtils.js';

// Plot heatmap by celltype as a callback for the AJAX request
// Use global variables to store persistent data
var plotData = {};

function plotMeasurementSpeciesComparison1Celltype(
    result,
    dataScale,
    tableOrder,
    heatDot) {

    let htmlElementId = 'plotDiv';
    let htmlElement = document.getElementById(htmlElementId);

    let nPlots = result['data'].length;
    let x_axis, y_axiss;
    if (tableOrder == "original") {
        y_axiss = result['features'];
        x_axis = result['speciess'];
    } else {
        y_axiss = result['features_hierarchical'];
        x_axis = result['speciess_hierarchical'];
    }

    let longestXlabel = 0, longestYlabel = 0;
    for (let i=0; i < x_axis.length; i++) {
        longestXlabel = Math.max(longestXlabel, result['speciess'][i].length);
    }
    for (let k=0; k < y_axiss.length; k++) {
        for (let i=0; i < y_axiss[k].length; i++) {
            longestYlabel = Math.max(longestYlabel, result['features'][k][i].length);
        }
    }

    let nfeatures = y_axiss.reduce((acc, a) => acc + a.length, 0);
    let nspecies = x_axis.length;
    let pxCell = 40, pxChar = 4.4, plotGap = 10;
    let ytickMargin = 85 + pxChar * longestYlabel;
    let xtickMargin = 15 + pxChar * longestXlabel;
    let graphWidth = ytickMargin + pxCell * nspecies + 60;
    let dendrographHeight = 50;
    let graphHeight = pxCell * nfeatures + dendrographHeight + (plotGap * nPlots) + xtickMargin;

    let yAxisDomains = getDomains(
        y_axiss, false,
        0,
        1.0 - 1.0 * dendrographHeight / graphHeight);

    // Add hyperlinks to feature names if they are genes
    let yticktexts = [];
    for (let k=0; k < nPlots; k++) {
        let yticktexts_k = getTickTexts(
            y_axiss[k],
            result['feature_type'][k],
            result['gene_ids'][k],
        );
        yticktexts.push(yticktexts_k);
    }

    // Layout for plotly
    let traces = [];
    let layout = {
        autosize: true,
        width: graphWidth,
        height: graphHeight,
        margin: {
            l: ytickMargin,
            r: 0,
            b: 0,
            t: 0,
            pad: 4,
        },
        xaxis: {
            automargin: true,
            tickangle: 270,
            type: 'category',
        },
    };
    for (let k=0; k < nPlots; k++) {
        traces.push({});
        let yaxisName = 'yaxis', yaxisShort = 'y';
        if (k != 0) {
            yaxisName += (k+1);
            yaxisShort += (k+1);
        }
        traces[k]['yaxis'] = yaxisShort;
        layout[yaxisName] = {
            autorange: "reversed",
            type: 'category',
            automargin: false,
            scaleanchor: 'x',
            scaleratio: 1,
            tickvals: y_axiss[k],
            ticktext: yticktexts[k],
            domain: yAxisDomains[k],
        };
    }

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

    // Fill trace data
    for (let k=0; k < nPlots; k++) {
        traces[k]['type'] = 'heatmap';
        traces[k]['hoverongaps'] = false;
        traces[k]['colorscale'] = 'Reds';

        let zs = [];
        let measurement;
        for (let i = 0; i < y_axiss[k].length; i++) {
            const feature = y_axiss[k][i];
            zs.push([]);
            for (let j = 0; j < x_axis.length; j++) {
                const spec = x_axis[j];
                measurement = result['data'][k][spec][feature];
                if (dataScale != "original") {
                    let pseudoCount = getPseudocount(result['feature_type']);
                    measurement = Math.log10(measurement + pseudoCount);
                }
                zs[i].push(measurement);
            }
        }
        traces[k]['z'] = zs;
        traces[k]['x'] = x_axis;
        traces[k]['y'] = y_axiss[k];
    }

    if ($('#'+htmlElementId).html() === "") {
        Plotly.newPlot(htmlElement, traces, layout, config); 
    } else {
        Plotly.react(htmlElement, traces, layout, config); 
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
    plotMeasurementSpeciesComparison1Celltype(
        plotData, 
        dataScale,
        tableOrder,
    );
}


// gene of interest: Car4,Vwf,Col1a1,Ptprc,Ms4a1
// Col1a1,Fsd1l
function AssembleAjaxRequest() {
    // Get celltype
    let celltype = $("#celltypeSuggestionActive").text();

    // Get the list of genes to plot from the search box
    let feature_names = $('#searchFeatures').val();
    // NOTE: you cannot cache the genes because the hierarchical clustering
    // will differ anyway

    let requestData = {
        feature_names: feature_names,
        tissue: tissue,
        species: species,
        celltype: celltype,
    }

    // sent conditions and gene names to the API
    $.ajax({
        type:'GET',
        url:'/data/speciescomparison/1celltype',
        data: $.param(requestData),
        dataType:'json',
        success: function(result) {
            plotData = result;

            updateCelltypeSuggestions();
            $('#searchFeatures').val(plotData['searchstring']);

            updatePlot();
        },
        error: function(e) {
            console.log(e);
            alert('Request data Failed');
        }
    });

};


function updateCelltypeSuggestions() {

    let htmlDiv = $('#celltypeSuggestionsDropdown');
    let suggestions = plotData['celltypes_tissue'];
    if (!suggestions) {
        return;
    }
    // Empty current suggestions
    htmlDiv.html("");
    for(let i=0; i < suggestions.length; i++) {
        if (i != 0) {
            htmlDiv.append("<hr class=\"dropdown-divider\">");
        }
        htmlDiv.append("<a href=\"#\" class=\"dropdown-item celltypeSuggestion\" id=\"celltypeSuggestion_" + suggestions[i] + "\">" + suggestions[i] + "</a>");
    }

    $("#celltypeSuggestionActive").text(plotData['celltype']);

    // Rebind the callback since the old elements are gone
    $(".celltypeSuggestion").click(onClickCelltypeSuggestions);

}


// Check out another tissue
function onClickTissueSuggestions() {
    let newTissue = $(this).text().trim();
    tissue = newTissue;
    $("#tissueSuggestionActive").text(tissue);
    AssembleAjaxRequest();
}

// Check out another cell type
function onClickCelltypeSuggestions() {
    let newCelltype = $(this).text().trim();
    $("#celltypeSuggestionActive").text(newCelltype);
    AssembleAjaxRequest();
}

// Check out a similar feature
function onClickFeatureSuggestion() {
    let oldFeature = $('#searchFeatures').val();
    let newFeature = $(this).text();
    $('#searchFeatures').val(newFeature);

    // Get new data and plot
    AssembleAjaxRequest();
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
$(".celltypeSuggestion").click(onClickCelltypeSuggestions);
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
