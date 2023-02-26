import { computeMarkerSizeOvertime } from './plotUtils.js';

var plotData = {};

function plotMeasurementOvertime1Celltype(
    result,
    scaleData,
    tableOrder) {

    let htmlDivId = 'plotDiv';
    let htmlDiv = document.getElementById(htmlDivId);

    let nPlots = result['data'].length;
    let y_axis = result['row_labels'];
    let x_axiss;
    if (tableOrder === "original") {
        x_axiss = result['features'];
    } else {
        x_axiss = [];
        for (let k=0; k < nPlots; k++) {
            x_axiss.push([]);
            for (let i = 0; i < result['features_hierarchical'][k].length; i++) {
                const feature = result['features_hierarchical'][k][i];
                x_axiss[k].push(feature);
            }
        }
    }

    let longestYlabel = 0, longestXlabel = 0;
    for (let i=0; i < y_axis.length; i++) {
        longestYlabel = Math.max(longestYlabel, result['row_labels'][i].length);
    }
    for (let k=0; k < nPlots; k++) {
        for (let i=0; i < x_axiss[k].length; i++) {
            longestXlabel = Math.max(longestXlabel, result['features'][k][i].length);
        }
    }

    let celltype = result['celltype'];
    let nfeatures = x_axiss.reduce((acc, a) => acc + a.length, 0);
    let ny = y_axis.length;
    let pxCell = 40, pxChar = 4.2, plotGap = 10;
    let ytickMargin = 75 + pxChar * longestYlabel;
    let xtickMargin = 15 + pxChar * longestXlabel;
    let graphWidth = ytickMargin + pxCell * nfeatures + plotGap * (nPlots - 1) + 60;
    let graphHeight = pxCell * ny + xtickMargin;

    let xAxisDomains = [];
    let acc = 0;
    for (let k=0; k < nPlots; k++) {
        let frac = 1.0 * result['features'][k].length / nfeatures;
        xAxisDomains.push([acc, acc + frac]);
        acc += frac;
    }

    let xs = [],
        ys = [],
        tooltipss = [],
        markersizes = [],
        markercolors = [];
    let measurement, ms;
    for (let k=0; k < nPlots; k++) {
        xs.push([]);
        ys.push([]);
        markercolors.push([]);
        markersizes.push([]);
        tooltipss.push([]);
        for (let i = 0; i < y_axis.length; i++) {
            const label = y_axis[i];
            for (let j = 0; j < x_axiss[k].length; j++) {
                const feature = x_axiss[k][j];
                xs[k].push(feature);
                ys[k].push(label);

                // this time point is missing from this feature type
                if (!(label in result['data'][k])) {
                    measurement = -1;
                    ms = 0;
                } else {
                    measurement = result['data'][k][label][feature];
                    // FIXME: change the pseudocount for ATAC
                    if (scaleData == "log10") {
                        measurement = Math.log10(measurement + 0.5);
                    }
                    ms = computeMarkerSizeOvertime(result['ncells'][k][label][celltype]);
                }
                markercolors[k].push(measurement);
                markersizes[k].push(ms);

                const labelArray = label.split("_");
                const tooltip = "Measurement: "+measurement+", Dataset: "+labelArray[1]+", Time point: "+labelArray[2];
                tooltipss[k].push(tooltip);
            }
        }
    }

    let layout = {
        grid: {
            rows: 1,
            columns: nPlots,
        },
        margin: {
            l: ytickMargin,
            r: 0,
            b: 0,
            t: 0,
            pad: 4,
        },
        autosize: true,
        width: graphWidth,
        height: graphHeight,
        yaxis: {
            automargin: true,
            autorange: 'reversed',
            type: 'category',
            tickvals: result['yticks'],
            ticktext: result['yticktext'],
        },
    };

    let traces = [];
    for (let k=0; k < nPlots; k++) {
        let xaxisName = 'xaxis', xaxisShort = 'x';
        if (k != 0) {
            xaxisName += (k+1);
            xaxisShort += (k+1);
        }
        layout[xaxisName] = {
            tickangle: 270,
            automargin: true,
            linewidth: 0,
            type: 'category',
            scaleanchor: 'y',
            scaleratio: 1,
            domain: xAxisDomains[k],
        };

        traces.push({
            mode: 'markers',
            x: xs[k],
            y: ys[k],
            text: tooltipss[k],
            marker: {
                symbol: 'square',
                colorscale: 'Reds',
                colorbar: {},
                color: markercolors[k],
                size: markersizes[k],
            },
            'hoverinfo': 'text',
            xaxis: xaxisShort,
        });
    }

    let config = {
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
          {
            name: 'Download data as CSV',
            icon: Plotly.Icons.disk,
            click: function(gd) {
                var text = '';
                let geneExps = gd['data'][0]['marker']['color'];
                const nct = x_axis.length;
                // Header with cell type names
                text += 'Feature';
                for(var i = 0; i < nct; i++){
                    text += ',' + gd['data'][0]['x'][i];
                };
                // Gene expression
                for (var i = 0; i < geneExps.length; i++) {
                    if (i % nct == 0) {
                        text += '\n' + gd['data'][0]['y'][i];
                    }
                    text += ',' + geneExps[i];
                }
                text += '\n';

                var blob = new Blob([text], {type: 'text/plain'});
                var a = document.createElement('a');
                const object_URL = URL.createObjectURL(blob);
                a.href = object_URL;
                a.download = 'measurement_table.csv';
                document.body.appendChild(a);
                a.click();
                URL.revokeObjectURL(object_URL);
            },
          }
        ], 
    }

    // Make new plot if none is present
    if ($('#measurementPlot').html() === "") {
        Plotly.newPlot(htmlDiv, traces, layout, config); 

    // Update existing plot if present
    } else {
        Plotly.react(htmlDiv, traces, layout, config); 
    }
}


function AssembleAjaxRequest( featurestring = "") {
    let featureNames;
    if (featurestring !== "") {
        featureNames = featurestring;
    } else {
        featureNames = $('#searchFeatures').val();
    }

    let requestData = {
        celltype: celltype,
        feature_names: featureNames,
        species: species,
        tissue: tissue,
    }
    $.ajax({
        type:'GET',
        url:'/data/overtime_1celltype',
        data: $.param(requestData),
        dataType:'json',
        success: function(result) {
            plotData = result;
            updatePlot();
        },
        error: function (e) {
          alert('Request data Failed')
        }
    });
}

// Check another species, same gene
function onClickSpeciesSuggestions() {
    const genestring = $('#searchFeatures').val();
    const newSpecies = this.id.slice("suggest".length);
    let requestData = {
        celltype: celltype,
        newSpecies: newSpecies,
        gene_names: genestring,
        species: species,
    }
    $.ajax({
        type:'GET',
        url:'/data/overtime_1celltype',
        data: $.param(requestData),
        success: function(result) {
            // Store global variable
            plotData = result;

            $("#suggest"+newSpecies).text(species.slice(0, 1).toUpperCase()+species.slice(1)).prop('id', "suggest"+species);
            species = newSpecies;

            // Update search box: corrected gene names, excluding missing features
            setSearchBox(result['gene']);

            updatePlot();
        },
        error: function (e) {
          alert('Error: Could not find orthologs for '+geneName+'.')
        }
    });
}


// Check out a similar gene
function onClickFeatureSuggestions(correlatesType) {
    var featureNames = $('#searchFeatures').val();

    let requestData = {
        feature_names: featureNames,
        species: species,
        correlates_type: correlatesType,
    }
    $.ajax({
        type:'GET',
        url:'/data/features_correlated',
        data: $.param(requestData),
        success: function(result) {
            // Update search box: corrected feature names, excluding missing ones
            setSearchBox(result);

            // Request data
            AssembleAjaxRequest();
        },
        error: function (e) {
          alert('Error: Could not find genes correlated with ' + featureNames + '.')
        }
    });

}


// Check out another cell type
function onClickCelltypeSuggestions() {
    var newCelltype = $(this).text().trim();
    celltype = newCelltype;
    $("#celltypeSuggestionActive").text(celltype);
    AssembleAjaxRequest();
}


function setSearchBox(text, gseaText = "") {
    $('#searchFeatures').val(text);
}


function updatePlot() {
    let scaleData, tableOrder;
    
    if ($("#cpmTab").hasClass('is-active')) {
      scaleData = "original";
    } else {
      scaleData = "log10";
    } 
    
    if ($("#originalOrderTab").hasClass('is-active')) {
        tableOrder = "original";
    } else {
        tableOrder = "hierarchical";
    }

    plotMeasurementOvertime1Celltype(plotData, scaleData, tableOrder);
}

// Normalise the data with log10 and generate a new plot (when user click the button)
$("#log10OnClick" ).click(function() {
    // if User has input their gene of interest, generate the heatmap with that value
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
    // if User has input their gene of interest, generate the heatmap with that value
    // otherwise use the default data
    $("#hierachicalTab").addClass('is-active');
    $("#originalOrderTab").removeClass('is-active');
    updatePlot();
});


$("#originalOnClick" ).click(function() {
    $("#originalOrderTab").addClass('is-active');
    $("#hierachicalTab").removeClass('is-active');
    updatePlot();
});

// Suggestions
$(".speciesSuggestion").click(onClickSpeciesSuggestions);
$(".celltypeSuggestion").click(onClickCelltypeSuggestions);
$("#geneSuggestions").click(function() {
    return onClickFeatureSuggestions("gene_expression");
});
$("#regionSuggestions").click(function() {
    return onClickFeatureSuggestions("chromatin_accessibility");
});

// Search
$("#searchOnClick").click(function() { AssembleAjaxRequest() });
$("body").keyup(function(event) {
    if (event.keyCode === 13) {
        $("#searchOnClick").click();
    }
});

$(document).ready(function() {
    AssembleAjaxRequest()
});
