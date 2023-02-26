import { computeMarkerSizeOvertime } from './plotUtils.js';

var plotData = {};

function plotMeasurementOvertime1Gene(
    result,
    scaleData,
    celltypeOrder,
    refresh=false) {

    let htmlElementId = 'plotDiv';
    let htmlElement = document.getElementById(htmlElementId);

    if ($('#expressionPlot').html() === "") {
        refresh = true;
    }

    let featureName = result['feature'];
    let x_axis;
    if (celltypeOrder === "original") {
        x_axis = result['celltypes'];
    } else {
        x_axis = result['celltypes_hierarchical'];
    }
    let y_axis = result['row_labels'];

    let longestXlabel = 0, longestYlabel = 0;
    for (let i=0; i < x_axis.length; i++) {
        longestXlabel = Math.max(longestXlabel, result['celltypes'][i].length);
    }
    for (let i=0; i < y_axis.length; i++) {
        longestYlabel = Math.max(longestYlabel, result['row_labels'][i].length);
    }

    let nx = x_axis.length;
    let ny = y_axis.length;
    let pxCell = 40, pxChar = 4.2, plotGap = 10;
    let ytickMargin = 75 + pxChar * longestYlabel;
    let xtickMargin = 15 + pxChar * longestXlabel;
    let graphWidth = ytickMargin + pxCell * nx + 60;
    let graphHeight = pxCell * ny + xtickMargin;

    let title;
    if (('geneId' in result) && (result['geneId'] !== "")) {
        let geneId = result['geneId'];
        if (geneId.startsWith('MGI')) {
            geneUrl = 'http://www.informatics.jax.org/marker/'+geneId;
        } else {
            geneUrl = 'https://www.genecards.org/cgi-bin/carddisp.pl?gene='+geneId;
        }
        title = '<a href="'+geneUrl+'">'+featureName+'</a> expression over time';
    } else {
        title = featureName + ' expression over time';
    }

    let trace = {
        mode: 'markers',
        x: [],
        y: [],
        marker: {
            symbol: 'square',
            colorscale: 'Reds',
            colorbar: {},
            color: [],
            size: [],
        },
        'hoverinfo': 'text',
        text: [],
    };
    for (let i = 0; i < y_axis.length; i++) {
        const label = y_axis[i];
        for (let j = 0; j < x_axis.length; j++) {
            const celltype = x_axis[j];
            trace['x'].push(celltype)
            trace['y'].push(label)

            let ge = result['measurement'][label][celltype];
            if (scaleData == "log10") {
                ge = Math.log10(ge + 0.5);
            }
            trace['marker']['color'].push(ge);

            const ms = computeMarkerSizeOvertime(result['ncells'][label][celltype]);
            trace['marker']['size'].push(ms);

            const labelArray = label.split("_");
            const tooltip = "Expression: "+ge+", Dataset: "+labelArray[0]+", Time point: "+labelArray[1];
            trace['text'].push(tooltip);
        }
    }

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
            tickangle: 270,
            automargin: true,
            linewidth: 0,
            type: 'category',
        },
        yaxis: {
            automargin: false,
            autorange: 'reversed',
            type: 'category',
            scaleanchor: 'x',
            scaleratio: 1,
        },
    };

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
            name: 'Download expression as CSV',
            icon: Plotly.Icons.disk,
            click: function(gd) {
                var text = '';
                let geneExps = gd['data'][0]['marker']['color'];
                const nct = x_axis.length;
                // Header with cell type names
                text += 'Gene';
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
        }],
    }

    if (refresh) {
        Plotly.newPlot(htmlElementId, [trace], layout, config);
    } else {
        Plotly.react(htmlElementId, [trace], layout, config);
    }
}


function AssembleAjaxRequest(errorCallback) {
    var feature = $('#searchFeatures').val();
    let requestData = {
        feature: feature,
        species: species,
        tissue: tissue,
    }
    $.ajax({
        type:'GET',
        url:'/data/overtime_1feature',
        data: $.param(requestData),
        dataType:'json',
        success: function(result) {
            plotData = result;
            
            updateSimilarFeatures();

            updatePlot(true);
        },
        error: function (e) {
          errorCallback(e);
          alert('Request data Failed')
        }
    });
}


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

// Check another species, same gene
function onClickSpeciesSuggestions() {
    var featureName = $('#searchFeatures').val();
    const newSpecies = this.id.slice("suggest".length);
    let requestData = {
        newSpecies: newSpecies,
        feature: featureName,
        species: species,
        tissue: tissue,
    }
    $.ajax({
        type:'GET',
        url:'/data/overtime_1feature',
        data: $.param(requestData),
        success: function(result) {
            // Store global variable
            plotData = result;

            species = newSpecies;
            $("#speciesSuggestionActive").text(species.charAt(0).toUpperCase() + species.slice(1));

            updateSimilarFeatures();

            // Update search box: corrected gene names, excluding missing genes
            setSearchBox(result['feature']);

            updatePlot(true);
        },
        error: function (e) {
          alert('Error: Could not find orthologs for '+featureName+'.')
        }
    });
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


function setSearchBox(text, gseaText = "") {
    $('#searchFeatures').val(text);
}



function updatePlot(refresh=false) {
    let scaleData, celltypeOrder;
    
    if ($("#cpmTab").hasClass('is-active')) {
      scaleData = "original";
    } else {
      scaleData = "log10";
    } 
    
    if ($("#originalOrderTab").hasClass('is-active')) {
        celltypeOrder = "original";
    } else {
        celltypeOrder = "hierarchical";
    }

    plotMeasurementOvertime1Gene(
        plotData, scaleData, celltypeOrder,
        refresh);
}

$("#searchOnClick").click(AssembleAjaxRequest);
$(".speciesSuggestion").click(onClickSpeciesSuggestions);
$(".tissueSuggestion").click(onClickTissueSuggestions);

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

$("body").keyup(function(event) {
    if (event.keyCode === 13) {
        $("#searchOnClick").click();
    }
});
$(document).ready(function() {
    AssembleAjaxRequest();
});
