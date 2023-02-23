import { dotPlotSizeToFrac, dotPlotFracToSize, getDomains, getPseudocount, getTickTexts } from './plotUtils.js';


// Plot heatmap by celltype as a callback for the AJAX request
// Use global variables to store persistent data
var plotData = {};

function plotMeasurementByCelltype(
    result,
    dataScale,
    tableOrder,
    heatDot,
    refresh=false) {

    let htmlElementId = 'plotDiv';
    let htmlElement = document.getElementById(htmlElementId);

    // Make new plot if none is present
    if (($('#'+htmlElementId).html() === ""))
        refresh = true;

    let nPlots = result['data'].length;
    let x_axis, y_axiss;
    if (tableOrder == "original") {
        x_axis = result['celltypes'];
        y_axiss = result['features'];
    } else {
        x_axis = [];
        for (let i = 0; i < result['celltypes_hierarchical'].length; i++) {
            const ct = result['celltypes'][result['celltypes_hierarchical'][i]];
            x_axis.push(ct);
        }
        y_axiss = [];
        for (let k=0; k < nPlots; k++) {
            y_axiss.push([]);
            for (let i = 0; i < result['features_hierarchical'][k].length; i++) {
                const feature = result['features'][k][result['features_hierarchical'][k][i]];
                y_axiss[k].push(feature);
            }
        }
    }

    let longestXlabel = 0, longestYlabel = 0;
    for (let i=0; i < x_axis.length; i++) {
        longestXlabel = Math.max(longestXlabel, result['celltypes'][i].length);
    }
    for (let k=0; k < nPlots; k++) {
        for (let i=0; i < y_axiss[k].length; i++) {
            longestYlabel = Math.max(longestYlabel, result['features'][k][i].length);
        }
    }

    let nfeatures = y_axiss.reduce((acc, a) => acc + a.length, 0);
    let ncelltypes = x_axis.length;
    let pxCell = 40, pxChar = 4.4, plotGap = 10;
    let ytickMargin = 85 + pxChar * longestYlabel;
    let xtickMargin = 15 + pxChar * longestXlabel;
    let graphWidth = ytickMargin + pxCell * ncelltypes + 60;
    let graphHeight = pxCell * nfeatures + plotGap * (nPlots - 1) + xtickMargin;

    // Height ratios for the plots
    let yAxisDomains = getDomains(y_axiss, true);

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
        grid: {
            rows: nPlots, columns: 1,
            roworder: "top to bottom",
        },
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
            autorange: true,
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

    // Heatmap
    if (heatDot == "heat") {
        // Layout options
        for (let k=0; k < nPlots; k++) {
            traces[k]['type'] = 'heatmap';
            traces[k]['hoverongaps'] = false;
            traces[k]['colorscale'] = 'Reds';
        }

        // Config options: download as CSV
        config["modeBarButtonsToAdd"].push({
            name: 'Download expression as CSV',
            icon: Plotly.Icons.disk,
            click: function(gd) {
                let text = 'Feature,' + gd['data'][0]['x'] + '\n';
                for (let k=0; k < nPlots; k++) {
                    let data = gd['data'][k]['z'];
                    for(let i=0; i < zs[k].length; i++){
                        text += gd['data'][0]['y'][i] + ',' + zs[k][i] + '\n';
                    }
                }
                var blob = new Blob([text], {type: 'text/plain'});
                var a = document.createElement('a');
                const object_URL = URL.createObjectURL(blob);
                a.href = object_URL;
                a.download = result['feature_type'] + '.csv';
                document.body.appendChild(a);
                a.click();
                URL.revokeObjectURL(object_URL);
            },
        });

        // Fill trace data
        let zs = [];
        for (let k=0; k < nPlots; k++) {
            zs.push([]);
            for (let i = 0; i < y_axiss[k].length; i++) {
                zs[k].push([]);
                for (let j = 0; j < x_axis.length; j++) {
                    let ii = i, jj = j;
                    if (tableOrder != "original") {
                        ii = result['features_hierarchical'][k][i];
                        jj = result['celltypes_hierarchical'][j];
                    }
                    let measurement = result['data'][k][ii][jj];
                    if (dataScale == "log10") {
                        let pseudoCount = getPseudocount(result['feature_type'][k]);
                        measurement = Math.log10(measurement + pseudoCount);
                    }
                    zs[k][i].push(measurement);
                }
            }
           traces[k]['z'] = zs[k];
           traces[k]['x'] = x_axis;
           traces[k]['y'] = y_axiss[k];
        }

        if (refresh == true) {
            Plotly.newPlot(htmlElement, traces, layout, config);
        } else {
            Plotly.react(htmlElement, traces, layout, config);
        }

    // Dot plot
    } else {
        // Config: modebar buttons
        config["modeBarButtonsToAdd"].push({
            name: 'Download average as CSV',
            icon: Plotly.Icons.disk,
            click: function(gd) {
                const nct = x_axis.length;
                // Header with cell type names
                let text = 'Feature';
                for(var i = 0; i < nct; i++){
                    text += ',' + gd['data'][0]['x'][i];
                };
                // Table data
                for (let k=0; k < nPlots; k++) {
                    let measurements = gd['data'][k]['marker']['color'];
                    for (var i = 0; i < measurements.length; i++) {
                        if (i % nct == 0) {
                            text += '\n' + gd['data'][k]['y'][i];
                        }
                        text += ',' + measurements[i];
                    }
                }
                text += '\n';

                var blob = new Blob([text], {type: 'text/plain'});
                var a = document.createElement('a');
                const object_URL = URL.createObjectURL(blob);
                a.href = object_URL;
                a.download = 'gene_expression.csv';
                document.body.appendChild(a);
                a.click();
                URL.revokeObjectURL(object_URL);
            },
        });
        config["modeBarButtonsToAdd"].push({
            name: 'Download fraction of cells as CSV',
            icon: Plotly.Icons.disk,
            click: function(gd) {
                const nct = x_axis.length;
                // Header with cell type names
                let text = 'Feature';
                for(let i=0; i < nct; i++){
                    text += ',' + gd['data'][0]['x'][i];
                };
                // Table data
                for (let k=0; k < nPlots; k++) {
                    let markerSizes = gd['data'][k]['marker']['size'];
                    for (let i=0; i < markerSizes.length; i++) {
                        if (i % nct == 0) {
                            text += '\n' + gd['data'][k]['y'][i];
                        }
                        let frac = dotPlotSizeToFrac(markerSizes[i], result['feature_types'][k]);
                        text += ',' + frac;
                    }
                }
                text += '\n';

                var blob = new Blob([text], {type: 'text/plain'});
                var a = document.createElement('a');
                const object_URL = URL.createObjectURL(blob);
                a.href = object_URL;
                a.download = 'fraction_of_cells.csv';
                document.body.appendChild(a);
                a.click();
                URL.revokeObjectURL(object_URL);
            },
        });

        // Fill trace data: x, y, markercolor/size
        for (let k=0; k < nPlots; k++) {
            traces[k]['mode'] = 'markers';
            traces[k]['x'] = [];
            traces[k]['y'] = [];
            traces[k]['text'] = [];
            traces[k]['marker'] = {
                symbol: 'circle',
                colorscale: 'Reds',
                colorbar: {},
                color: [],
                size: [],
            }
            traces[k]['hoverinfo'] = 'text';

            for (let i = 0; i < y_axiss[k].length; i++) {
                for (let j = 0; j < x_axis.length; j++) {
                    let ii = i, jj = j;
                    if (tableOrder != "original") {
                        ii = result['features_hierarchical'][k][i];
                        jj = result['celltypes_hierarchical'][j];
                    }
                    traces[k]['x'].push(x_axis[jj]);
                    traces[k]['y'].push(y_axiss[k][ii]);

                    let measurement = result['data'][k][ii][jj];
                    if (dataScale == "log10") {
                        let pseudocount = getPseudocount(result['feature_type'][k]);
                        measurement = Math.log10(measurement + pseudocount);
                    }
                    traces[k]['marker']['color'].push(measurement);

                    // ATAC-Seq has generally smaller fractions, so highlight more
                    let frac, markersize;
                    if (result['feature_type'][k] == 'gene_expression') {
                        frac = result['data_fractions'][k][ii][jj];
                        markersize = 2 + 18 * Math.sqrt(frac);
                    } else {
                        frac = result['data'][k][ii][jj];
                        markersize = 2 + 60 * Math.sqrt(frac);
                    }
                    traces[k]['marker']['size'].push(markersize);

                    const tooltip = "Average: "+measurement+", Fraction of cells: "+parseInt(100 * frac)+"%";
                    traces[k]['text'].push(tooltip);
                }
            }
        }

        if (refresh) {
            Plotly.newPlot(htmlElement, traces, layout, config);
        } else {
            Plotly.react(htmlElement, traces, layout, config);
        }
    }

    connectTooltip(
        result['feature_type'],
        result['GO_terms'],
        result['feature_coords']);

}


function connectTooltip(feature_types, goTerms, feature_coords) {
    let iGeneExpression = -1;
    for (let k=0; k < feature_types.length; k++) {
        if (feature_types[k] == "gene_expression") {
            iGeneExpression = k;
        }
    }
    if (iGeneExpression == -1) {
        return;
    }

    let ytickName;
    if (iGeneExpression == 0) {
        ytickName = 'ytick';
    } else {
        ytickName = 'y' + (iGeneExpression + 1) + 'tick';
    }

    // Add tooltips to gene names
    $("." + ytickName + " > text > tspan").click(function(evt) {
        // If already active, it's a second click
        let wasActive = this.classList.contains("is-active");

        // Deactivate all, and if this was not active activate it
        hideTooltip();
        $("." + ytickName + " > text > tspan").map(function(elem) {
            this.classList.remove("is-active");
            let gene = this.parentElement.getElementsByTagName("a")[0].innerHTML;
        })
        if (wasActive) {
            return;
        }

        // Activate and show pop-up window
        this.classList.add("is-active");
        let gene = this.parentElement.getElementsByTagName("a")[0].innerHTML;
        let goTermsGene = goTerms[iGeneExpression][gene];
        if (goTermsGene === undefined) {
            return;
        }
        let text = "<b>Name:</b> " + gene + "</br>";
        const geneCoordsParts = feature_coords[iGeneExpression][gene].split("-");
        const chrom = geneCoordsParts[0];
        const start = geneCoordsParts[1];
        const end = geneCoordsParts[2];
        const featureUrl = 'https://genome.ucsc.edu/cgi-bin/hgTracks?db=mm10&position=' + chrom + ':' + start + '-' + end;
        text += '<b>Position:</b> <a href="' + featureUrl + '" target="_blank">' + geneCoordsParts.join("-") + "</a></br>";
        text += '<b>Protein atlas:</b> <a href="https://www.proteinatlas.org/search/' + gene + '" target="_blank">' + gene + "</a></br>";
        text += "<b>GO terms:</b></br>";
        for (let i = 0; i < goTermsGene.length; i++) {
            text += '<div><a class="goHyperlink">'+goTermsGene[i]+"</a></div>";
        }
        showTooltip(evt, text);
    });
}

function showTooltip(evt, text) {
    let tooltip = document.getElementById("tooltip");
    tooltip.innerHTML = text;
    tooltip.style.background = "white";
    tooltip.style.border = "1px solid black";
    tooltip.style.borderRadius = "5px";
    tooltip.style.padding = "5px";
    tooltip.style.left = evt.pageX + 25 + 'px';
    tooltip.style.top = evt.pageY - 20 + 'px';
    tooltip.style.display = "block";

    $(".goHyperlink").click(onClickGOTermSuggestions);
}

function hideTooltip() {
    var tooltip = document.getElementById("tooltip");
    tooltip.style.display = "none";
}


// NOTE: this is why react was invented...
function updatePlot(refresh=false) {
    let dataScale = "original";
    if (!$("#cpmTab").hasClass('is-active')) {
        dataScale = "log10";
    }
    let tableOrder = "original";
    if (!$("#originalOrderTab").hasClass('is-active')) {
        tableOrder = "hierarchical";
    }
    let heatDot = "heat";
    if (!$("#heatTab").hasClass('is-active')) {
        heatDot = "dot";
    }

    // NOTE: plotData is the global persistent object
    if (!plotData['result']) {
        alert("Error: no data to plot");
        return;
    }

    // Plot inside plotDiv
    plotMeasurementByCelltype(
        plotData['result'],
        dataScale,
        tableOrder,
        heatDot,
        refresh,
    );

    if (refresh) {
        updatePathwayAnalysis();
    }
}

function updatePathwayAnalysis() {
    // Update pathway analysis
    let genes;
    let nPlots = plotData['result']['feature_type'].length;
    let found = false;
    for (let k=0; k < nPlots; k++) {
        if (plotData['result']['feature_type'][k] == 'gene_expression') {
            genes = plotData['result']['features'][k].join(',');
            $('#suggestGO > a').attr(
                'href', '/barplot_gsea?species='+species+'&gene_set=GO&genes='+genes);
            $('#suggestKEGG > a').attr(
                'href', '/barplot_gsea?species='+species+'&gene_set=KEGG&genes='+genes);
            found = true;
            break;
        }
    }
    if (found) {
        document.getElementById("pathwaySuggestions").style.display = "inline";
    } else {
        document.getElementById("pathwaySuggestions").style.display = "none";
    }
}

function AssembleAjaxRequest( featurestring = "" ) {
    // Get the list of features to plot from the search box
    // If too long for the search box, it should be supplied as a specialGenestring
    let featureNames;
    if (featurestring !== "") {
        featureNames = featurestring;
    } else {
        featureNames = $('#searchFeatures').val();
    }

    let requestData = {
        feature_names: featureNames,
        species: species,
    }
    // HTML GET method length is capped at 4,000, but search box might be shorter
    let htmlVerb = (featureNames.length > 500) ? 'POST' : 'GET';

    // sent feature names to the API
    // FIXME: this fails at random times with large payloads?
    $.ajax({
        type: htmlVerb,
        url:'/data/by_celltype',
        data: $.param(requestData),
        success: function(result) {
            // Store global variable
            plotData = {
                'result': result,
            };

            // Update search box: corrected feature names, excluding missing features
            setSearchBox(result['features']);

            // Create plot. Always refresh when you are asking for new data, specifically
            // to address a bug with plot domains (e.g. chromatin accessibility AND gene
            // expression)
            updatePlot(true);
        },
        error: function (e) {
            console.log(e);
            alert('Error: Could not find some feature names.')
        },
    });
};

// Check another species, same features
function onClickSpeciesSuggestions() {
    var featureNames = $('#searchFeatures').val();
    const newSpecies = this.id.slice("suggest".length);
    let requestData = {
        newSpecies: newSpecies,
        feature_names: featureNames,
        species: species,
    }
    $.ajax({
        type:'GET',
        url:'/data/by_celltype',
        data: $.param(requestData),
        success: function(result) {
            // Store global variable
            plotData = {
                'result': result,
            };
            $("#suggest"+newSpecies).text(species.slice(0, 1).toUpperCase()+species.slice(1)).prop('id', "suggest"+species);
            species = result['species'];

            // Update search box: corrected feature names, excluding missing features
            setSearchBox(result['features']);

            // Create heatmap
            updatePlot(true);
        },
        error: function (e) {
          alert('Error: Could not find orthologs for '+featureNames+'.')
        }
    });
}

function onClickFeatureSimilarSuggestions(targetType) {
    var featureNames = $('#searchFeatures').val();
    let requestData = {
        feature_names: featureNames,
        species: species,
        correlates_type: targetType,
        n_correlates: 5,
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
          alert('Error: Could not find features correlated with ' + featureNames + '.')
        }
    });
}


function onClickFeatureNearbySuggestions(targetType) {
    var featureNames = $('#searchFeatures').val();
    let requestData = {
        feature_names: featureNames,
        species: species,
        target_type: targetType,
        distance_max: 50000,
    }
    $.ajax({
        type:'GET',
        url:'/data/features_nearby',
        data: $.param(requestData),
        success: function(result) {
            // Update search box: corrected feature names, excluding missing ones
            setSearchBox(result);

            // Request data
            AssembleAjaxRequest();
        },
        error: function (e) {
          alert('Error: Could not find features near ' + featureNames + '.')
        }
    });
}


// Request genes in a clicked GO term and update plot
function onClickGOTermSuggestions () {
    let goTerm = $(this).text();
    let requestData = {
        goTerm: goTerm,
        species: species,
    }

    hideTooltip();

    $.ajax({
        type:'GET',
        url:'/data/genes_in_go_term',
        data: $.param(requestData),
        success: function(result) {

            // Update search box: corrected gene names, excluding missing genes
            // NOTE: this is a "short text" HTML element, so it can crash if we have
            // many genes...
            if (result.length > 200) {
                setSearchBox(
                    result.slice(0, 10)+"...",
                    result.split(',').slice(0, 10).join(','),
                );
                // Request data
                AssembleAjaxRequest(result);
            } else {
                setSearchBox(result);
                AssembleAjaxRequest();
            }
        },
        error: function (e) {
          alert('Error: Could not find genes within GO term '+goTerm+'.')
        }
    });
}

function setSearchBox(text) {
    $('#searchFeatures').val(text);
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
$(".speciesSuggestion").click(onClickSpeciesSuggestions);
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
    $('#pathwaySuggestion > a').click(function() {
        $("body").addClass("loading");
    });
    AssembleAjaxRequest();
});

