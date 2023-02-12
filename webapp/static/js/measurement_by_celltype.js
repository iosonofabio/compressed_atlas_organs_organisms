// Plot heatmap by celltype as a callback for the AJAX request
// Use global variables to store persistent data
var heatmapData = {};

function HeatmapByCelltype(
    result_input,
    htmlElementId,
    dataIndex,
    dataScale,
    celltypeOrder,
    heatDot) {

    // Only the bottom plot requires xtick labels
    let isBottomPlot = (dataIndex + 1) == result_input['data'].length;

    let result = {
        'data': result_input['data'][dataIndex],
        'features': result_input['features'][dataIndex],
        'feature_type': result_input['feature_type'][dataIndex],
        'features_hierarchical': result_input['features_hierarchical'][dataIndex],
        'celltypes': result_input['celltypes'],
        'celltypes_hierarchical': result_input['celltypes_hierarchical'],
        'species': result_input['species'],
    }
    if (result['feature_type'] == 'gene_expression') {
        result['data_fractions'] = result_input['data_fractions'][dataIndex]
        result['gene_ids'] = result_input['gene_ids'][dataIndex]
        result['GO_terms'] = result_input['GO_terms'][dataIndex]
    }

    let x_axis, y_axis;
    let longestXlabel = 0, longestYlabel = 0;
    if (celltypeOrder == "original") {
        x_axis = result['celltypes'];
        y_axis = result['features'];
    } else {
        x_axis = [];
        for (let i = 0; i < result['celltypes_hierarchical'].length; i++) {
            const ct = result['celltypes'][result['celltypes_hierarchical'][i]];
            longestXlabel = Math.max(longestXlabel, ct.length);
            x_axis.push(ct);
        }
        y_axis = [];
        for (let i = 0; i < result['features_hierarchical'].length; i++) {
            const feature = result['features'][result['features_hierarchical'][i]];
            y_axis.push(feature);
            longestYlabel = Math.max(longestYlabel, feature.length);
        }
    }
    let ytickMargin = 200;

    let nfeatures =  y_axis.length;
    let ncelltypes = x_axis.length;
    let pxCell = 40, pxChar = 15;
    let graph_width = pxChar * (5 + longestYlabel) + pxCell * ncelltypes + 120;
    let graph_height = pxCell * nfeatures + pxChar * longestXlabel + 20;

    // Add hyperlinks to feature names if they are genes
    let yticktext = [];
    for (let i = 0; i < y_axis.length; i++) {
        const feature = y_axis[i];
        let tickText;
        if (result['feature_type'] == 'gene_expression') {
            const geneId = result['gene_ids'][feature];
            if (geneId !== "") {
                let geneUrl = feature;
                if (geneId.startsWith('MGI')) {
                    geneUrl = 'http://www.informatics.jax.org/marker/'+geneId;
                } else {
                    geneUrl = 'https://www.genecards.org/cgi-bin/carddisp.pl?gene='+geneId;
                }
                tickText = '<a href="'+geneUrl+'">'+feature+'</a> - <span><b>GO</b></span>';
            }
        } else {
            tickText = feature.split('-');
        }
        yticktext.push(tickText);
    }

    // Add SVG download button
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
        // Fill data
        let data_content = [];
        if (celltypeOrder == "original") {
            for (let i = 0; i < y_axis.length; i++) {
                data_content.push([]);
                for (let j = 0; j < x_axis.length; j++) {
                    let measurement = result['data'][i][j];
                    if (dataScale == "log10") {
                        let pseudocount = 0.5;
                        if (result['feature_type'] != 'gene_expression') {
                            pseudocount = 0.01;
                        }
                        measurement = Math.log10(measurement + pseudocount);
                    }
                    data_content[i].push(measurement);
                }
            }
        } else {
            for (let i = 0; i < y_axis.length; i++) {
                const ii = result['features_hierarchical'][i];
                data_content.push([]);
                for (let j = 0; j < x_axis.length; j++) {
                    const jj = result['celltypes_hierarchical'][j];
                    let measurement = result['data'][ii][jj];
                    if (dataScale == "log10") {
                        let pseudocount = 0.5;
                        if (result['feature_type'] != 'gene_expression') {
                            pseudocount = 0.01;
                        }
                        measurement = Math.log10(measurement + pseudocount);
                    }
                    data_content[i].push(measurement);
                }
            }
        }
        var data = {
            type: 'heatmap',
            hoverongaps: false,
            colorscale: 'Reds',
        };

        // Make new plot if none is present
        if (($('#'+htmlElementId).html() === ""))
            plotForceRefresh = true;
        // If it's a new plot or a refresh was forced
        if (plotForceRefresh == true) {
            data['z'] = data_content;
            data['x'] = x_axis;
            data['y'] = y_axis;

            var layout = {
                autosize: true,
                width: graph_width,
                height: graph_height,
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
                    constraintoward: 'left',
                },
                yaxis: {
                    automargin: false,
                    autorange: "reversed",
                    type: 'category',
                    scaleratio: 1,
                    scaleanchor: 'x',
                    tickvals: y_axis,
                    ticktext: yticktext,
                },
            };

            config["modeBarButtonsToAdd"].push({
                name: 'Download expression as CSV',
                icon: Plotly.Icons.disk,
                click: function(gd) {
                    var text = '';
                    let data = gd['data'][0]['z'];
                    text += 'Gene,' + gd['data'][0]['x'] + '\n';
                    for(var i = 0; i < data.length; i++){
                        text += gd['data'][0]['y'][i] + ',' + data[i] + '\n';
                    };
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

            Plotly.newPlot(
                document.getElementById(htmlElementId),
                [data],
                layout,
                config,
            );

        // Update existing plot if present
        } else {
            data['z'] = [data_content];
            data['x'] = [x_axis];
            data['y'] = [y_axis];
            Plotly.update(
                document.getElementById(htmlElementId),
                data,
                {
                    height: graph_height,
                    yaxis: {
                        autorange: "reversed",
                        automargin: false,
                        type: 'category',
                        scaleratio: 1,
                        scaleanchor: 'x',
                        tickvals: y_axis,
                        ticktext: yticktext,
                    },
                },
                [0],
            );
        }

    // Dot plot
    } else {

        let x = [], y = [], tooltips = [], markersize = [], markercolor = [];
        if (celltypeOrder == "original") {
            for (let i = 0; i < y_axis.length; i++) {
                for (let j = 0; j < x_axis.length; j++) {
                    const feature = y_axis[i];
                    const celltype = x_axis[j];
                    let measurement = result['data'][i][j];
                    if (dataScale == "log10") {
                        let pseudocount = 0.5;
                        if (result['feature_type'] != 'gene_expression') {
                            pseudocount = 0.01;
                        }
                        measurement = Math.log10(measurement + pseudocount);
                    }
                    let frac;
                    if (result['feature_type'] == 'gene_expression') {
                        frac = result['data_fractions'][i][j];
                    } else {
                        frac = result['data'][i][j];
                    }

                    // ATAC-Seq has generally smaller fractions, so highlight more
                    let ms;
                    if (result['feature_type'] == 'gene_expression') {
                        ms = 2 + 18 * Math.sqrt(frac);
                    } else {
                        ms = 2 + 60 * Math.sqrt(frac);
                    }

                    const tooltip = "Average: "+measurement+", Fraction of cells: "+parseInt(100 * frac)+"%";
                    x.push(celltype);
                    y.push(feature);
                    markercolor.push(measurement);
                    markersize.push(ms);
                    tooltips.push(tooltip);
                }
            }
        } else {
            for (let i = 0; i < y_axis.length; i++) {
                const ii = result['features_hierarchical'][i];
                for (let j = 0; j < x_axis.length; j++) {
                    const jj = result['celltypes_hierarchical'][j];
                    const feature = y_axis[i];
                    const celltype = x_axis[j];
                    let measurement = result['data'][ii][jj];
                    if (dataScale == "log10") {
                        let pseudocount = 0.5;
                        if (result['feature_type'] != 'gene_expression') {
                            pseudocount = 0.01;
                        }
                        measurement = Math.log10(measurement + pseudocount);
                    }
                    let frac;
                    if (result['feature_type'] == 'gene_expression') {
                        frac = result['data_fractions'][ii][jj];
                    } else {
                        frac = result['data'][ii][jj];
                    }
                    let ms;
                    if (result['feature_type'] == 'gene_expression') {
                        ms = 2 + 18 * Math.sqrt(measurement);
                    } else {
                        ms = 2 + 60 * Math.sqrt(frac);
                    }
                    const tooltip = "Average: "+measurement+", Fraction of cells: "+parseInt(100 * frac)+"%";
                    x.push(celltype);
                    y.push(feature);
                    markercolor.push(measurement);
                    markersize.push(ms);
                    tooltips.push(tooltip);
                }
            }
        }
        var data = {
            mode: 'markers',
            marker: {
                symbol: 'circle',
                colorscale: 'Reds',
                colorbar: {},
            },
            'hoverinfo': 'text',
        };

        if (($('#'+htmlElementId).html() === "") || (plotForceRefresh == true)) {
            data['x'] = x;
            data['y'] = y;
            data['text'] = tooltips;
            data['marker']['color'] = markercolor;
            data['marker']['size'] = markersize;

            var layout = {
                autosize: true,
                width: graph_width,
                height: graph_height,
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
                yaxis: {
                    autorange: "reversed",
                    type: 'category',
                    automargin: false,
                    scaleanchor: 'x',
                    scaleratio: 1,
                    tickvals: y_axis,
                    ticktext: yticktext,
                },
            };

            config["modeBarButtonsToAdd"].push({
                name: 'Download average as CSV',
                icon: Plotly.Icons.disk,
                click: function(gd) {
                    var text = '';
                    let measurements = gd['data'][0]['marker']['color'];
                    const nct = x_axis.length;
                    // Header with cell type names
                    text += 'Feature';
                    for(var i = 0; i < nct; i++){
                        text += ',' + gd['data'][0]['x'][i];
                    };
                    // Gene expression
                    for (var i = 0; i < measurements.length; i++) {
                        if (i % nct == 0) {
                            text += '\n' + gd['data'][0]['y'][i];
                        }
                        text += ',' + measurements[i];
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
                    var text = '';
                    let markerSizes = gd['data'][0]['marker']['size'];
                    const nct = x_axis.length;
                    // Header with cell type names
                    text += 'Feature';
                    for(var i = 0; i < nct; i++){
                        text += ',' + gd['data'][0]['x'][i];
                    };
                    // Gene expression
                    for (var i = 0; i < markerSizes.length; i++) {
                        if (i % nct == 0) {
                            text += '\n' + gd['data'][0]['y'][i];
                        }
                        let frac = (markerSizes[i] - 2) / 18.0;
                        // The marker radius is prop to the sqrt
                        frac *= frac;
                        text += ',' + frac;
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

            Plotly.newPlot(
                document.getElementById(htmlElementId),
                [data],
                layout,
                config,
            );

        // Update existing plot
        } else {
            data['x'] = [x];
            data['y'] = [y];
            data['text'] = [tooltips];
            data['marker']['color'] = markercolor;
            data['marker']['size'] = markersize;
            Plotly.update(
                document.getElementById(htmlElementId),
                data,
                {
                    height: graph_height,
                    yaxis: {
                        autorange: "reversed",
                        type: 'category',
                        automargin: false,
                        scaleanchor: 'x',
                        scaleratio: 1,
                        tickvals: result['yticks'],
                        ticktext: result['yticktext'],
                    },
                    xaxis: {
                        autorange: true,
                        automargin: true,
                        type: 'category',
                        tickangle: 270,
                    },

                },
                [0],
            );
        }
    }

    // Add tooltips to gene names
    $(".ytick > text > tspan").click(function(evt) {
        // If already active, it's a second click
        let wasActive = this.classList.contains("active");

        // Deactivate all, and if this was not active activate it
        hideTooltip();
        $(".ytick > text > tspan").map(function(elem) {
            this.classList.remove("active");
            let gene = this.parentElement.getElementsByTagName("a")[0].innerHTML;
        })
        if (wasActive) {
            return;
        }

        // Activate and show pop-up window
        this.classList.add("active");
        let gene = this.parentElement.getElementsByTagName("a")[0].innerHTML;
        let goTerms = result['GO_terms'][gene];
        if (goTerms === undefined) {
            return;
        }
        let text = "<b>GO terms:</b></br>";
        for (let i = 0; i < goTerms.length; i++) {
            text += '<div><a class="goHyperlink">'+goTerms[i]+"</a></div>";
        }
        showTooltip(evt, text);
    });
    //$(".ytick > text > tspan").mouseout(function(evt) { hideTooltip(); });
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
function updatePlot() {
    let dataScale = "original";
    if (!$("#cpmTab").hasClass('is-active')) {
        dataScale = "log10";
    }
    let celltypeOrder = "original";
    if (!$("#originalOrderTab").hasClass('is-active')) {
        celltypeOrder = "hierarchical";
    }
    let heatDot = "heat";
    if (!$("#heatTab").hasClass('is-active')) {
        heatDot = "dot";
    }

    // NOTE: heatmapData is the global persistent object
    if (!heatmapData['result']) {
        alert("Error: no data to plot");
        return;
    }

    // Empty divs that are not plotted anymore
    const plotTypes = ['gene_expression', 'chromatin_accessibility'];
    let plotActive = {}
    for (let i=0; i < plotTypes.length; i++) {
        plotActive[plotTypes[i]] = false;
    }
    for (let i=0; i < heatmapData['result']['n_feature_types']; i++) {
        plotActive[heatmapData['result']['feature_type'][i]] = true;
    }
    for (let i=0; i < plotTypes.length; i++) {
        if (!plotActive[plotTypes[i]]) {
            $('#plot_'+plotTypes[i]).html("");
        }
    }

    // Fill divs that are plotted
    for (let i=0; i < heatmapData['result']['n_feature_types']; i++) {
        HeatmapByCelltype(
            heatmapData['result'],
            'plot_'+heatmapData['result']['feature_type'][i],
            i,
            dataScale,
            celltypeOrder,
            heatDot,
        );
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
            heatmapData = {
                'result': result,
            };

            // Update search box: corrected feature names, excluding missing features
            setSearchBox(result['features']);

            // Create heatmap
            updatePlot();
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
            heatmapData = {
                'result': result,
            };
            $("#suggest"+newSpecies).text(species.slice(0, 1).toUpperCase()+species.slice(1)).prop('id', "suggest"+species);
            species = result['species'];

            // Update search box: corrected feature names, excluding missing features
            setSearchBox(result['features']);

            // Create heatmap
            updatePlot();
        },
        error: function (e) {
          alert('Error: Could not find orthologs for '+featureNames+'.')
        }
    });
}

// SuggestGenes: create a div with a "suggest" button
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

function onClickGeneSuggestions() {
    return onClickFeatureSuggestions(correlatesType="gene_expression");
}
function onClickRegionSuggestions() {
    return onClickFeatureSuggestions(correlatesType="chromatin_accessibility");
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

function setSearchBox(text, gseaText = "") {
    $('#searchFeatures').val(text);
    // Sync with GSEA box
    if (gseaText == "") {
        gseaText = text;
    }
    $('#suggestGO > a').attr(
        'href', '/barplot_gsea?species='+species+'&genes='+gseaText);
    $('#suggestKEGG > a').attr(
        'href', '/barplot_gsea?species='+species+'&gene_set=KEGG&genes='+gseaText);
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
    plotForceRefresh = false;
    updatePlot();
});

$("#CPMOnClick" ).click(function() {
    $("#logTab").removeClass('is-active');
    $("#cpmTab").addClass('is-active');
    plotForceRefresh = false;
    updatePlot();
});

// Second set of buttons
$("#hClusterOnClick" ).click(function() {
    // if User has input their features of interest, generate the heatmap with that value
    // otherwise use the default data
    $("#hierachicalTab").addClass('is-active');
    $("#originalOrderTab").removeClass('is-active');
    plotForceRefresh = false;
    updatePlot();
});


$("#originalOnClick").click(function() {
    $("#originalOrderTab").addClass('is-active');
    $("#hierachicalTab").removeClass('is-active');
    plotForceRefresh = false;
    updatePlot();
});

// Third set of buttons
$("#heatOnClick").click(function() {
    $("#heatTab").addClass('is-active');
    $("#dotTab").removeClass('is-active');
    plotForceRefresh = true;
    updatePlot();
});

$("#dotOnClick").click(function() {
    $("#dotTab").addClass('is-active');
    $("#heatTab").removeClass('is-active');
    plotForceRefresh = true;
    updatePlot();
});

// Both on click and load, plot the heatmap
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
$("#geneSuggestions").click(onClickGeneSuggestions);
$("#regionSuggestions").click(onClickRegionSuggestions);
$(".speciesSuggestion").click(onClickSpeciesSuggestions);
