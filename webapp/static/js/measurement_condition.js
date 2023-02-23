import { dotPlotSizeToFrac, dotPlotFracToSize, getDomains, getPseudocount, getTickTexts } from './plotUtils.js';

// Plot heatmap by celltype as a callback for the AJAX request
// Use global variables to store persistent data
var plotData = {};

function plotCondition(
    result,
    dataScale,
    tableOrder) {

    if ((!result) | (result.length == 0)) {
        alert("Error: Nothing to plot or feature names invalid")
        return;
    }

    const htmlElementId = "plotDiv";
    let htmlElement = document.getElementById(htmlElementId);

    let nPlots = result.length;
    let x_axis, y_axiss;
    if (tableOrder == "original") {
        x_axis = result[0]['celltypes'];
        y_axiss = [];
        for (let k=0; k < nPlots; k++) {
            y_axiss.push(result[k]['features']);
        }
    } else {
        // FIXME: this is not correct
        x_axis = result[0]['celltypes_hierarchical'];
        y_axiss = [];
        for(let k=0; k < nPlots; k++) {
            y_axiss.push(result[k]['features_hierarchical']);
        }
    }

    let longestXlabel = 0, longestYlabel = 0;
    for (let i=0; i < x_axis.length; i++) {
        longestXlabel = Math.max(longestXlabel, x_axis[i].length);
    }
    for (let k=0; k < nPlots; k++) {
        for (let i=0; i < y_axiss[k].length; i++) {
            longestYlabel = Math.max(longestYlabel, y_axiss[k][i].length);
        }
    }

    let nfeatures = y_axiss.reduce((acc, a) => acc + a.length, 0);
    let ncelltypes = x_axis.length;
    let pxCell = 40, pxChar = 4.2, plotGap = 10;
    let ytickMargin = 75 + pxChar * longestYlabel;
    let xtickMargin = 15 + pxChar * longestXlabel;
    let graphWidth = ytickMargin + pxCell * ncelltypes + 60;
    let graphHeight = pxCell * nfeatures + plotGap * (nPlots - 1) + xtickMargin;

    // Height ratios for the plots
    let yAxisDomains = getDomains(y_axiss, true);

    // Add hyperlinks to feature names if they are genes
    let yticktexts = [];
    for (let k=0; k < nPlots; k++) {
        let yticktexts_k = getTickTexts(
            y_axiss[k], result[k]['feature_type'], result[k]['gene_ids']);
        yticktexts.push(yticktexts_k);
    }

    var layout = {
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
    let traces = [];
    for (let k=0; k < nPlots; k++) {
        let yaxisName = 'yaxis', yaxisShort = 'y';
        if (k != 0) {
            yaxisName += (k+1);
            yaxisShort += (k+1);
        }
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

        // Fill data for this trace (feature type)
        let z = [];
        let measurement;
        const pseudoCount = getPseudocount(result[k]['feature_type']);
        for (let i=0; i < y_axiss[k].length; i++) {
            z.push([]);
            const feature = y_axiss[k][i];
            for (let j=0; j < x_axis.length; j++) {
                const celltype = x_axis[j];
                if (dataScale == "log2FC") {
                    measurement = result[k]['data'][feature][celltype];
                    const measBaseline = result[k]['data_baseline'][feature][celltype];
                    measurement = Math.log2(measurement + pseudoCount) - Math.log2(measBaseline + pseudoCount);
                } else {
                    if (dataScale.endsWith("baseline")) {
                        measurement = result[k]['data_baseline'][feature][celltype];
                    } else {
                        measurement = result[k]['data'][feature][celltype];
                    }
                    if (dataScale.startsWith("log10")) {
                        measurement = Math.log10(measurement + pseudoCount);
                    }
                }
                z[i].push(measurement);
            }
        }

        let trace = {
            z: z,
            x: x_axis,
            y: y_axiss[k],
            yaxis: yaxisShort,
            type: 'heatmap',
            hoverongaps: false,
        }
        if (dataScale == "log2FC") {
            trace['colorscale'] = 'RdBu';
            trace['zmid'] = 0;
        } else {
            trace['colorscale'] = 'Reds';
            trace['zmid'] = '';
        }
        traces.push(trace);
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
        
    // Create or update plot
    if (($('#'+htmlElementId).html() === "")) {
        Plotly.newPlot(htmlElement, traces, layout, config);
    } else {
        Plotly.react(htmlElement, traces, layout, config);
    }
} 


// gene of interest: Car4,Vwf,Col1a1,Ptprc,Ms4a1
// Col1a1,Fsd1l
function AssembleAjaxRequest() {

    // Get the list of genes to plot from the search box
    var feature_names = $('#searchFeatures').val();
    let requestData = {
        feature_names: feature_names,
        species: species,
    }

    // sent gene names to the API
    $.ajax({
        type:'GET',
        url:'/data/condition',
        data: $.param(requestData),
        dataType:'json',
        success: function(result) {
            plotData = result;
            // NOTE: do not plot already, let the user choose the time point
            //updatePlot();
        },
        error: function (e) {
            console.log(e);
            alert('Request data Failed');
        }
    });

};


function updatePlot() {
    // NOTE: plotData is the global persistent object
    let dataScale, tableOrder;
    if ($("#log2FCTab").hasClass('is-active')) {
        dataScale = 'log2FC';
    } else if ($("#cpmTab").hasClass('is-active')) {
        dataScale = "cpm";
    } else if ($("#cpmBaselineTab").hasClass('is-active')) {
        dataScale = "baseline";
    } else if ($("#logBaselineTab").hasClass('is-active')) {
        dataScale = "log10baseline";
    } else {
        dataScale = "log10";
    }
    
    if ($("#originalOrderTab").hasClass('is-active')) {
      tableOrder = "original";
    } else {
      tableOrder = "hierachical";
    }

    // Plot the time point chosen in the dropdown
    // (they are all cached for convenience)
    // There should be one item per feature type (e.g. chromatin accessibility)
    let itemsToPlot = [];
    for(let i = 0; i < plotData.length; i++) {
        let item = plotData[i];
        let dsTp = item['dataset']+'_'+item['timepoint'];
        let htmlElement = document.getElementById("dropdownItem_"+dsTp);
        if (!htmlElement)
            continue;
        if (htmlElement.classList.contains('is-active')) {
            itemsToPlot.push(item);
        }
    }
    plotCondition(
          itemsToPlot, 
          dataScale,
          tableOrder,
        );
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

// dataScale, tableOrder button callbacks
function triggerButtons() {
    let buttonClasses = ["dataScaleButton", "tableOrderButton"];
    for (let i=0; i < buttonClasses.length; i++) {
        $("."+buttonClasses[i]).click(function() {
            $("."+buttonClasses[i]).removeClass('is-active');
            $(this).addClass('is-active');
            updatePlot();
        });
    }
}
triggerButtons();
