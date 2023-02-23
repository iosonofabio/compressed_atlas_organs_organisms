import { getDomains, getPseudocount, getTickTexts } from './plotUtils.js';

// This is injected into the template for some reason, perhaps optimization
//let plotData = {};

// Plot heatmap by celltype as a callback for the AJAX request
// Use global variables to store persistent data
function plotDifferential(
    result,
    dataScale,
    tableOrder) {

    if (!result) {
        alert("Error: Nothing to plot or gene names invalid")
        return;
    }

    const htmlElementId = "plotDiv";
    let htmlElement = document.getElementById(htmlElementId);

    let nPlots = result['feature_type'].length;
    let x_axis, y_axiss;
    if (tableOrder == "original") {
        x_axis = result['celltypes'];
        y_axiss = result['features'];
    } else {
        x_axis = result['celltypes_hierarchical'];
        y_axiss = result['features_hierarchical'];
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
    let pxCell = 40, pxChar = 4.2, plotGap = 10;
    let ytickMargin = 85 + pxChar * longestYlabel;
    let xtickMargin = 15 + pxChar * longestXlabel;
    let graphWidth = ytickMargin + pxCell * ncelltypes + 60;
    let graphHeight = pxCell * nfeatures + plotGap * (nPlots - 1) + xtickMargin;

    // Height ratios for the plots
    let yAxisDomains = getDomains(y_axiss, true);

    // Add hyperlinks to gene names
    let yticktexts = [];
    for (let k=0; k < nPlots; k++) {
        let yticktext_k = getTickTexts(y_axiss[k], result['feature_type'][k], result['gene_ids'][k]);
        yticktexts.push(yticktext_k);
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

        // Fill heatmap data
        let z = [];
        let measurement;
        const pseudoCount = getPseudocount(result['feature_type'][k]);
        for (let i=0; i < y_axiss[k].length; i++) {
            z.push([]);
            const feature = y_axiss[k][i];
            for (let j=0; j < x_axis.length; j++) {
                const celltype = x_axis[j];
                measurement = result['data'][k][feature][celltype];
                if (dataScale == "log2FC") {
                    measurement = result['data'][k][feature][celltype];
                    let measBaseline;
                    if (result['comparison'] == 'celltypes') {
                        measBaseline = result['data_baseline'][k][feature][celltype2];
                    } else {
                        measBaseline = result['data_baseline'][k][feature][celltype];
                    }
                    measurement = Math.log2(measurement + pseudoCount) - Math.log2(measBaseline + pseudoCount);
                } else {
                    if (dataScale.endsWith("baseline")) {
                        measurement = result['data_baseline'][k][feature][celltype];
                    } else {
                        measurement = result['data'][k][feature][celltype];
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

    //let title;
    //if (dataScale == "log2FC") {
    //    title = 'Differential expression ';
    //    if (plotData['comparison'] == 'hyperoxia/normal') {
    //        title += 'hyperoxia vs normal at '+plotData['timepoint'];
    //    } else if (plotData['comparison'] == 'timepoints') {
    //        title += plotData['timepoint']+' vs '+plotData['timepoint_baseline'];
    //    } else if (plotData['comparison'] == 'celltypes') {
    //        title += ' at '+plotData['timepoint']+' in '+plotData['celltype']+' vs '+plotData['celltype_baseline'];
    //    }
    //}

    if ($('#'+htmlElementId).html() === "") {
        Plotly.newPlot(htmlElement, traces, layout, config);
    } else {
        Plotly.react(htmlElement, traces, layout, config);
    }
} 


// NOTE: this is why react was invented...
function updatePlot() {
    let dataScale, tableOrder;
    if ($("#log2FCTab").hasClass('is-active')) {
        dataScale = "log2FC";
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

    // NOTE: plotData is the global persistent object
    plotDifferential(
        plotData, 
        dataScale,
        tableOrder,
    );
}


// gene of interest: Car4,Vwf,Col1a1,Ptprc,Ms4a1
// Col1a1,Fsd1l
function AssembleAjaxRequest() {

    // Get the list of genes to plot from the search box
    let featureNames = $('#searchFeatures').val();
    // NOTE: you cannot cache the genes because the hierarchical clustering
    // will differ anyway

    let requestData = {
        comparison: plotData['comparison'],
        ct1: plotData['celltype'],
        ct2: plotData['celltype_baseline'],
        ds1: plotData['dataset'],
        ds2: plotData['dataset_baseline'],
        tp1: plotData['timepoint'],
        tp2: plotData['timepoint_baseline'],
        dis1: plotData['disease'],
        dis2: plotData['disease_baseline'],
        featurestr: featureNames,
        species: species,
    }

    // sent conditions and gene names to the API
    $.ajax({
        type:'GET',
        url:'/data/differential',
        data: $.param(requestData),
        dataType:'json',
        success: function(result) {
            console.log(result);

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
  AssembleAjaxRequest();
});

// On load, the heatmap data are already embedded in the template
$(document).ready(function() {
    if (plotData == {}) {
        AssembleAjaxRequest();
    } else {
        updatePlot();
    }
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
