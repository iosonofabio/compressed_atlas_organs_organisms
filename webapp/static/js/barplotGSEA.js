// Plot heatmap by celltype as a callback for the AJAX request
// Use global variables to store persistent data
var plotData = {};

function barplotGSEA(result, html_element_id, xaxisType) {
    if (!result) {
        alert("Error: no data to plot");
        return;
    }

    let y_axis = result['pathways'];
    var npathways =  y_axis.length;
    var graph_width = 1300;
    var graph_height = 270 + 20 * npathways;

    // Add hyperlinks to heatmap of all genes in each pathway
    let yticktext = [];
    for (let i = 0; i < y_axis.length; i++) {
        const pathway = y_axis[i];
        let url;
        if (result['pathways_urls'].length > 0) {
            url = result['pathways_urls'][i];
        } else {
            url = '/heatmap_by_celltype?'+$.param({
                species: species,
                pathway: pathway,
            })
        }
        const tickText = '<a href="'+url+'">'+pathway+'</a>'
        yticktext.push(tickText);
    }

    // Fill heatmap data
    let data_content;
    if (xaxisType == "overlap") {
        data_content = [];
        for (let i = 0; i < result['overlap'].length; i++) {
            let ovArray = result['overlap'][i].split('/');
            let fraction = parseInt(ovArray[0]) * 1.0 / parseInt(ovArray[1]);
            data_content.push(fraction);
        }
    } else {
        data_content= result['neglog10_p_value'];
    }

    let data = {
        type: 'bar',
        orientation: 'h',
        marker: {
            colorscale: 'Reds',
        },
    };

    // Make new plot if none is present
    if ($('#'+html_element_id).html() === "") {

        data['x'] = data_content;
        data['y'] = y_axis;
        data['marker']['color'] = data_content;

        var layout = {
            autosize: true,
            width: graph_width,
            height: graph_height,
            title: 'Gene set enrichment analysis',
            xaxis: {
                automargin: true,
                tickangle: 70,
            },
            yaxis: {
                automargin: true,
                autorange: "reversed",
                type: 'category',
                tickvals: y_axis,
                ticktext: yticktext,
            },
        };
            
        Plotly.newPlot(
            document.getElementById(html_element_id),
            [data],
            layout,
        ); 

    // Update existing plot if present
    } else {
        data['x'] = [data_content];
        data['y'] = [y_axis];
        data['marker']['color'] = data_content;

        Plotly.update(
            document.getElementById(html_element_id),
            data,
            {
                height: graph_height,
                yaxis: {
                    automargin: true,
                    autorange: "reversed",
                    tickvals: y_axis,
                    ticktext: yticktext,
                },
            },
            [0],
        ); 
    }
} 


// NOTE: this is why react was invented...
function updatePlot() {
    let xaxisType = "p_value";

    if ($("#overlapTab").hasClass('is-active')) {
        xaxisType = "overlap";
    }

    // NOTE: plotData is the global persistent object
    barplotGSEA(
        plotData, 
        plotData['div'],
        xaxisType,
    );
}


////////////////////
// EVENTS
////////////////////
$(".xaxisType" ).click(function() {
    $(".xaxisType").removeClass('is-active');
    $(this).addClass('is-active');
    updatePlot();
});

$(document).ready(function() {
    plotData['div'] = 'data_plot';
    updatePlot()
});
