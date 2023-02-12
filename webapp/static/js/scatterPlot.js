function ScatterPlot(result) {
    // compare the expression level of 2 specific genes across all cell types
    let expr1 = result["gene1_expr"];
    let expr2 = result["gene2_expr"];
    let cell_types =result["cell_types"];
    var trace1 = {
        x: expr1,
        y: expr2,
        mode: 'markers+text',
        type: 'scatter',
        text: cell_types,
        textposition: 'bottom center',
        textfont: {
          family:  'Raleway, sans-serif'
        },
        marker: { size: 12 }
    };
        
    var data = [ trace1 ];
    
    var layout = {
        height: 800,
        width: 800,
        xaxis: {
            title: {
              text: result['gene1_name']+' expression',
            },
            //range: [0, 8]
        },
        yaxis: {
            title: {
              text: result['gene2_name']+' expression',
          },
          //range: [0, 8]
        },
        title:'Correlation between ' + result['gene1_name'] + ' and ' + result['gene2_name']
    };
    
    Plotly.newPlot(document.getElementById('scatter_plot'), data, layout);
}
