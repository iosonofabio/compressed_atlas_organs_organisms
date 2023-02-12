// Plot heatmap by celltype as a callback for the AJAX request
function barplotCelltype(result, html_element_id) {
    //TODO
        if (!result) {
            alert("Error:Input gene name is invalid, please make sure you type in the corrent gene names")
        } else {
            // x-axis: genes of interest
            var x_axis = Object.keys(result[Object.keys(result)[0]]);
            var y_axis = Object.keys(result);
            var ngenes =  y_axis.length;
            var graph_width = 1300;
            var graph_height = 270 + 26 * ngenes;
            // y-axis:41 cell types
            var data_content = [];
            for (var i = 0; i < Object.keys(result).length; i++) {
                cell_type = Object.keys(result)[i] // get the cell_type name as a string
                all_gene_expression = result[cell_type]         // find it from the dictionary as a key
                
                data_content.push(Object.values(all_gene_expression))
            }
            var data = [
                {
                    z: data_content,
                    x: x_axis,
                    y: y_axis,
                    type: 'heatmap',
                    hoverongaps: false,
                    colorscale: 'Reds',
                }
                ];
            var layout = {
                autosize: true,
                width: graph_width,
                height: graph_height,
                title: 'Heatmap of gene expression level in selected cell types',
                xaxis: {
                    //title: 'Cell types',
                    automargin: true,
                    tickangle: 70,
                    scaleanchor: 'y',
                    scaleratio: 1,
                },
                yaxis: {
                    //title: 'Genes',
                    automargin: true,
                    autorange: "reversed",
                },
            };
                
            Plotly.newPlot(
                document.getElementById(html_element_id),
                data,
                layout,
            ); 
        };
    } 


// SuggestGenes: create a div with a "suggest" button
function onClickSuggestions() {
    var gene_names = $('#searchGeneName').val();
    $.ajax({
        type:'GET',
        url:'/gene_friends',
        data: "gene_names=" + gene_names,
        success: function(result) {
            $('#searchGeneName').val(result);
            //if (("#"+html_element_id).text() == "Suggest similar genes") {
            //    $("#"+html_element_id).text("Suggest MORE similar genes");
            //    // NOTE: the right function is already bound
            //} else if (("#"+html_element_id).text() == "Suggest MORE similar genes") {
            //    // Stop the search after two layers
            //    $("#"+html_element_id).text("");
            //    // FIXME: there must be a more elegant way to unbind...
            //    $("#"+html_element_id).click(function() {});
            //}
            AssembleAjaxRequest();
        },
        error: function (e) {
          alert('Error: Could not find gene friends for '+gene_names+'.')
        }
    });
}

function SuggestGenes(html_element_id) {

    // Fill div
    $("#"+html_element_id).empty();
    $("#"+html_element_id).text('Suggest similar genes');

    // Add link
    $("#"+html_element_id).click(onClickSuggestions);
}

// gene of interest: Car4,Vwf,Col1a1,Ptprc,Ms4a1
// Col1a1,Fsd1l
function AssembleAjaxRequest() {
  if(! $('#scatter_plot').is('empty')) {
    $('#scatter_plot').empty();
  }

  // Get the list of genes to plot from the search box
  var gene_names = $('#searchGeneName').val();
  if (gene_names === '') {
    gene_names = "Car4,Vwf,Col1a1,Ptprc,Ms4a1";
  }

  // When doing the search gene name action, we want it to be change immediatly without switching back to the original heatmap,
  //  for example, if we are looking at a log10 plot,and we do the search action, the tab stays at the log10 
  cpm_is_active = $("#cpmTab").hasClass('is-active');
  orginal_is_active = $("#originalOrderTab").hasClass('is-active')
  var plot_type = 'original';
  var data_type = 'original';
  
  if (!cpm_is_active) {
    data_type = 'log10';
  } 
  
  if (!orginal_is_active) {
    plot_type = "hierachical";
  }
  const gene_array = gene_names.split(",")
  if (gene_array.length == 2) {
    $.ajax({
      type:'GET',
      url:'/2_genes',
      data: "gene_names=" + gene_names,
      success: ScatterPlot,
      error: function (e) {
        alert('Request data Failed')
      }
    });
  }
    // sent gene names to the API
  $.ajax({
    type:'GET',
    url:'/data/by_celltype',
    data: "gene_names=" + gene_names + "&plot_type=" + plot_type + "&data_type=" + data_type,
    success: function(result) {  
      // Create heatmap
      HeatMap(result, "h5_data_plot");

      // Create gene suggestions DOM element
      SuggestGenes("gene_suggestions")

      //FIXME: ScatterPlot should happen here
    },
    error: function (e) {
      alert('Error:Input gene name is invalid, please make sure you type in the corrent gene names.')
    }
    });
};

// Both on click and load, plot the heatmap
$("#searchOnClick").click(function() {
  // action here when clicking the search button
  AssembleAjaxRequest();
});

$(document).ready(function() {
  AssembleAjaxRequest();
});
