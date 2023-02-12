function AssembleAjaxRequestTimepoint() {

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
    plot_type = "hieracical";
  }
  
  // action here when clicking the search button
  var gene_name = $('#searchGeneName').val();
  if (gene_name === '') {
    gene_name = "Col1a1";
  }

  const gene_array = gene_name.split(",")
    // sent gene names to the API
  $.ajax({
    type:'GET',
    url:'/data_timepoint',
    data: "gene=" + gene_name + "&plottype=" + plot_type + "&datatype=" + data_type,
    success: function(result) {
        num = 1
        for (key in Object.keys(result)) {
            dataset = Object.keys(result)[key]
            let div_id = 'dataset_' + num
            HeatMapTimepoint(result[dataset], div_id,dataset);
            num++
        }
    },
    error: function (e) {
        alert('Request data Failed(in assemble timepoint)')
    }
    });
  }
$("#searchOnClick" ).click(AssembleAjaxRequestTimepoint)
