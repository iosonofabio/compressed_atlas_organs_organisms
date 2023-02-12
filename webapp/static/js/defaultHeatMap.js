$(document).ready(function() {
    $.ajax({
        type:'GET',
        url:'/data/by_celltype',
        data: "gene_names=Car4,Vwf,Col1a1,Ptprc,Ms4a1&plot_type=original&data_type=original",
        dataType:'json',
        success: function(result) {
            HeatMap(result, "h5_data_plot");
        },
        error: function (e) {
        alert('Request data Failed')
        }
    });
    $("#originalTab").addClass('is-active');
})
