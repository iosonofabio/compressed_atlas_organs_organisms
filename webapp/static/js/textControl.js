function AssembleAjaxRequest() {
    var command = $('#textCommand').val();

    $.ajax({
        type: 'POST',
        url: '/submit_text',
        data: "text="+command,
        dataType: 'json',
        success: function(result) {
            if (result['outcome'] === 'success') {
                window.location.href = result['url'];
            } else {
                console.log(result);
            }
        },
        error: function(e) {
            console.log(e);
            alert('Text command not understood');
        }
    })
}


$("#askOnClick").click(function() {
  // action here when clicking the search button
  AssembleAjaxRequest();
});
$("body").keyup(function(event) {
    if (event.keyCode === 13) {
        $("#askOnClick").click();
    }
});
