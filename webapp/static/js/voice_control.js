let chunks = [];
var rec;
var stream;
let input;
let AudioContext = window.AudioContext || window.webkitAudioContext;
let audioContext = new AudioContext;

//add events to the button
$("#recordButton").mousedown(() => {
    $(window).off();
    $(window).mouseup(stopRecording);
    console.log("recordButton clicked");

    // Switch to stop button until we hear back from getUserMedia
    $("#recordButton").attr("src", "/static/images/stop_button_navbar.png")
                      .attr("alt", "stop button");

    //// Create audio context and processor to resample at 16kHz, to make google happy
    //AudioContext = window.AudioContext || window.webkitAudioContext;
    //context = new AudioContext({
    //  // if Non-interactive, use 'playback' or 'balanced' // https://developer.mozilla.org/en-US/docs/Web/API/AudioContextLatencyCategory
    //  latencyHint: 'interactive',
    //});
    //processor = context.createScriptProcessor(bufferSize, 1, 1);
    //processor.connect(context.destination);
    //context.resume();

    var constraints = { audio: true, video:false };
    navigator.mediaDevices.getUserMedia(constraints).then(localStream => {
        console.log("getUserMedia() success, stream created, ready for recording ...");

        // make stream and recorder global
        stream = localStream;

        input = audioContext.createMediaStreamSource(stream);
        rec = new Recorder(input);
        rec.record()

        console.log("Recording started");

    }).catch(function(err) {
        //enable the record button if getUserMedia() fails
        $("#recordButton").attr("src", "/static/images/rec_button_navbar.png")
                          .attr("alt", "rec button");
    });
});


function stopRecording() {
    console.log("recordButton released");

    // Remove the mouseup event handler
    $(window).off();

    // enable the record button back to signal the user we got the message
    $("#recordButton").attr("src", "/static/images/rec_button_navbar.png")
                      .attr("alt", "rec button");

    //tell the recorder to stop the recording
    rec.stop();
    stream.getAudioTracks().forEach(track => track.stop());

    // export to WAV format (Google is picky)
    rec.exportWAV(postAudio);

    console.log("end of stopRecording");
}


function clarifyResponse(response) {
    // Gene string unclear, open an input popup asking to correct spelling
    if (response['question'] === 'gene_string') {
        let genestr = response['gene_string'];

        // Open confirmation popup
        genestr = window.prompt("Confirm or correct gene names:", genestr);
        if (genestr == null || genestr == "") {
            return;
        }

        // Validate new gene string
        $.ajax({
            type:'GET',
            url:'/check_genenames',
            data: "gene_names="+genestr,
            success: function(result) {
                if (result['outcome'] == 'fail') {
                    console.log("check for genestring found no match");
                    return;
                }

                let url = response['url_prefix'] + result['genenames'];
                window.location.href = url;
            },
            error: function (e) {
               console.log("check for genestring failed");
            }
        });
    } else if (response['question'] === 'celltype_string') {
        let ctstr = response['celltype_string'];

        // Open confirmation popup
        ctstr = window.prompt("Confirm or correct cell type names:", ctstr);
        if (ctstr == null || ctstr == "") {
            return;
        }

        // Validate new cell type string
        $.ajax({
            type:'GET',
            url:'/data/marker_genes',
            data: "celltype_names="+ctstr,
            success: function(result) {
                if (result['outcome'] == 'fail') {
                    console.log("check for ctstr found no match");
                    return;
                }

                let url = response['url_prefix'] + result['genenames'];
                window.location.href = url;
            },
            error: function (e) {
               console.log("check for genestring failed");
            }
        });
    }
}


function postAudio(blob) {
    // Do this with jQuery?
    var xhr=new XMLHttpRequest();
    xhr.onload=function(e) {
        if(this.readyState === 4) {
            var response = JSON.parse(e.target.responseText);
            console.log("Server returned: ", response);

            // Empty response, failed
            if (response['outcome'] == "fail") {
                alert("Voice command not understood.");

            // Request for clarification from the user
            } else if (response['outcome'] == "question") {
                clarifyResponse(response);

            // Successful response, go there
            } else {
                window.location.href = response['url'];
            }
        }
    };

    // This uses a form with POST, a little dangerous...
    var fd = new FormData();
    fd.append("audio_data", blob);
    xhr.open("POST","/submit_audio", true);
    xhr.send(fd);

}
