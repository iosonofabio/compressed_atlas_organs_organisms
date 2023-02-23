function dotPlotFracToSize (frac, assay) {
    if (assay == "gene_expression")
        return 2 + 18 * Math.sqrt(frac);
    else
        return 2 + 60 * Math.sqrt(frac);
}

function dotPlotSizeToFrac(markersize, assay) {
    if (assay == "gene_expression")
        frac = (markersize - 2) / 18.0;
    else
        frac = (markersize - 2) / 60.0;
    return frac * frac;
}

function computeMarkerSizeOvertime(nc) {
    let ms;
    if (nc == 0) {
        ms = 2;
    } else if (nc < 5) {
        ms = 8;
    } else if (nc < 40) {
        ms = 13;
    } else {
        ms = 20;
    }
    return ms;
}

function getDomains(axiss, invert=false) {
    let nPlots = axiss.length;
    let ntot = axiss.reduce((acc, a) => acc + a.length, 0);
    let domains = [];
    let acc;

    if (invert)
        acc = 1;
    else
        acc = 0;
    for (let k=0; k < nPlots; k++) {
        let frac = 1.0 * axiss[k].length / ntot;
        let domain;
        if (invert) {
            domain = [acc - frac, acc];
            acc -= frac;
        } else {
            domain = [acc, acc + frac];
            acc += frac;
        }
        domains.push(domain);
    }
    return domains
}

function getPseudocount(featureType) {
    switch (featureType) {
        case "gene_expression":
            return 0.5;
        case "chromatin_accessibility":
            return 0.01;
        default:
            return 1;
    }
}

function getTickTexts(axis, feature_type, geneIds) {
    let tickText, featureUrl;
    let ticktexts = [];
    for (let i=0; i < axis.length; i++) {
        const feature = axis[i];
        if (feature_type == "gene_expression") {
            const geneId = geneIds[feature];
            if (geneId !== "") {
                featureUrl = feature;
                if (geneId.startsWith('MGI')) {
                    featureUrl = 'http://www.informatics.jax.org/marker/'+geneId;
                } else {
                    featureUrl = 'https://www.genecards.org/cgi-bin/carddisp.pl?gene='+geneId;
                }
                tickText = '<a href="'+featureUrl+'">'+feature+'</a> - <span><b>O</b></span>';
            } else {
                tickText = feature;
            }
        } else if (feature_type == "chromatin_accessibility") {
            const featureParts = feature.split('-');
            const chrom = featureParts[0];
            const start = featureParts[1];
            const end = featureParts[2];

            // UCSC used chr4:445-566, colon is urlencode-d by the browser automatically
            featureUrl = 'https://genome.ucsc.edu/cgi-bin/hgTracks?db=mm10&position=' + chrom + ':' + start + '-' + end;
            tickText = '<a href="'+featureUrl+'">'+feature+'</a>';
        } else {
            tickText = "" + feature;
        }
        ticktexts.push(tickText);
    }
    return ticktexts;
}

export {
    dotPlotSizeToFrac, dotPlotFracToSize, computeMarkerSizeOvertime,
    getDomains, getPseudocount, getTickTexts,
};
