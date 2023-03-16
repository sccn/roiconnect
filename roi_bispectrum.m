function EEG = roi_bispectrum(EEG, filt)

    data = EEG.roi.source_roi_data;
    [nchan, ndat, nepo] = size(data);
    
    nROI = EEG.roi.nROI;
    segleng = ndat;
    segshift = floor(ndat/2);
    epleng = ndat;
    
    fres = EEG.pnts/2;
    frqs = sfreqs(fres, EEG.srate);
    freqinds_low = [find(frqs==mean(filt.low)) find(frqs==(mean(filt.high)-mean(filt.low)))];
    freqinds_up = [find(frqs==mean(filt.low)) find(frqs==mean(filt.high))];

    for iroi = 1:nROI
        for jroi = iroi:nROI
            X = data([iroi jroi],:,:);
            [bs_up,~] = data2bs_event(X(:,:)', segleng, segshift, epleng, freqinds_up);
        end
    end
end

