function EEG = roi_bispectrum(EEG, filt)
    % filt: frequenzpaar
    % test: filt.low < filt.high

    data = EEG.roi.source_roi_data;
    [nchan, ndat, nepo] = size(data);
    
    nROI = EEG.roi.nROI;
    segleng = ndat;
    segshift = floor(ndat/2);
    epleng = ndat;

    params = [];
    params.nROI = nROI;
    params.segleng = seglent;
    params.segshift = segshift;
    params.epleng = epleng;
    
%     [b_orig, b_anti, b_orig_norm,b_anti_norm] = pac_bispec(data,fs,filt);
    [b_orig, b_anti, b_orig_norm,b_anti_norm] = pac_bispec(data, params, filt); % suggestion, unpack params in pac_bispec

end

