function measure = get_connect_mat( measureOri, nROI, signVal)
    % create a ROI x ROI connectivity matrix, if needed
    % TRGCmat(f, ii, jj) is net TRGC from jj to ii
    measure = [];
    iinds = 0;
    for iroi = 1:nROI
        for jroi = (iroi+1):nROI
            iinds = iinds + 1;
            measure(:, iroi, jroi) = signVal * measureOri(:, iinds);
            measure(:, jroi, iroi) = measureOri(:, iinds);
        end
    end
end