function measure = rm_components(measure, nPCA)
    % only keep the first PC
    % measure has the size (n_freq, nROI*nPCA, nROI*nPCA)
    if nPCA > 1
        measure = measure(:, 1:nPCA:end, 1:nPCA:end);
    end
end