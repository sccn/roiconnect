% compute ROI activity

function source_roi_data = roi_act(source_voxel_data, ind_roi, nPCA)

if nargin < 3
    nPCA = 1;
end

% optional z-scoring, this makes the PCA independent of the power in each
% voxel, and favors to find components that are consistently expressed in
% many voxels rather than only in a few voxels with strong power (which
% may leak from a neighboring region)
data_  = source_voxel_data(:, ind_roi, :);
data_(:, :) = zscore(data_(:, :));
if nPCA == 1
    [source_roi_data, ~, ~] = svds(double(data_(:, :)), nPCA); 
else
    % old code
    [data_, ~, ~] = svd(data_(:, :), 'econ'); % WARNING SHOULD USE SVDS for SPEED
                                               % however, sometime polarity
                                               % inverted compared to svds
                                               % Should not have an
                                               % incidence but keeping the
                                               % code for now
    source_roi_data = data_(:, 1:nPCA);
end
