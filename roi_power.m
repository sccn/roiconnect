% compute ROI power

function [source_roi_power, source_roi_power_norm] = roi_power(source_voxel_data, ind_roi)

data_  = source_voxel_data(:, ind_roi, :);
source_roi_power = sum(var(data_(:, :)))';
source_roi_power_norm = source_roi_power/length(ind_roi);

