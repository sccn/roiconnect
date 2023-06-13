% Function to test if the ROI selection works by comparing the values of the respective data arrays. 
% Let 'n_roisel' be the number of selected ROIs. Used in 'test_roi_selection'.
%
% Inputs:
%   roi_selection - [cell array of integers] Cell array of ROI indices {1, 2, 3, ...} indicating for which regions (ROIs) connectivity was computed. 
%   conn1         - (nfreq, nROI, nROI) connectivity tensor, original data array without ROI selection
%   conn2         - (nfreq, n_roisel, n_roisel) connectivity tensor, data array with ROI selection (connectivity was computed for a selected
%                   number of ROIs)

function cmp_roi_selection(roi_selection, conn1, conn2)
% generate all possible combinations of ROI combinations

    [p,q] = meshgrid(cell2mat(roi_selection), cell2mat(roi_selection));
    roi_pairs = [p(:) q(:)]; 
    
    roi_inds = 1:1:length(roi_selection);
    [p,q] = meshgrid(roi_inds, roi_inds);
    roi_inds_pairs = [p(:) q(:)]; 
    
    for ipair = 1:size(roi_pairs, 1)
        cmp1 = conn1(roi_pairs(ipair, 1), roi_pairs(ipair, 2));
        cmp2 = conn2(roi_inds_pairs(ipair, 1), roi_inds_pairs(ipair, 2));
        if ~isequal(cmp1, cmp2)
            disp(cmp1)
            disp(cmp2)
            fprintf('Values for the ROI combination [%d %d] are not exactly the same.\n', roi_pairs(ipair, :))
        end
    end
end