% eegplugin_roiconnect() - Plugin to perform connectivity analysis between
%                          pairs of region of interest (ROIs)
% Usage:
%   >> eegplugin_roiconnect(fig, trystrs, catchstrs);
%
% Inputs:
%   fig        - [integer] eeglab figure.
%   trystrs    - [struct] "try" strings for menu callbacks.
%   catchstrs  - [struct] "catch" strings for menu callbacks.
%
% Authors: Stefan Haufe, Guido Nolte, Arnaud Delorme

function vers = eegplugin_roiconnect(fig, trystrs, catchstrs)

vers = 'roiconnect1.0';
if nargin < 3
    error('eegplugin_roiconnect requires 3 arguments');
end

% find tools menu
% ---------------
tool_m = findobj(fig, 'tag', 'tools');

% command to check that the '.source' is present in the EEG structure
% -------------------------------------------------------------------
cb_process     = [ 'try, [EEG, LASTCOM] = pop_roi_connectivity_process(EEG);' catchstrs.new_and_hist ]; 
cb_plot        = [ 'try, LASTCOM = pop_roi_connectivity_plot(EEG);'           catchstrs.add_to_hist  ]; 

roi_m = uimenu( tool_m, 'label', 'ROI connectivity analysis');
uimenu( roi_m, 'Label', 'Compute ROI connectivity', 'CallBack', cb_process);
uimenu( roi_m, 'Label', 'Plot ROI connectivity', 'CallBack', cb_plot);


