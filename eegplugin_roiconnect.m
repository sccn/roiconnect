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
% Authors: Arnaud Delorme, Stefan Haufe, Guido Nolte

% Copyright (C) Arnaud Delorme, arnodelorme@gmail.com
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
% THE POSSIBILITY OF SUCH DAMAGE.

function vers = eegplugin_roiconnect(fig, trystrs, catchstrs)

vers = 'roiconnect1.0';
if nargin < 3
    error('eegplugin_roiconnect requires 3 arguments');
end

% use global variable to assess ROI status
% ----------------------------------------
try
    curroiFlag = evalin('base', 'unique(cellfun(@(x)getfield(x, ''eeglab_using_roi''), { ALLEEG(:).roi }));');
    if length(curroiFlag) ~= 1 && curroiFlag
        disp('Warning: different settings for using ROIs detected in different datasets');
    end
catch
    curroiFlag = 0;
end

% find tools menu
% ---------------
tool_m = findobj(fig, 'tag', 'tools');

% command to check that the '.source' is present in the EEG structure
% -------------------------------------------------------------------
cb_toggle      = [ 'try, [ALLEEG, LASTCOM] = pop_toggleroiandica(ALLEEG); EEG = ALLEEG(CURRENTSET);' catchstrs.add_to_hist '; eeglab(''redraw'');' ]; 
cb_act         = [ 'try, [EEG, LASTCOM] = pop_roi_activity(EEG);' catchstrs.new_and_hist ]; 
cb_connect     = [ 'try, [EEG, LASTCOM] = pop_roi_connect(EEG);'  catchstrs.new_and_hist ]; 
cb_plot        = [ 'try, LASTCOM = pop_roi_connectplot(EEG);'     catchstrs.add_to_hist  ]; 

roi_m = uimenu( tool_m, 'label', 'ROI connectivity analysis', 'userdata', 'startup:off;study:on');
if curroiFlag
    uimenu( roi_m, 'Label', 'Have EEGLAB use ICA instead of ROIs', 'CallBack', cb_toggle, 'tag', 'toggleroi', 'userdata', 'startup:off;study:on');
else
    uimenu( roi_m, 'Label', 'Have EEGLAB use ROIs instead of ICA', 'CallBack', cb_toggle, 'tag', 'toggleroi', 'userdata', 'startup:off;study:on');
end
uimenu( roi_m, 'Label', 'Compute ROI activity', 'CallBack', cb_act, 'userdata', 'startup:off;study:on', 'separator', 'on');
uimenu( roi_m, 'Label', 'Compute ROI connectivity', 'CallBack', cb_connect, 'userdata', 'startup:off;study:on');
uimenu( roi_m, 'Label', 'Plot ROI connectivity', 'CallBack', cb_plot);


