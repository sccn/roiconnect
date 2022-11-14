% roi_plotnetwork() - multiple calls to plotconnectivity
%
% Usage:
%     [imgFile, txtFile] = plotconnectivity(networkfile, array, 'key', val);
%
% Input:
%     networkfile - [string|struct] network file or structure containing the 
%                   fields:
%                    - chanlocs (from EEG.chanlocs)
%                    - loreta_P (from EEG.roi.P_eloreta)
%                    - loreta_ROIS (structure from Scouts EEG.roi.atlas.Scouts(:))  
%                    - loreta_networks(x).name (name of network x)
%                    - loreta_networks(x).ROI_inds (indices of ROIs for network x)
%     array       - [n x n] square array indicating which cell are connected
%
% Optional inputs:
%  'subplots'   - ['on'|'off'] create subplots (when more than one plot).
%                 Default is 'off'.
%  'exporttxt'  - ['on'|'off'] export results as text files (see filename)
%                 Default is 'off'.
%  'title'      - ['string'] figure title. Default none.
%  'plotmode'   - ['2D'|'3D'|'both'] figure plotting mode. Default is '2D'
%  'filename'   - ['string'] base file name (without extension). This is
%                 used to save a variety of files with different postfix. 
%                 Default is empty and no file are saved.
%
% Outputs:
%     imgFile - Image file list
%     txtFile - Text file list
%
% Example:
% EEG = pop_dipfit_settings( EEG ); % select boundary element model
% EEG = pop_leadfield(EEG, 'sourcemodel','dipfit/LORETA-Talairach-BAs.mat','sourcemodel2mni',[],'downsample',1);
% EEG = pop_roi_activity(EEG, 'resample','on','regepochs','on','leadfield',EEG.dipfit.sourcemodel,...
%       'model','LCMV','modelparams',{0.05},'atlas','LORETA-Talairach-BAs','nPCA',3);
% EEG = pop_roi_connect(EEG, 'methods', { 'CS'});
%
% ns.chanlocs = EEG.chanlocs;
% ns.loreta_ROIS = EEG.roi.atlas.Scouts(:);
% ns.loreta_P = EEG.roi.P_eloreta;
% ns.loreta_Networks(1).name     = 'DNM';
% ns.loreta_Networks(1).ROI_inds = [129 127 128];
%
% roi = EEG.roi;
% freqind = 10; % 10th bin (usually 5 Hz)
% connectArray = abs(imag(cs2coh(roi.CS)));
% connectArray = squeeze(mean(mean(reshape({connectArray}, roi.srate+1, roi.nPCA, roi.nROI, roi.nPCA, roi.nROI), 2), 4));
% connectArray = squeeze(mean(connectArray(freqind, :, :),1));
%
% roi_plotnetwork(ns, connectArray)
%
%
% See also: roi_plotbrainmovie() and plotconnectivity()
%
% Author: Arnaud Delorme

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

function [imgFileName,txtFileName] = roi_plotnetwork(networkfile, connectSpecSelect, varargin)

if nargin < 2
    help plotconnectivitymultiple;
    return
end

g = finputcheck(varargin, { ...
    'subplots'    'string'    {'on' 'off'}    'off';
    'exporttxt'   'string'    {'on' 'off'}    'on';
    'title'       'string'    {}              '';
    'plotmode'    'string'    {'2D' '3D' 'both' }  '2D';
    'filename'    'string'    {}              '';
    }, 'roi_network');
if isstr(g)
    error(g);
end
if ~isempty(g.filename) % remove file extension
    [filePath,g.filename] = fileparts(g.filename);
    g.filename = fullfile(filePath, g.filename);
end
if strcmpi(g.subplots, 'on')
    if ~strcmpi(g.plotmode, '2D')
        error('When using subplots, you cannot use 3-D plots');
    end
end
if ~strcmpi(g.plotmode, '2D')
    if ~exist('roi_plotbrainmovie')
        error('You need to install the Brainmovie plugin to plot sources in 3-D');
    end
end
if nargin < 5
    % create a single packet
    singlefig = true;
end

% decode network param
if ischar(networkfile)
    networkfile = load('-mat', networkfile);
end
loreta_Networks = networkfile.loreta_Networks;
loreta_ROIS     = networkfile.loreta_ROIS;

if strcmpi(g.subplots, 'on')
    figure('position', [100 100 1000 700], 'paperpositionmode', 'auto');
    ncol = ceil(sqrt(length(loreta_Networks)));
    nrow = ceil(length(loreta_Networks)/ncol);
end
imgFileName = {};
txtFileName = {};
for iNet = 1:length(loreta_Networks)
    if strcmpi(g.subplots, 'on')
        subplot(nrow, ncol, iNet)
    else
        figure('position', [100 100 400 700], 'paperpositionmode', 'auto');
    end
    
    labels = { loreta_ROIS(loreta_Networks(iNet).ROI_inds).Label };
    tmpTitle = loreta_Networks(iNet).name;
    tmpTitle(tmpTitle == '_') = ' ';
    tmpTitle = [ tmpTitle ' ' g.title ];
    
    % 2-D plot
    if strcmpi(g.plotmode, '2D') || strcmpi(g.plotmode, 'both')
        plotconnectivity(connectSpecSelect{iNet}(:,:), 'labels', labels, 'axis', gca, 'threshold', 0.1);
        h = title(tmpTitle, 'interpreter', 'none');
        pos = get(h, 'position');
        set(h, 'position', pos + [0 0.1 0]);
    
        % save individual plots
        if ~strcmpi(g.subplots, 'on') && ~isempty(g.filename)
            set(h, 'fontsize', 16, 'fontweight', 'bold');
            tmpFileName = [ g.filename '_' loreta_Networks(iNet).name '.jpg' ];
            imgFileName{end+1} = tmpFileName;
            print('-djpeg', tmpFileName );
            close
        end
    end
    
    % 3-D plot
    if strcmpi(g.plotmode, '3D') || strcmpi(g.plotmode, 'both')
        options = {'brainmovieopt' { 'moviename' '' } };
        if ~strcmpi(g.subplots, 'on') && ~isempty(g.filename)
            tmpFileName = [ g.filename '_' loreta_Networks(iNet).name '_3d' ];
            options = { options{:} 'filename'  tmpFileName };
            imgFileName{end+1} = [ tmpFileName '.xhtml' ];
        end
        roi_plotbrainmovie(connectSpecSelect{iNet}(:,:), 'labels', labels, 'threshold', 0.1, options{:});
    end
    
    % save text
    if strcmpi(g.exporttxt, 'on') && ~isempty(g.filename)
        tmptable = array2table(connectSpecSelect{iNet}(:,:), 'variablenames', labels, 'rownames', labels);
        tmpFileName = [ g.filename '_' loreta_Networks(iNet).name '.txt' ];
        txtFileName{end+1} = tmpFileName;
        writetable(tmptable, tmpFileName,'WriteRowNames', true );
    end
    
end

if strcmpi(g.subplots, 'on')
    h = textsc(titl, 'title');
    set(h, 'fontsize', 16, 'fontweight', 'bold');
    print('-djpeg', g.filename);
    close
end