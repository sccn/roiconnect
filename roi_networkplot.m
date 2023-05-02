% roi_networkplot() - multiple calls to plotconnectivity
%
% Usage:
%     [imgFile, txtFile] = roi_networkplot(EEG, networkdefs, measure, 'key', val);
%
% Input:
%     EEG         - EEGLAB structure with connectivity pre-calculated
%                   using pop_roi_connectplot
%     networkdefs - [string|struct] network definition file (same input as
%                   roi_definenetwork()). Alternatively network structure
%                   (second output of roi_definenetwork())
%                   attempts to use EEG.roi.atlas.networks
%     measure     - [string] measure to plot. Same as pop_roi_connectplot()
%                   Alternatively, cell array of connectivity matrices one
%                   per network.
%
% Optional inputs:
%  'freqrange'  - [min max] frequency range in Hz. Default is to plot broadband power.
%  'subplots'   - ['on'|'off'] create subplots (when more than one plot).
%                 Default is 'off'.
%  'exportmean' - ['on'|'off'] export results as text files (when filename defined)
%                 Default is 'on'.
%  'exportstd'  - ['on'|'off'] export results as text files (when filename defined)
%                 Default is 'on'.
%  'title'      - ['string'] figure title. Default none.
%  'plotmode'   - ['2D'|'3D'|'both'] figure plotting mode. Default is '2D'
%  'filename'   - ['string'] base file name (without extension). This is
%                 used to save a variety of files with different postfix. 
%                 Default is empty and no file are saved.
%
% Outputs:
%     EEG     - modified EEG structure
%     imgFile - Image file list
%     txtFile - Text file list
%
% Example:
% % Compute ROI connectivity
% EEG = pop_dipfit_settings( EEG ); % select boundary element model
% EEG = pop_leadfield(EEG, 'sourcemodel','dipfit/LORETA-Talairach-BAs.mat','sourcemodel2mni',[],'downsample',1);
% EEG = pop_roi_activity(EEG, 'resample','on','regepochs','on','leadfield',EEG.dipfit.sourcemodel,...
%       'model','LCMV','modelparams',{0.05},'atlas','LORETA-Talairach-BAs','nPCA',3);
% EEG = pop_roi_connect(EEG, 'methods', { 'CS'});
%
% % Define a single network and plot it
% DNM = { '25L' '25R' '32L' '32R' '33L' '33R' '23L' '23R' '31L' '31R' '39L' '39R' }';
% [EEG, net] = roi_definenetwork(EEG, table(DNM));
% roi_networkplot(EEG, net, 'crossspecimag', 'freqrange', [8 12]);
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


% or structure containing the 
%                   fields:
%                    - chanlocs (from EEG.chanlocs)
%                    - loreta_P (from EEG.roi.P_eloreta)
%                    - roiStruct (structure from Scouts EEG.roi.atlas.Scouts(:))  
%                    - networks(x).name (name of network x)
%                    - networks(x).ROI_inds (indices of ROIs for network x)

function [EEG,imgFileName,txtFileName,measures] = roi_networkplot(EEG, networks, measure, varargin)

if nargin < 2
    help roi_networkplot;
    return
end

[g, addopts] = finputcheck(varargin, { ...
    'subplots'    'string'    {'on' 'off'}    'off';
    'exportmean'  'string'    {'on' 'off'}    'on';
    'exportstd'   'string'    {'on' 'off'}    'off';
    'title'       'string'    {}              '';
    'netstat'     'string'    {'on' 'off'}    'off';
    'addrois'     ''          {}              [];
    'columns'     'integer'   {}              [];
    'limits'      'float'     {}              [];
    'plotmode'    'string'    {'2D' '3D' 'both' 'none' }  '2D';
    'plotopt'     'cell'      {}              {};
    'filename'    'string'    {}              '';
    'threshold'   'float'     {}              0.1;
    'precomputed'    'struct'    {}      struct([]);
    }, 'roi_networkplot', 'ignore');
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
if ~isfield(EEG, 'roi')
    error('"roi" field not present in EEG structure, use pop_roi_activity to compute ROI activity');
end
if ~isfield(EEG.roi, 'atlas')
    error('"roi.atlas" field not present in EEG structure, use pop_leadfield to use a source model with an atlas');
end
% decode network param
if isempty(networks)
    if ~isfield(EEG.roi.atlas, 'networks')
        error('"roi.atlas.networks" field not found in EEG structure, use roi_definenetworks to define networks');
    end
    networks = EEG.roi.atlas.networks;
end

% legacy code reading reading saved network
if ischar(networks)
    try
        networks = load('-mat', networks);
    catch
    end
end

% get value of matrix based on measure for the frequency of interest
% ------------------------------------------------------------------
fprintf('Thresold of %s (all connectivity values below the threshold are removed)\n', num2str(g.threshold));
if ischar(measure) % measure contains the name of the measure
    matrix = pop_roi_connectplot(EEG, 'measure', measure, 'noplot', 'on', addopts{:});
elseif iscell(measure) % measure contains connectivity matrices
    if length(measure) ~= length(networks)
        error('When a cell array, "measure" should have as many element as networks');
    end
    if ~isempty(addopts)
        error('Unknown option "%s"', addopts{1});
    end
    matrix = measure;
else
    matrix = measure;
end

% get network and convert if necessary
% ------------------------------------
if isstruct(networks) % in case network is already converted from a table to a structure
    [EEG,~,matrix] = roi_definenetwork(EEG, [], 'addrois', g.addrois, 'connectmat', matrix, 'ignoremissing', 'on');
else
    [EEG,networks,matrix] = roi_definenetwork(EEG, networks, 'addrois', g.addrois, 'connectmat', matrix, 'ignoremissing', 'on');
end

if strcmpi(g.subplots, 'on')
    if isempty(g.columns)
        ncol = ceil(sqrt(length(networks)));
    else
        ncol = g.columns;
    end
    nrow = ceil(length(networks)/ncol);
    figure('position', [100 100 350*ncol 350*nrow], 'paperpositionmode', 'auto');
end

imgFileName = {};
txtFileName = {};
roiStruct =  EEG.roi.atlas.Scouts(:);
if length(g.threshold) < length(networks)
    g.threshold = g.threshold*ones(1,length(networks));
end
sumVals = {};
for iNet = 1:length(networks)
    if length(networks(iNet).ROI_inds) < 2
        error('Cannot plot network %s: you need at least two brain areas to make a network', length(networks(iNet).ROI_inds));
    end
    % create structure containing connectivity for the network of interest
    if iscell(matrix)
        networkMatSubj = matrix{iNet};
        if size(networkMatSubj,1) ~= size(networks(iNet).ROI_inds,1)
            try
                networkMatSubj = networkMatSubj(networks(iNet).ROI_inds, networks(iNet).ROI_inds, :);
            catch
                error('When a cell array, "measure" should contain matrices which have the same number of elements as the corresponding network');
            end
        end
    else
        networkMatSubj = matrix(networks(iNet).ROI_inds, networks(iNet).ROI_inds, :);
    end
    if 0
        % sum
        sumVals{iNet} = nansum(nansum(networkMatSubj,1),2);
    elseif 0
        % max
        sumVals{iNet} = nanmax(nanmax(networkMatSubj,[],1),[],2);
    else
        % median
        tmp = networkMatSubj;
        for i1 = 1:size(tmp,1) 
            for i2 = 1:size(tmp,3)
                tmp(i1,i1,i2) = NaN;
            end
        end
        tmp = reshape(tmp, size(tmp,1)*size(tmp,2), size(tmp,3));
        sumVals{iNet} = nanmedian(tmp, 1);
    end
    networkMat    = mean(networkMatSubj,3);
    networkMatStd = std(networkMatSubj,[],3);

    if ~strcmpi(g.plotmode, 'none')
        if strcmpi(g.subplots, 'on')
            subplot(nrow, ncol, iNet)
        else
            figure('position', [100 100 400 700], 'paperpositionmode', 'auto');
        end
    end

    labels = { roiStruct(networks(iNet).ROI_inds).Label };
    tmpTitle = networks(iNet).name;
    tmpTitle(tmpTitle == '_') = ' ';
    tmpTitle = [ tmpTitle ' ' g.title ];
    
    % 2-D plot
    if strcmpi(g.plotmode, '2D') || strcmpi(g.plotmode, 'both')
        plotconnectivity(networkMat(:,:), 'labels', labels, 'axis', gca, 'threshold', g.threshold(iNet), 'limits', g.limits, g.plotopt{:});
        h = title(tmpTitle, 'interpreter', 'none');
        pos = get(h, 'position');
        set(h, 'position', pos + [0 0.1 0]);
    
        % save individual plots
        if ~strcmpi(g.subplots, 'on') && ~isempty(g.filename)
            set(h, 'fontsize', 16, 'fontweight', 'bold');
            tmpFileName = [ g.filename '_' networks(iNet).name '.jpg' ];
            imgFileName{end+1} = tmpFileName;
            print('-djpeg', tmpFileName );
            close
        end
    end
    
    % 3-D plot
    if strcmpi(g.plotmode, '3D') || strcmpi(g.plotmode, 'both')
        options = {'brainmovieopt' { 'moviename' '' } g.plotopt{:} };
        if ~strcmpi(g.subplots, 'on') && ~isempty(g.filename)
            tmpFileName = [ g.filename '_' networks(iNet).name '_3d' ];
            options = { options{:} 'filename'  tmpFileName };
            imgFileName{end+1} = [ tmpFileName '.xhtml' ];
        end
        roi_plotbrainmovie(networkMat(:,:), 'labels', labels, 'threshold', g.threshold(iNet), options{:});
    end
    
    % save text
    if strcmpi(g.exportmean, 'on') && ~isempty(g.filename)
        tmptable = array2table(networkMat(:,:), 'variablenames', labels, 'rownames', labels);
        tmpFileName = [ g.filename '_' networks(iNet).name '_mean.txt' ];
        txtFileName{end+1} = tmpFileName;
        writetable(tmptable, tmpFileName,'WriteRowNames', true );
    end
    if strcmpi(g.exportstd, 'on') && ~isempty(g.filename)
        tmptable = array2table(networkMatStd(:,:), 'variablenames', labels, 'rownames', labels);
        tmpFileName = [ g.filename '_' networks(iNet).name '_std.txt' ];
        txtFileName{end+1} = tmpFileName;
        writetable(tmptable, tmpFileName,'WriteRowNames', true );
    end

    if nargin > 2
        tmpTitle(tmpTitle == ' ') = '_';
        tmpTitle(tmpTitle == '.') = '_';
        tmpTitle(tmpTitle == '-') = '_';
        measures.(tmpTitle).mean   = networkMat;
        measures.(tmpTitle).labels = labels;
    end
    
end

% perform some stats
if length(networks) > 1 && length(sumVals) == 2 && length(sumVals{1}) == length(sumVals{2})
    [H,P,CI,STATS] = ttest(sumVals{1}, sumVals{2});
    if P< 0.05
        fprintf(2, 'Network %s effect size %1.2f p-value: %1.3f\n', networks(end).name, STATS.tstat, P);
    else
        fprintf('Network %s effect size %1.2f p-value: %1.3f\n', networks(end).name, STATS.tstat, P);
    end
    %[ef,~,pvals] = statcond(sumVals(1:2), 'mode', 'perm');
    %fprintf('Network %s effect size %1.2f p-value: %1.3f\n', networks(end).name, ef, pvals);
end

% if strcmpi(g.subplots, 'on') && ~isempty(g.title)
%     h = textsc(g.title, 'title');
%     set(h, 'fontsize', 16, 'fontweight', 'bold');
% end

if strcmpi(g.subplots, 'on') && ~isempty(g.filename)
    print('-djpeg', g.filename);
    close
end
