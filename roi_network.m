% roi_network() - compute connectivity between ROIs
%
% Usage:
%  EEG = roi_network(EEG, 'key', 'val', ...); 
%
% Inputs:
%  EEG - EEGLAB dataset with ROI activity computed
%
% Optional inputs (choose at least one):
%  'nfft'        - [integer] FFT padding. Default is twice the sampling rate.
%  'freqdb'      - ['on'|'off'] compute spctral activity in dB. Default is 'on'.
%  'freqrange'   - [cell] frequency ranges. Default is { [4 6] [ 8 12] [18 22] }
%                  for theta (4 to 6 Hz), alpha and beta.
%  'processfreq' - [struct of func] how to process spectral data. Default is
%                   processfreq.theta = @(x)x(:,1);
%                   processfreq.alpha = @(x)x(:,2);
%                   processfreq.beta  = @(x)x(:,3);
% 'processconnect' - [struct of func] how to process connectivity data.
%                   Default is (the diverder is the number of non zero values)
%                   processconnect.theta = @(x)sum(sum(x(:,:,1)))/((size(x,1).^2)-size(x,1));
%                   processconnect.alpha = @(x)sum(sum(x(:,:,2)))/((size(x,1).^2)-size(x,1));
%                   processconnect.beta  = @(x)sum(sum(x(:,:,3)))/((size(x,1).^2)-size(x,1));
% 'plotmode'       - ['2D'|'3D'|'both'] plot in 2-D, 3-D or both. Default
%                   is 2D.
%
% Output:
%   EEG - EEG structure with EEG.roi field updated and now containing
%         connectivity information.
%   results - result with the fields defined as input.

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

function [EEG,results,loretaFile,imgFileName,txtFileName] = roi_network(EEG, varargin)

if EEG.trials == 1 % fast call
    opt = struct(varargin{:});
    if ~isfield(opt, 'networkfile')
        error('No field networkfile given as input');
    end
    
    % other paratemers
    if ~isfield(opt, 'roilist'),       opt.roilist = []; end 
    if ~isfield(opt, 'nfft'),          opt.nfft = EEG.srate*2; end
    if ~isfield(opt, 'freqdb'),        opt.freqdb = 1; end
    if ~isfield(opt, 'freqrange')      opt.freqrange = { [4 6] [ 8 12] [18 22] }; end
    if ~isfield(opt, 'processfreq')    opt.processfreq = []; end
    if ~isfield(opt, 'processconnect') opt.processconnect = [];  end
    if ~isfield(opt, 'plotnetworkfile') opt.plotnetworkfile = '';  end
else
    opt = finputcheck( varargin, { ...
        'networkfile'    ''          {}      '';
        'nfft'           'integer'   {}      EEG.srate*2;
        'freqdb'         'integer'   {}      1;
        'measureoutput'  'string'    {'on' 'off'}      'off';
        'roilist'        'integer'   {}      []; 
        'freqrange'      'cell'      {}      { [4 6] [ 8 12] [18 22] };
        'freqname'       'cell'      {}      { 'theta' 'alpha' 'beta' };
        'processfreq'    ''          {}      [];
        'precomputed'    'struct'    {}      struct([]);
        'leadfield'      'string'    {}      '';
        'sourcemodel'    'string'    {}      '';
        'plotnetworkfile' ''         {}      '';
        'plotmode'        'string'    {'2D' '3D' 'both' 'off' }  '2D';
        'plotloretafile'  ''         {}      '';
        'loretalimits'    ''         {}      [];
        'processconnect'  ''         {}      [] }, 'roi_network');
end     
if ischar(opt), error(opt); end
if isempty(opt.processfreq)
    for iFreq = 1:length(opt.freqrange)
        opt.processfreq.(opt.freqname{iFreq}) = @(x)x(:,iFreq);
    end
end
if isempty(opt.processconnect)
    for iFreq = 1:length(opt.freqrange)
        opt.processconnect.(opt.freqname{iFreq}) = @(x)sum(sum(x(:,:,iFreq)))/((size(x,1).^2)-size(x,1)); % the diverder is the number of non zero values
    end
end

if isempty(opt.networkfile)
    dipfitdefs;
    EEG = pop_dipfit_settings( EEG, 'hdmfile', template_models(2).hdmfile,'coordformat','MNI', ...
        'mrifile', template_models(2).mrifile,'chanfile',template_models(2).chanfile,...
        'coord_transform',[0 0 0 0 0 -1.5708 1 1 1] ,'chansel',[1:EEG.nbchan] );
    chans = { 'FP1','FP2','F3','F4','C3','C4','P3','P4','O1','O2','F7','F8','T3','T4','T5','T6','FZ','CZ','PZ' };
    if ischar(opt.leadfield)
        tmpLeadfield = load('-mat', opt.leadfield);
        if isfield(tmpLeadfield, 'label') && ~isequal(upper({EEG.chanlocs.labels}), upper(tmpLeadfield.label))
            error('Electrode name inconsistency');
        end
    end
    options = { ...
        'headmodel' template_models(2).hdmfile ...
        'leadfield' opt.leadfield ...
        'sourcemodel' opt.sourcemodel ... % 'sourcemodel' '/data/matlab/eeglab/plugins/roiconnect/LORETA-Talairach-BAs.mat' ...
        'sourcemodel2mni' [] ...
        'sourcemodelatlas' 'BrainDx' ... %    'sourcemodelatlas' 'LORETA-Talairach-BAs' ...
        'downsample' 1 ...
        'nPCA' 1 ...
        'model' 'eLoreta' ...
        'roiactivity' 'on', ...
        'exportvoxact' 'on'
        };
    %             'trgc' 'off' ...
    %             'crossspec' 'off' ...
    %             'morder' 20 ...
    EEG = pop_roi_activity(EEG, options{:});
    source_voxel_data = EEG.roi.source_voxel_data;
    if isempty(opt.roilist)
        opt.roilist = 1:EEG.roi.nROI; % list of ROI necessary to compute connectivity
    end
    tmp = load('-mat', 'supportfiles\BrainDx_sourcemodel.mat');
    loreta_ROIS = tmp.Atlas.Scouts;
    loreta_Networks = [];
else
    if ischar(opt.networkfile)
        opt.networkfile = load('-mat', opt.networkfile);
    end
    loreta_P        = opt.networkfile.loreta_P;
    loreta_Networks = opt.networkfile.loreta_Networks;
    loreta_ROIS     = opt.networkfile.loreta_ROIS;
    if isempty(opt.roilist)
        opt.roilist = unique([loreta_Networks.ROI_inds]); % list of ROI necessary to compute connectivity
    end

    % project to source space
    source_voxel_data = reshape(EEG.data(:, :)'*loreta_P(:, :), size(EEG.data,2)*size(EEG.data,3), size(loreta_P,2), 3);
end

% Computing spectrum
% ALSO IMPLEMENT USING ROI_ACTIVITY
sz = size(source_voxel_data);
tmpdata = reshape(source_voxel_data, sz(1), sz(2)*sz(3)); % THIS IS MOSTLY WRONG HERE AS EPOCHS ARE CONCATENATED
source_voxel_spec = pwelch(tmpdata, EEG.srate, EEG.srate/2, opt.nfft, EEG.srate); % assuming 1 second of data
source_voxel_spec = reshape(source_voxel_spec, size(source_voxel_spec,1), sz(2), sz(3));
source_voxel_spec = mean(source_voxel_spec(2:size(source_voxel_spec,1),:,:,:),length(sz)); % frequency selection 2 to 31 (1Hz to 30Hz)
freqs  = linspace(0, EEG.srate/2, floor(opt.nfft/2)+1);
freqs  = freqs(2:end); % remove DC (match the output of PSD)

% Plot loreta file
if ~isempty(opt.plotloretafile)
    loretaFile = opt.plotloretafile;
    options = { 'freqrange', opt.freqrange, 'limits', opt.loretalimits, 'saveasfile',  opt.plotloretafile, 'precomputed', opt.precomputed };
    if strcmpi(opt.plotmode, 'off') options = [ options { 'noplot' 'on' } ]; end
    loretaMeasures = roi_sourceplot(freqs, source_voxel_spec', opt.sourcemodel, options{:});
end

% Compute ROI activity
for ind_roi = opt.roilist
    % data used for connectivity analysis
    spatiallyFilteredDataTmp = roi_getact( source_voxel_data, loreta_ROIS(ind_roi).Vertices, 1, 0); % Warning no zscore here; also PCA=1 is too low
    spatiallyFilteredSpecTmp = roi_getact( source_voxel_spec, loreta_ROIS(ind_roi).Vertices, 1, 0);
    if ind_roi == 1
        spatiallyFilteredData = zeros(max(opt.roilist), length(spatiallyFilteredDataTmp));
        spatiallyFilteredSpec = zeros(max(opt.roilist), length(spatiallyFilteredSpecTmp));
    end
    spatiallyFilteredData(ind_roi,:) = spatiallyFilteredDataTmp;
    spatiallyFilteredSpec(ind_roi,:) = spatiallyFilteredSpecTmp;
end
loretaSpec = spatiallyFilteredSpec';

% select frequency bands
for iSpec = 1:length(opt.freqrange)
    freqRangeTmp = intersect( find(freqs >= opt.freqrange{iSpec}(1)), find(freqs <= opt.freqrange{iSpec}(2)) );
    loretaSpecSelect(:,iSpec) = mean(abs(loretaSpec(freqRangeTmp,:)).^2,1); % mean power in frequency range
    if opt.freqdb
        loretaSpecSelect(:,iSpec) = 10*log10(abs(loretaSpecSelect(:,iSpec)).^2);
    end
end

% compute metric of interest
processfreqFields = fieldnames(opt.processfreq);
for iProcess = 1:length(processfreqFields)
    results.(['loreta_regions_' processfreqFields{iProcess}]) = feval(opt.processfreq.(processfreqFields{iProcess}), loretaSpecSelect);
end

% compute cross-spectral density for each network
% -----------------------------------------------
if ~isempty(opt.processconnect)
    for iNet = 1:length(loreta_Networks)
        if 1
            % ALSO IMPLEMENT USING ROI_CONNECT
            [restmp,connectSpecSelect{iNet}] = roi_csnetworkact( spatiallyFilteredData, loreta_Networks(iNet).ROI_inds, 'nfft', opt.nfft, 'postprocess', opt.processconnect, 'freqranges', opt.freqrange);
            % copy results
            fields = fieldnames(restmp);
            for iField = 1:length(fields)
                meanField   = [ loreta_Networks(iNet).name '_' fields{iField} ];
                detailField = [ loreta_Networks(iNet).name '_' fields{iField} '_details' ];
                results.(meanField) = restmp.(fields{iField});
                results.(detailField) = connectSpecSelect{iNet}(:,:,iField);
                
                % reuse data
                if isfield(opt.precomputed, meanField)
                    restmp.(fields{iField}) = opt.precomputed.(meanField);
                end           
                if isfield(opt.precomputed, detailField)
                    connectSpecSelect{iNet}(:,:,iField) = opt.precomputed.(detailField);
                end           
            end
        else
            networkData = spatiallyFilteredData(loreta_Networks(iNet).ROI_inds,:);
            S = cpsd_welch(networkData,size(networkData,2),0,g.measure.nfft);
            [nchan, nchan, nfreq] = size(S);
            
            % imaginary part of cross-spectral density
            % ----------------------------------------
            absiCOH = S;
            for ifreq = 1:nfreq
                absiCOH(:, :, ifreq) = squeeze(S(:, :, ifreq)) ./ sqrt(diag(squeeze(S(:, :, ifreq)))*diag(squeeze(S(:, :, ifreq)))');
            end
            absiCOH = abs(imag(absiCOH));
            
            % frequency selection
            % -------------------
            connectSpecSelect = zeros(size(absiCOH,1), size(absiCOH,2), length(opt.freqrange));
            for iSpec = 1:length(g.measure.freqrange)
                freqRangeTmp = intersect( find(freqs >= opt.freqrange{iSpec}(1)), find(freqs <= opt.freqrange{iSpec}(2)) );
                connectSpecSelect(:,:,iSpec) = mean(absiCOH(:,:,freqRangeTmp),3); % mean power in frequency range
            end
            
            connectprocessFields = fieldnames(opt.processconnect);
            for iProcess = 1:length(connectprocessFields)
                results.([ loreta_Networks(iNet).name '_' connectprocessFields{iProcess} ]) = feval(opt.processconnect.(connectprocessFields{iProcess}), connectSpecSelect);
            end
        end
    end
    if ~isempty(loreta_Networks) && ~isempty(opt.plotnetworkfile) && ~strcmpi(opt.plotmode, 'off')
        imgFileName = {};
        txtFileName = {};
        for iField = 1:length(fields)
            connectTmp = cellfun(@(x)x(:,:,iField), connectSpecSelect, 'uniformoutput', false);
            [imgFileNameTmp,txtFileNameTmp] = roi_plotnetwork(opt.networkfile, connectTmp, 'title', fields{iField}, 'filename' ,[opt.plotnetworkfile '_' fields{iField} ], 'plotmode', opt.plotmode);
            imgFileName = [ imgFileName imgFileNameTmp ];
            txtFileName = [ txtFileName txtFileNameTmp ];
        end
    else
        imgFileName = {};
        txtFileName = {};
    end
    
end

if strcmpi(opt.measureoutput, 'on')
    out.measures = [];
    fields = fieldnames(results);
    for iField = 1:length(fields)
        out.measures.(fields{iField}).mean = results.(fields{iField});
    end
    fields = fieldnames(loretaMeasures);
    for iField = 1:length(fields)
        out.measures.(fields{iField}).mean = loretaMeasures.(fields{iField});
    end
    results = out;
end