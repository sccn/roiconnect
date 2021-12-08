% roi_activity - call roi_connectivity_process to compute
%                                connectivity between ROIs
% Usage:
%  EEG = roi_activity(EEG, 'key', 'val', ...);
%
% Inputs:
%  EEG - EEGLAB dataset
%
% Required inputs:
%  'leadfield'   - [string or struct] leadfield structure (Fieldtrip)
%                  or file containing Fieldtrip or Brainstrom leadfield
%                  matrix.
%  'sourcemodel' - [string] source model file also containing Atlas info.
% 
% Optional inputs:
%  'sourcemodel2mni' - [9x float] homogeneous transformation matrix to convert
%                  sourcemodel to MNI space.
%  'sourcemodelatlas' - [string] name of Atlas to use (must be contained
%                       in Atlas field of the sourcemodel file.
%  'morder'    - [interger] Autoregressive model order (default is 20)
%  'nPCA'      - [interger] Number of PCA component for each ROI. Each ROI
%                is made of many voxel. Instead of averaging their activity,
%                this function takes the x first PCA components, then use
%                these to compute connectivity (default is 3)
%  'naccu'     - [interger] For bootstrap, number of accumulation. Default is 
%                none.
%  'eloretareg' - [float] regularization term for eLoreta. Default is 0.05.
%  'trgc'      - ['on'|'off'] compute time-reverse Granger Causality. Default
%                is 'on'.
%  'crossspec' - ['on'|'off'] compute cross-spectrum from which coherence can
%                be derived. Default is 'on'.
%  'roiactivity'  - ['on'|'off'] compute ROI activity. Default is on. If
%                you just need voxel activity, you can set this option to
%                'off'.
%  'exportvoxact' - ['on'|'off'] export voxel activity in EEG.roi.source_voxel_data
%                These array are huge, so the default is 'off'.
%
% Output:
%  EEG - EEGLAB dataset with field 'roi' containing connectivity info.
%  source_voxel_data - voxel data activity (voxels x times x trials).
%                      Usually several Gb in size.
%
% Author: Stefan Haufe and Arnaud Delorme
%
% Example: call pop_roi_activity instead because it will
% compute the leadfield matrix automatically using Dipfit information.

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

function [EEG, source_voxel_data] = roi_activity(EEG, varargin)

if nargin < 2
    help roi_activity;
    return
end

% decode input parameters
% -----------------------
g = finputcheck(varargin, { ...
    'leadfield'  {'struct' 'string'}  {{} {}}             '';
    'headmodel'   'string'  { }             ''; % sometimes useful when loading volume to see which voxels are inside/outside
    'sourcemodel' 'string'  { }             '';
    'sourcemodel2mni' 'real'    { }         [];
    'sourcemodelatlas' 'string'    { }      '';
    'modelparams'   'cell'    { }          { 0.05 };
    'model'          'string'    { 'eLoretaFieldtrip' 'lcmvFieldtrip' 'eLoreta' 'lcmv' } 'eLoreta';
    'nPCA'        'integer' { }              3;
    'downsample'  'integer' { }              1;
    'roiactivity' 'string' { 'on' 'off' }    'on';
    'exportvoxact' 'string' { 'on' 'off' }   'off';
    'outputdir'   'string'  { }              '' }, 'roi_activity');
if ischar(g), error(g); end
if isempty(g.leadfield), error('Leadfield is mandatory parameter'); end

%%% Creating result folder
if ~isempty(g.outputdir)
    mkdir(fullfile( g.outputdir, 'data'));
end

% Cortex mesh or volume
% ---------------------
[~,~,ext] = fileparts(g.sourcemodel);

if strcmpi(ext, '.head')
    [~, grid, labels, strlabels ] = load_afni_atlas(g.sourcemodel, g.headmodel, g.sourcemodel2mni, g.downsample);
    uniqueROIs = unique(labels);
    nROI = length(uniqueROIs);
    cortex.Atlas(1).Name = g.sourcemodelatlas;
    for iROI = 1:nROI
        indVertices = find(labels == uniqueROIs(iROI));
        cortex.Atlas(1).Scouts(iROI).Label    = strlabels{iROI};
        cortex.Atlas(1).Scouts(iROI).Vertices = indVertices;
    end
    cortex.Vertices = grid;
else
    cortex = load(g.sourcemodel);
    if isfield(cortex, 'Faces')
        % make brainstorm coordinate system consistent with MNI coordinates for
        % plotting (in terms of axis directions)
        disp('Brainstorm cortex mesh detected - transforming to MNI coordinates');
        tf = traditionaldipfit(g.sourcemodel2mni);
        pos      = tf*[cortex.Vertices ones(size(cortex.Vertices,1),1)]';
        pos      = pos';
        cortex.Vertices = pos(:,1:3);
    elseif isfield(cortex, 'cortex') && isfield(cortex, 'atlas')
        hm = cortex;
        clear cortex;
        % align with MNI coordinates
        if ~isempty(g.sourcemodel2mni)
            tf = traditionaldipfit(g.sourcemodel2mni);
            pos      = tf*[hm.cortex.vertices ones(size(hm.cortex.vertices,1),1)]';
            pos      = pos';
        else
            pos = hm.cortex.vertices;
        end
        cortex.Vertices = pos(:,1:3);
        cortex.Faces = hm.cortex.faces;

        % make Alejandro Atlas definition compatible with Brainstrom one
        nROI = length(hm.atlas.label);
        cortex.Atlas(1).Name = g.sourcemodelatlas;
        for iROI = 1:nROI
            indVertices = find(hm.atlas.colorTable == iROI);
            cortex.Atlas(1).Scouts(iROI).Label    = hm.atlas.label{iROI};
            cortex.Atlas(1).Scouts(iROI).Vertices = indVertices;
        end
    elseif isnumeric(cortex) && mod(size(cortex,1),3) == 0 && size(cortex,2) == 6
        % NFT matrix
        cortextmp = cortex;
        clear cortex;
        cortex.Vertices = cortextmp(:,1:3);
        cortex.Atlas(1).Name = g.sourcemodelatlas;
    elseif ~isfield(cortex, 'Vertices')
        % code below is functional to load a mesh
        % However, need to align with an Atlas
        % This can be achieve with Fieldtrip functions
        sourcemodelOriOld = ft_read_headshape(fullfile(ftPath, 'template', 'sourcemodel', 'cortex_20484.surf.gii'));

        error('Unknown mesh format')
    end
end

% Select Atlas
% ------------
found = false;
for iAtlas = 1:length(cortex.Atlas)
    if strcmpi(cortex.Atlas(iAtlas).Name, g.sourcemodelatlas)
        cortex.Atlas = cortex.Atlas(iAtlas);
        found = true;
        break
    end
end
if ~found
    error('Atlas not found');
end

% leadfield matrix (Brainstorm or Fieldtrip)
% ------------------------------------------
if ~isstruct(g.leadfield)
    leadfield = load(g.leadfield, '-mat');
else
    leadfield = g.leadfield;
end
if isstruct(leadfield) && isfield(leadfield, 'roiconnectleadfield') 
    leadfield = leadfield.roiconnectleadfield;
elseif isstruct(leadfield) && isfield(leadfield, 'Gain') 
    % brainstorm
    % make format compatible with Stefan's routines
    leadfield = permute(reshape(leadfield.Gain, [], 3, nvox), [1 3 2]);
elseif isstruct(leadfield) && isfield(leadfield, 'leadfield') 
    % fieldtrip
    oldLeadfield = leadfield;
    leadfield.gain = reshape( [ leadfield.leadfield{:} ], [length(leadfield.label) 3 length(leadfield.leadfield)]);
    leadfield.gain = permute(leadfield.gain, [1 3 2]);
    leadfield = leadfield.gain;
elseif isfield(leadfield, 'LFM')
    % NFT
    leadfield = leadfield.LFM;
else
    disp('Warning: unknown leadfield matrix format, assuming array of gain values');
end

nvox = size(cortex.Vertices, 1);
nvox2 = size(leadfield,2);
if ~isequal(nvox, nvox2)
    error('There must be the same number of vertices/voxels in the leadfield and source model');
end
if ~isequal(size(leadfield,1), EEG.nbchan)
    error('There must be the same number of channels in the leadfield and in the dataset');
end

% use frequency resolution of 0.5 Hz
fres = EEG.srate;

% from the MVGC toolbox, compute frequencies in Hz for a
frqs = sfreqs(fres, EEG.srate);

%% source reconstruction
if strcmpi(g.model, 'eLoreta')
    % common average reference transform
    H = eye(EEG.nbchan) - ones(EEG.nbchan) ./ EEG.nbchan;

    % apply to data and leadfield
    EEG.data = reshape(H*EEG.data(:, :), EEG.nbchan, EEG.pnts, EEG.trials);
    leadfield = reshape(H*leadfield(:, :), EEG.nbchan, nvox, 3);

    % eLORETA inverse projection kernel
    disp('Computing eLoreta...');
    P_eloreta = mkfilt_eloreta_v2(leadfield, g.modelparams{:});
    
    % project to source space
    source_voxel_data = reshape(EEG.data(:, :)'*P_eloreta(:, :), EEG.pnts*EEG.trials, nvox, 3);
elseif strcmpi(g.model, 'LCMV')
    C = cov(EEG.data(:, :)');
    if length(g.modelparams) == 1
        lcmv_reg = g.modelparams{1};
    end
    alpha = lcmv_reg*trace(C)/length(C);
    Cr = C + alpha*eye(EEG.nbchan);
    [~, P_eloreta] = lcmv(Cr, leadfield, struct('alpha', 0, 'onedim', 0));
    source_voxel_data = reshape(EEG.data(:, :)'*P_eloreta(:, :), EEG.pnts*EEG.trials, nvox, 3);
    source_voxel_data = 10^3*source_voxel_data; % the units are nA*m
else
    % transform the data to continuous so we can get an estimate for each sample
    EEG2 = EEG;
    EEG2.data = EEG2.data(:,:);
    EEG2.pnts = size(EEG2.data,2);
    EEG2.trials = 1;
    EEG2 = eeg_checkset(EEG2);
    dataPre = eeglab2fieldtrip(EEG2, 'preprocessing', 'dipfit');  
    
    % prepare data
    cfg = [];
    cfg.channel = {'all', '-EOG1'};
    cfg.reref = 'yes';
    cfg.refchannel = {'all', '-EOG1'};
    dataPre = ft_preprocessing(cfg, dataPre);

    % load head model and prepare leadfield matrix
    vol = load('-mat', g.headmodel);

    % source reconstruction
    cfg             = [];
    if lower(g.model(1)) == 'e'
        cfg.method      = 'eLoreta';
    else
        cfg.method      = 'lcmv';
    end
    try
        cfg.(g.sourcemethod) = struct(g.modelparams{:});
    catch, end
    cfg.sourcemodel = oldLeadfield;
    cfg.headmodel   = vol.vol;
    cfg.keeptrials  = 'yes';
    source          = ft_sourceanalysis(cfg, dataPre);  % compute the source
    
    % reformat for ROI analysis below
    source_voxel_data = reshape([ source.avg.mom{:} ], 3, size(source.avg.mom{1},2), length(source.avg.mom));
    source_voxel_data = permute(source_voxel_data, [2 3 1]);
end
    
% number of ROIs in the Desikan-Killiany Atlas
nROI  = length(cortex.Atlas.Scouts);
nPCAs = zeros(1, nROI);

% ROI labels
labels = {cortex.Atlas.Scouts.Label};

% keep only the first nPCA strongest components for each ROI
if strcmpi(g.roiactivity, 'on')
    disp('Computing ROI activity...');
    source_roi_data = [];
    
    for iROI = 1:nROI
        ind_roi = cortex.Atlas.Scouts(iROI).Vertices;
        [source_roi_power(iROI), source_roi_power_norm(iROI)] = roi_getpower(source_voxel_data, ind_roi);
        [source_roi_data_tmp, nPCAs(iROI)] = roi_getact(source_voxel_data, ind_roi, g.nPCA);
        source_roi_data = cat(2, source_roi_data, source_roi_data_tmp);
    end

    % version with nPCA components
    source_roi_data = permute(reshape(source_roi_data, EEG.pnts, EEG.trials, []), [3 1 2]);
else
    source_roi_data = [];
    source_roi_power_norm = [];
end
disp('Done');

% Output paramters
EEG.roi.cortex    = cortex;
EEG.roi.atlas     = cortex.Atlas.Scouts;
if strcmpi(g.exportvoxact, 'on')
    EEG.roi.source_voxel_data     = source_voxel_data; % large (takes lots of RAM)
end
EEG.roi.source_roi_data       = single(source_roi_data);
EEG.roi.source_roi_power_norm = source_roi_power_norm; % used for cross-sprectum
EEG.roi.freqs     = frqs;
EEG.roi.nPCA      = g.nPCA;
EEG.roi.nROI      = nROI;
EEG.roi.atlas     = cortex.Atlas;
EEG.roi.srate     = EEG.srate;
EEG.roi.leadfield = g.leadfield;
EEG.roi.headmodel = g.headmodel;
EEG.roi.parameters = varargin;
if exist('P_eloreta', 'var')
    EEG.roi.P_eloreta = single(P_eloreta);
end



return;







nPCAs = zeros(1, length(nROI));
nPCAs(:) = 3;
beg_inds = cumsum([1 nPCAs(1:end-1)]);
end_inds = cumsum([nPCAs]);

for iroi = 1:nROI
  PCA_inds{iroi} = beg_inds(iroi):end_inds(iroi);
end

%% ROI analysis using nPCA components per ROI

% test TRGC between ROIs (that is, pairs of nPCA-dimensional spaces)
% compute indices
% new: 
inds = {}; ninds = 0;
for iroi = 1:nROI
  for jroi = (iroi+1):nROI
    inds{ninds+1} = {PCA_inds{iroi}, PCA_inds{jroi}};    
%     inds{ninds+2} = {PCA_inds{jroi}, PCA_inds{iroi}};  
    ninds = ninds + 1;
  end
end


% compute TRGC
% note, this is not yet the net score i->j - j->i , 
% only the difference forward-backward


% TRGC_ = data2sctrgc(source_roi_data, fres, morder, 0, nbootstrap, [], inds);
% new FC measures: MIC and MIM (Ewald and Nolte) are multivariate undirected (symmetric) FC
% metrics. They represent the first eigenvalue and the sum of eigenvalues
% of a maximization of iCOH. 
conn_mult = data2sctrgcmim(source_roi_data, fres, morder, 0, nbootstrap, [], inds, {'TRGC', 'GC', 'MIM', 'MIC'});

% compute cross-spectrum
conn_uni = data2spwctrgc2(source_roi_data, fres, morder, 0, nbootstrap, [], {'CS'});

%% TRGC plotting

% calculation of net TRGC scores (i->j minus j->i), recommended procedure
% TRGCnet = TRGC_(:, 1:2:end)-TRGC_(:, 2:2:end);
% new way to compute net scores
TRGCnet = conn_mult.TRGC(:, :, 1) - conn_mult.TRGC(:, :, 2);

% create a ROI x ROI connectivity matrix, if needed
% TRGCmat(f, ii, jj) is net TRGC from jj to ii
TRGC = [];
iinds = 0;
for iroi = 1:(nROI)
  for jroi = (iroi+1):(nROI)
    iinds = iinds + 1;
    TRGC(:, iroi, jroi) = -TRGCnet(:, iinds);
    TRGC(:, jroi, iroi) = TRGCnet(:, iinds);
  end
end

% plot overall TRGC spectrum, sanity check for numerical instabilities
figure; plot(frqs, TRGC(:, :))

% TRGC per frequency band
TRGCband = [];
for iba = 1:length(bands)
  TRGCband(iba, :, :) = sum(TRGC(frq_inds{iba}, :, :));
end

% plot alpha-band net TRGC for each region
roi_ = sum(TRGCband(2, :, :), 3);
allplots_cortex_BS(cortex_highres, roi_, [-max(abs(roi_)) max(abs(roi_))], cm17, 'TRGC', smooth_cortex, ...
  [fig_folder sbj '/TRGC_net_alpha']);

figure; imagesc(squeeze(TRGCband(2, :, :)))

% use left inferior parietal region as sender and plot TRGC to all other
% areas
% positive means seed sends more than it receives
roi_ = squeeze(TRGCband(2, :, 15));
allplots_cortex_BS(cortex_highres, roi_, [-max(abs(roi_)) max(abs(roi_))], cm17, 'TRGC', smooth_cortex, ...
  [fig_folder sbj '/TRGC_net_alpha_left_inferior_parietal_sender']);


%% iCOH plotting

% note: as an alternative to iCOH, the MIM/MIC scores could be brought into
% matrix form and be plotted.

% compute (absolute value of) the imaginary part of coherence
iCOH_ = abs(imag(cs2coh(conn_uni.CS)));

clear iCOH
for iroi = 1:nROI
  for jroi = 1:nROI
    iCOH(:, iroi, jroi) = mean(mean(iCOH_(:, PCA_inds{iroi}, PCA_inds{jroi}), 2), 3);
  end
end


% iCOH per frequency band
iCOHband = [];
for iba = 1:length(bands)
  iCOHband(iba, :, :) = sum(iCOH(frq_inds{iba}, :, :));
end

%plot alpha-band net iCOH for each region
roi_ = mean(iCOHband(2, :, :), 3);
allplots_cortex_BS(cortex_highres, roi_, [min(abs(roi_)) max(abs(roi_))], cm17a, '|iCOH|', smooth_cortex, ...
  [fig_folder sbj '/iCOH_net_alpha']);

figure; imagesc(squeeze(iCOHband(2, :, :)))

% use right paracentral region as seed and plot iCOH with all other
% areas
% positive means more interaction
roi_ = squeeze(iCOHband(2, :, 34));
allplots_cortex_BS(cortex_highres, roi_, [min(abs(roi_)) max(abs(roi_))], cm17a, '|iCOH|', smooth_cortex, ...
  [fig_folder sbj '/iCOH_net_alpha_right_paracentral_seed']);


%%

% create a ROI x ROI MIC connectivity matrix
MIC = [];
iinds = 0;
for iroi = 1:(nROI)
  for jroi = (iroi+1):(nROI)
    iinds = iinds + 1;
    MIC(:, iroi, jroi) = conn_mult.MIC(:, iinds);
    MIC(:, jroi, iroi) = conn_mult.MIC(:, iinds);
  end
end

% MIC per frequency band
MICband = [];
for iba = 1:length(bands)
  MICband(iba, :, :) = sum(MIC(frq_inds{iba}, :, :));
end

%plot alpha-band net MIC for each region
roi_ = mean(MICband(2, :, :), 3);
allplots_cortex_BS(cortex_highres, roi_, [min(abs(roi_)) max(abs(roi_))], cm17a, 'MIC', smooth_cortex, ...
  [fig_folder sbj '/MIC_net_alpha']);

% use right paracentral region as seed and plot MIC with all other
% areas
% positive means more interaction
roi_ = squeeze(MICband(2, :, 34));
allplots_cortex_BS(cortex_highres, roi_, [min(abs(roi_)) max(abs(roi_))], cm17a, 'MIC', smooth_cortex, ...
  [fig_folder sbj '/MIC_net_alpha_right_paracentral_seed']);



%% power plotting

% plot regionwise broadband power on smoothed surface
voxel_ = nan*ones(size(cortex.Vertices, 1), 1);
voxel_(ind_cortex) = 10*log10(reshape(sum(sum(source_voxel_data.^2), 3), [], 1));
allplots_cortex_BS(cortex_highres, voxel_(in_normal_to_high), [min(voxel_) max(voxel_)], ...
  cm17a, 'power [dB]', smooth_cortex, [fig_folder sbj '/voxel_broadband_voxel_power_db_smooth']);


% % plot regionwise broadband power
% source_roi_power_norm_dB = 10*log10(source_roi_power_norm);
% allplots_cortex_BS(cortex_highres, source_roi_power_norm_dB, [min(source_roi_power_norm_dB) max(source_roi_power_norm_dB)], ...
%   cm17a, 'power [dB]', 0, ...
%   [fig_folder sbj '/roi_broadband_power_db']);

% plot regionwise broadband power on smoothed surface
source_roi_power_norm_dB = 10*log10(source_roi_power_norm);
allplots_cortex_BS(cortex_highres, source_roi_power_norm_dB, [min(source_roi_power_norm_dB) max(source_roi_power_norm_dB)], ...
  cm17a, 'power [dB]', smooth_cortex, [fig_folder sbj '/roi_broadband_voxel_power_db_smooth']);


clear PS
PS_ = cs2psd(conn_uni.CS);
for iroi = 1:nROI
  PS(:, iroi) = sum(PS_(:, PCA_inds{iroi}), 2);
end

% another way to plot region-wise broadband power 
% todo: need to check for the proper scales after the Fourier transform
roi_ = 10*log10(sum(PS)./nvoxel_per_roi)';
allplots_cortex_BS(cortex_highres, roi_, [min(roi_) max(roi_)], cm17a, 'power [dB]', smooth_cortex, ...
  [fig_folder sbj '/roi_broadband_power_db2']);
      

% power per frequency band
PSband = [];
for iba = 1:length(bands)
  PSband(iba, :) = sum(PS(frq_inds{iba}, :));
end

% 
roi_ = 10*log10(PSdBband(2, :)./nvoxel_per_roi);
allplots_cortex_BS(cortex_highres, roi_, [min(roi_) max(roi_)], cm17a, 'power [dB]', smooth_cortex, ...
  [fig_folder sbj '/power_dB_alpha']);

