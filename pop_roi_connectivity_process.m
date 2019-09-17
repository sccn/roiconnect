% (C) Stefan Haufe 2018
% Modification to use the standar Desikan-Kiliany atlas

% 'morder' model order for the spectral/connectivity analysis

function S = pop_roi_connectivity_process(EEG, varargin)

if nargin < 1
    help pop_roi_connectivity_process;
    return
end

% Source model and ROI
%   - Select file
%   - Select atlas in file
%   - Transformation to align source and leadfield matrix (plot both on top of each other)
% Select leadfield matrix
%   - Use leadfield Matrix from DIPFIT
%   - Use custom leadfield matrix
% Model order
% Bootstrap

S = roi_connectivity_process(EEG, 'leadfield', 'leadfield_tmp.mat', 'cortexfile', cortexfile);


if nargin < 1
    roi(1).label = 'Use Desikan-Kilianny ROI in Colin27 template';
    roi(1).file  = 'head_modelColin27_5003_Standard-10-5-Cap339.mat';
    roi(1).align = [0 -26.6046230000 -46 0.1234625600 0 -1.5707963000 1000 1000 1000];
    roi(1).enable = 'off';
    
    roi(2).label = 'Use Desikan-Kilianny ROI in ICBM152 template (Brainstrom)';
    roi(2).file  = 'tess_cortex_mid_low_2000V.mat';
    roi(2).align = [0 0 0 0 0 -1.5707963000 1 1 1];
    roi(2).enable = 'off';
    
    uigeom = { [1 1] [1 1] 1 [0.3 1.1 1 1] };
    uilist = {{ 'style' 'text' 'string' 'ROI connectivity analysis' 'fontweight' 'bold'} ...
              { 'style' 'popupmenu' 'string' { roi.label } } ...
              {} { 'style' 'text' 'string' '' 'tag' 'strfile' } ...
              {} { 'style' 'text' 'string' '' 'tag' 'transform' } ...
              { 'style' 'popupmenu' 'string' {'Use DIPFIT Leadfield matrix' 'Use custom Leadfield matrix' } 'tag' 'selection' } ...
              {} 
              { 'style' 'text' 'string' 'Model order for AR model' } ...
              { 'style' 'edit' 'string' '20' } ...
              {} ...
              { 'style' 'text' 'string' 'Bootstrap if any (n)' } ...
              { 'style' 'edit' 'string' '' } ...
              {} ...
              { 'style' 'checkbox' 'string' 'Compute TRGC' 'tag' 'trgc' 'value' 1 } ...
              { 'style' 'checkbox' 'string' 'Compute cross-spectrum' 'tag' 'spec' 'value' 1 } ...
              };
              
    [result,~,~,outs] = inputgui('geometry', uigeom, 'uilist', uilist, 'helpcom', 'pophelp(''pop_loadbv'')', ...
        'title', 'Load a Brain Vision Data Exchange format dataset', 'userdata', splot, 'eval', cb_select);
    if isempty(result), return, end


% decode input parameters
% -----------------------
g = finputcheck(varargin, { ...
    'leadfield'  'string'  { }             '';
    'cortexfile' 'string'  { }             '';
    'morder'     'integer' { }              20;
    'eloretareg' 'real'    { }              0.05;
    'nPCA'       'integer' { }              3;
    'naccu'       'integer' { }             1;
    'outputdir'  'string'  { }              'analysis_output' }, 'roi_connectivity_process');
if ischar(g), error(g); end
if isempty(g.leadfield), error('Leadfield is mandatory parameter'); end

% GO TO BRAINSTORM-MASTER3 folder AND START BRAINSTORM
addpath('libs/Daniele_ARMA');
addpath('libs/export_fig');
addpath('libs/haufe');
addpath('libs/mvgc_v1.0');
addpath('libs/mvgc_v1.0/core');
addpath('libs/mvgc_v1.0/stats');
addpath('libs/mvgc_v1.0/utils');
addpath('libs/nolte');
addpath('libs/ssgc_v1.0');
addpath('libs/brainstorm');

leadfieldFlag = 'brainstrom';
leadfieldFlag = 'fieldtrip';

% colormap
load cm17;

%%% Creating result folder
mkdir(fullfile( g.outputdir, 'data'));

% Cortex mesh
% -----------
cortex = load(g.cortexfile);
if isfield(cortex, 'Faces')
    % make brainstorm coordinate system consistent with MNI coordinates for
    % plotting (in terms of axis directions)
    disp('Brainstorm cortex mesh detected - transforming MNI coordinates');
    cortex.Vertices = cortex.Vertices(:, [2 1 3]);
    cortex.Vertices(:, 1) = -cortex.Vertices(:, 1);
elseif isfield(cortex, 'cortex') && isfield(cortex, 'atlas')
    hm = cortex;
    clear cortex;
    % align with MNI coordinates
    tf = traditionaldipfit([0.0000000000 -26.6046230000 -46.0000000000 0.1234625600 0.0000000000 -1.5707963000 1000.0000000000 1000.0000000000 1000.0000000000]);
    pos      = tf*[hm.cortex.vertices ones(size(hm.cortex.vertices,1),1)]';
    pos      = pos';
    cortex.Vertices = pos(:,1:3);
    cortex.Faces = hm.cortex.faces;
    
    % make Alejandro Atlas definition compatible with Brainstrom one
    nROI = length(hm.atlas.label);
    cortex.Atlas(3).Name = 'Desikan-Killiany';
    for iROI = 1:nROI
        indVertices = find(hm.atlas.colorTable == iROI);
        cortex.Atlas(3).Scouts(iROI).Label    = hm.atlas.label{iROI};
        cortex.Atlas(3).Scouts(iROI).Vertices = indVertices;
    end
else
    % code below is functional to load a mesh
    % However, need to align with an Atlas
    % This can be achieve with Fieldtrip functions
    sourcemodelOriOld = ft_read_headshape(fullfile(ftPath, 'template', 'sourcemodel', 'cortex_20484.surf.gii'));
    
    error('Unknown mesh format')
end
 
% leadfield matrix (Brainstorm or Fieldtrip)
% ------------------------------------------
leadfield = load(g.leadfield, '-mat');
if isstruct(leadfield) && isfield(leadfield, 'Gain') % brainstorm
    % make format compatible with Stefan's routines
    leadfield = permute(reshape(leadfield.Gain, [], 3, nvox), [1 3 2]);
elseif isstruct(leadfield) && isfield(leadfield, 'leadfield') % fieldtrip
    leadfield.gain = reshape( [ leadfield.leadfield{:} ], [length(leadfield.label) 3 length(leadfield.leadfield)]);
    leadfield.gain = permute(leadfield.gain, [1 3 2]);
    leadfield = leadfield.gain;
end

nvox = size(cortex.Vertices, 1);
nvox2 = size(leadfield,2);
if ~isequal(nvox, nvox2)
    error('There must be the same number of vertices/voxels in the leadfield and cortex mesh');
end

% use frequency resolution of 0.5 Hz
fres = EEG.srate;

% from the MVGC toolbox, compute frequencies in Hz for a
frqs = sfreqs(fres, EEG.srate);

%% source reconstruction

% common average reference transform
H = eye(EEG.nbchan) - ones(EEG.nbchan) ./ EEG.nbchan;

% apply to data and leadfield
EEG.data = reshape(H*EEG.data(:, :), EEG.nbchan, EEG.pnts, EEG.trials);
leadfield = reshape(H*leadfield(:, :), EEG.nbchan, nvox, 3);

% eLORETA inverse projection kernel
P_eloreta = mkfilt_eloreta_v2(leadfield, g.eloretareg);

% project to source space
source_voxel_data = reshape(EEG.data(:, :)'*P_eloreta(:, :), EEG.pnts*EEG.trials, nvox, 3);

% number of ROIs in the Desikan-Killiany Atlas
nROI = length(cortex.Atlas(3).Scouts);

% ROI labels
labels = {cortex.Atlas(3).Scouts.Label};

% keep only the first nPCA strongest components for each ROI
source_roi_data = [];
for iROI = 1:nROI
    ind_roi = cortex.Atlas(3).Scouts(iROI).Vertices;
    data_  = source_voxel_data(:, ind_roi, :);
    source_roi_power(iROI) = sum(var(data_(:, :)))';
    source_roi_power_norm(iROI) = source_roi_power(iROI)/length(ind_roi);
    
    % optional z-scoring, this makes the PCA independent of the power in each
    % voxel, and favors to find components that are consistently expressed in
    % many voxels rather than only in a few voxels with strong power (which
    % may leak from a neighboring region)
    data_(:, :) = zscore(data_(:, :));
    [data_, ~, ~] = svd(data_(:, :), 'econ');
    source_roi_data(:, :, iROI) = data_(:, 1:g.nPCA);
end

% version with nPCA components
source_roi_data = permute(reshape(source_roi_data, EEG.pnts, EEG.trials, []), [3 1 2]);

%% spectral and connectivity analysis

% to test TRGC between ROIs (that is, pairs of nPCA-dimensional spaces), we
% need to compute these indices
inds = {}; ninds = 0;
for iroi = 1:nROI
    for jroi = (iroi+1):nROI
        inds{ninds+1} = {(iroi-1)*g.nPCA + [1:g.nPCA], (jroi-1)*g.nPCA + [1:g.nPCA]};
        inds{ninds+2} = {(jroi-1)*g.nPCA + [1:g.nPCA], (iroi-1)*g.nPCA + [1:g.nPCA]};
        ninds = ninds + 2;
    end
end

if 0
    TRGC    = [];
    TRGCnet = [];
    TRGCmat = [];
else
    % compute time reversed spectral Granger causality between all pairs of ROIs
    TRGC = data2sctrgc(source_roi_data, fres, g.morder, 0, g.naccu, [], inds);
    
    % calculation of net TRGC scores (i->j minus j->i), recommended procedure
    TRGCnet = TRGC(:, 1:2:end)-TRGC(:, 2:2:end);
    
    % create a ROI x ROI connectivity matrix, if needed
    % TRGCmat(f, ii, jj) is net TRGC from jj to ii
    TRGCmat = [];
    iinds = 0;
    for iroi = 1:nROI
        for jroi = (iroi+1):nROI
            iinds = iinds + 1;
            TRGCmat(:, iroi, jroi) = -TRGCnet(:, iinds);
            TRGCmat(:, jroi, iroi) = TRGCnet(:, iinds);
        end
    end
end

% compute cross-spectrum, takes a while
conn = data2spwctrgc(source_roi_data, fres, g.morder, 0, g.naccu, [], {'CS'});

% Output paramters
S.cortex  = cortex;
S.source_voxel_data     = source_voxel_data;
S.source_roi_data       = source_roi_data;
S.source_roi_power_norm = source_roi_power_norm; % used for cross-sprectum
S.freqs   = frqs;
S.TRGCmat = TRGCmat;
S.CS      = conn.CS;
S.nPCA    = g.nPCA;
S.nROI    = nROI;
S.atlas   = cortex.Atlas(3);
S.srate   = EEG.srate;

