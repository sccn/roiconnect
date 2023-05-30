% pop_roi_connect - call roi_connect to connectivity between ROIs
%
% Usage:
%  EEG = pop_roi_connect(EEG, 'key', 'val', ...);
%
% Inputs:
%  EEG - EEGLAB dataset containing ROI activity
%
% Optional inputs:
%  'morder'   - [integer]  Order of autoregressive model. Default is 20.
%  'nepochs'  - [integer] number of data epoch. This is useful when
%               comparing conditions. if not enough epochs can be extracted
%               an error is returned. If there are too many, the first ones
%               are selected (selecting the first epochs ensure they are mostly
%               contiguous and that the correlation between them is similar 
%               accross conditions).
%  'naccu'          - [integer]  Number of accumulation for stats. Default is 0.
%  'methods'        - [cell] Cell of strings corresponding to methods.
%                       'CS'    : Cross spectrum
%                       'COH'   : Coherence
%                       'GC'    : Granger Causality
%                       'TRGC'  : Time-reversed Granger Causality
%                       'wPLI'  : Weighted Phase Lag Index
%                       'PDC'   : Partial directed coherence
%                       'TRPDC' : Time-reversed partial directed coherence
%                       'DTF'   : Directed transfer entropy
%                       'TRDTF' : Time-reversed directed transfer entropy
%                       'MIM'   : Multivariate Interaction Measure for each ROI
%                       'MIC'   : Maximized Imaginary Coherency for each ROI
%                       'PAC'   : Phase-amplitude coupling between ROIs
%  'snippet'        - ['on'|off]  Option to compute connectivity over snippets. Default is 'off'. 
%  'snip_length'    - ['on'|'off']  Length of the snippets. Default is 60 seconds.
%  'fcsave_format'  - ['mean_snips'|'all_snips']  Option to save mean over snippets 
%                     (shape: 101,68,68) or all snippets (shape: n_snips,101,68,68). Default is 'mean_snips.'
%  'freqresolution' - [integer] Desired frequency resolution (in number of frequencies). 
%                     If specified, the signal is zero padded accordingly. Default is 0 (means no padding).
%  'fcomb'          - [struct] Frequency combination for which PAC is computed (in Hz). Must have fields 'low' and 
%                     'high' with fcomb.low < fcomb.high. For example, fcomb.low = 10 and fcomb.high = 50 if single 
%                     frequencies are used. fcomb.low = [4 8] and fcomb.high = [48 50] if frequency bands are used 
%                     (might take a long time to compute so use with caution). Default is {} (this will cause an error).
%  'bs_outopts'     - [integer] Option which bispectral tensors should be stored in EEG.roi.PAC. Default is 1.
%                          1 - store all tensors: b_orig, b_anti, b_orig_norm, b_anti_norm
%                          2 - only store: b_orig, b_anti
%                          3 - only store: b_orig_norm, b_anti_norm
%                          4 - only store: b_orig, b_orig_norm
%                          5 - only store: b_anti, b_anti_norm
%  'roi_selection'  - [cell array of integers] Cell array of ROI indices {1, 2, 3, ...} indicating for which regions (ROIs) connectivity should be computed. 
%                     Default is all (set to EEG.roi.nROI).
%
% Output:
%  EEG - EEGLAB dataset with field 'roi' containing connectivity info.
%
% Note: Optional inputs to roi_connectivity_process() are also accepted.
%
% Author: Arnaud Delorme, UCSD, 2019
%
% Example
%   p = fileparts(which('eeglab')); % path
%   EEG = pop_roi_connect(EEG, 'headmodel', ...
%   EEG.dipfit.hdmfile, 'elec2mni', EEG.dipfit.coord_transform, ...
%   'sourcemodel', fullfile(p, 'functions', 'supportfiles', ...
%   'head_modelColin27_5003_Standard-10-5-Cap339.mat'), 'sourcemodel2mni', ...
%   [0 -26.6046230000 -46 0.1234625600 0 -1.5707963000 1000 1000 1000]);
%
% Use pop_roi_connectivity(EEG) to conectivity

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

% TO DO - Arno
% - Centralize reading head mesh and Atlas (there might be a function in
% Fieldtrip to do that) ft_read_volume ft_read_mesh
% - Make compatible with all Fieldtrip and FSL Atlases
% - Downsampling of Atlas - check bug submitted to Fieldtrip
% - Plot inside(blue) vs outside(red) voxels for source volume

function [EEG,com] = pop_roi_connect(EEG, varargin)

com = '';
if nargin < 1
    help pop_roi_connect;
    return
end

if ~isfield(EEG(1), 'roi') || ~isfield(EEG(1).roi, 'source_roi_data')
    error('Cannot find ROI data - ROI data first');
end

if nargin < 2

    rowg = [0.1 0.6 1 0.2];
    % uigeom = { 1 1 rowg rowg 1 rowg rowg [0.1 0.6 0.9 0.3] 1 rowg 1 [0.5 1 0.35 0.5] [0.5 1 0.35 0.5] [0.5 1 0.35 0.5] [1] [0.9 1.2 1] };
    uigeom = { [1] [1.2 1] [1.2 1] [1.2 1] [1.2 1] [1.2 1] [1.2 1] [1] [0.2 1 0.35 0.8] [0.2 1 0.35 0.8] };
    uilist = { { 'style' 'text' 'string' 'Select connectivity measures' 'fontweight' 'bold' } ...
               { 'style' 'checkbox' 'string' 'Cross-spectrum'               'tag' 'cs' 'value' 1  } {} ...
               { 'style' 'checkbox' 'string' 'Coherence'                    'tag' 'coh' 'value' 0  } ... 
               { 'style' 'checkbox' 'string' 'Weighted Phase Lag Index'   'tag' 'wpli' 'value' 0  } ...
               { 'style' 'checkbox' 'string' 'Granger Causality (GC)'          'tag' 'gc' 'value' 0   } ...
               { 'style' 'checkbox' 'string' 'Time-reversed GC'                'tag' 'trgc' 'value' 0   } ...
               { 'style' 'checkbox' 'string' 'Partial Directed Coherence (PDC)' 'tag' 'pdc' 'value' 0   } ...
               { 'style' 'checkbox' 'string' 'Time-reversed PDC'               'tag' 'trpdc' 'value' 0   } ...
               { 'style' 'checkbox' 'string' 'Directed Transfer Entropy (DTF)'  'tag' 'dtf' 'value' 0   } ...
               { 'style' 'checkbox' 'string' 'Time-reversed DTF'               'tag' 'trdtf' 'value' 0   } ...
               { 'style' 'checkbox' 'string' 'Multivariate Interaction Measure'       'tag' 'mim' 'value' 0  } ...
               { 'style' 'checkbox' 'string' 'Maximized Imaginary Coherency'       'tag' 'mic' 'value' 0   } ...
               {} ...
               {} { 'style' 'text' 'string' 'Autoregressive model order'   } { 'style' 'edit' 'string' '20' 'tag' 'morder' }  {} ...
               {} { 'style' 'text' 'string' 'Bootstrap if any (n)'         } { 'style' 'edit' 'string' '' 'tag' 'naccu2' }  {} };
                ...
    [result,~,~,out] = inputgui('geometry', uigeom, 'uilist', uilist, 'helpcom', 'pophelp(''pop_roi_connect'')', 'title', 'pop_roiconnect - connectivity');
    if isempty(result), return, end
                      
    % check we have the same naccu
    methods = {};
    if out.cs,    methods = [ methods { 'CS' } ]; end
    if out.coh,   methods = [ methods { 'COH' } ]; end
    if out.gc  ,  methods = [ methods { 'GC' } ]; end
    if out.trgc,  methods = [ methods { 'TRGC' } ]; end
    if out.wpli,  methods = [ methods { 'wPLI' } ]; end
    if out.pdc  , methods = [ methods { 'PDC' } ]; end
    if out.trpdc, methods = [ methods { 'TRPDC' } ]; end
    if out.dtf  , methods = [ methods { 'DTF' } ]; end
    if out.trdtf, methods = [ methods { 'TRDTF' } ]; end
    if out.mim  , methods = [ methods { 'MIM' } ]; end
    if out.mic,   methods = [ methods { 'MIC' } ]; end
    options = {  ...
        'morder' str2num(out.morder) ...
        'naccu' str2num(out.naccu2) ...
        'methods' methods };
else
    options = varargin;
end

% decode input parameters
% -----------------------
g = finputcheck(options, ...
    { 'morder'         'integer'  { }                           20;
      'naccu'          'integer'  { }                           0;
      'methods'        'cell'     { }                           { };
      'snippet'        'string'   { 'on', 'off' }               'off';
      'nepochs'        'real'                {}                 [ ];
      'snip_length'    'integer'  { }                           60; 
      'fcsave_format'  'string'   { 'mean_snips', 'all_snips'}  'mean_snips';
      'freqresolution' 'integer'  { }                           0; 
      'fcomb'          'struct'   { }                           struct; 
      'bs_outopts'     'integer'  { }                           1; 
      'roi_selection'  'cell'     { }                           {EEG.roi.nROI} }, 'pop_roi_connect');
if ischar(g), error(g); end

% process multiple datasets
% -------------------------
if length(EEG) > 1
    eeglab_options;
    if nargin < 2
        if option_storedisk
            res = questdlg2( [ 'Data files on disk will be automatically overwritten.' 10 ...
                                'Are you sure you want to proceed with this operation?' ], ...
                            'Confirmation', 'Cancel', 'Proceed', 'Proceed');
            switch lower(res)
                case 'cancel', return;
                case 'proceed'
            end
        end
    end
    % find common datasets accross subjects
    % and limit the number of epochs
    allsubjects = { EEG.subject };
    for index = 1:length(allsubjects)
        optionsSubj = options;
        allinds = strmatch(allsubjects{index}, allsubjects, 'exact');

        % datasets from the same subjects need to have the same number of epochs
        TMPEEG = EEG(allinds);
        if TMPEEG(1).trials > 1
            optionsSubj = [ optionsSubj { 'nepochs' min([TMPEEG.trials]) }];
        end
        if nargin < 2
            [ TMPEEG, com ] = eeg_eval( 'pop_roi_connect', TMPEEG, 'warning', 'off', 'params', optionsSubj );
        else
            [ TMPEEG, com ] = eeg_eval( 'pop_roi_connect', TMPEEG, 'params', optionsSubj );
        end
        EEG = eeg_store(EEG, TMPEEG, allinds);
    end
    return
end

if ~isempty(g.nepochs) 
    if EEG.trials < g.nepochs
        error('Not enough data epochs')
    elseif EEG.trials > g.nepochs
        % epochSelect = shuffle(1:EEG.trials);
        % better select contiguous epochs
        EEG = pop_select(EEG, 'trial', 1:g.nepochs);
    end
end

% compute connectivity over snippets
n_conn_metrics = length(options{2}); % number of connectivity metrics
conn_matrices_snips = {};
if strcmpi(g.snippet, 'on')
    snippet_length = g.snip_length; % seconds
    snip_eps = snippet_length/(size(EEG.data,2)/EEG.srate); % n epochs in snippet
    nsnips = floor(EEG.trials/snip_eps);
    if nsnips < 1
        error('Snippet length cannot exceed data length.')
    end
    diff = (EEG.trials * EEG.pnts/EEG.srate) - (nsnips * EEG.pnts/EEG.srate * snip_eps);
    if diff ~= 0
        warning(strcat(int2str(diff), ' seconds are thrown away.'));
    end

    source_roi_data_save = EEG.roi.source_roi_data;
    for isnip = 1:nsnips
        roi_snip = source_roi_data_save(:,:,(isnip-1)* snip_eps + 1 : (isnip-1)* snip_eps + snip_eps); % cut source data into snippets
        EEG.roi.source_roi_data = single(roi_snip);
        EEG = roi_connect(EEG, 'morder', g.morder, 'naccu', g.naccu, 'methods', g.methods,'freqresolution',g.freqresolution); % compute connectivity over one snippet
        for fc = 1:n_conn_metrics 
            fc_name = options{2}{fc};
            fc_matrix = EEG.roi.(fc_name);
            conn_matrices_snips{isnip,fc} = fc_matrix; % store each connectivity metric for each snippet in separate structure
        end
    end
    
    % compute mean over connectivity of each snippet
    for fc = 1:n_conn_metrics
        fc_name = options{2}{fc};
        [first_dim, second_dim, third_dim] = size(conn_matrices_snips{1,fc});

        conn_cell = conn_matrices_snips(:,fc); % store all matrices of one metric in a cell
        mat = cell2mat(conn_cell);
        reshaped = reshape(mat, first_dim, nsnips, second_dim, third_dim);
        reshaped = squeeze(permute(reshaped, [2,1,3,4]));
        if strcmpi(g.fcsave_format, 'all_snips')
            EEG.roi.(fc_name) = reshaped;
        else
            if nsnips > 1
                mean_conn = squeeze(mean(reshaped, 1)); 
            else
                mean_conn = reshaped;
            end
            EEG.roi.(fc_name) = mean_conn; % store mean connectivity in EEG struct
        end
    end
else
    EEG = roi_connect(EEG, 'morder', g.morder, 'naccu', g.naccu, 'methods', g.methods,'freqresolution',g.freqresolution);
end

if ~isempty(intersect(g.methods, {'PAC'}))
    EEG = roi_pac(EEG, g.fcomb, g.bs_outopts, g.roi_selection);
end

if nargout > 1
    com = sprintf( 'EEG = pop_roi_connect(EEG, %s);', vararg2str( options ));
end

