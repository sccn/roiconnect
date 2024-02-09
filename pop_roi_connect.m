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
%                       'aCOH'  : Coherence
%                       'cCOH'  : (complex-valued) coherency
%                       'iCOH'  : absolute value of the imaginary part of coherency
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
%  'firstsnippet'   - ['on'|off]  Only use the first snippet (useful for fast computation). Default is 'off'. 
%  'snip_length'    - ['on'|'off'] Length of the snippets. Default is 60 seconds.
%  'errornosnippet' - ['on'|'off'] Error if snippet too short. Default 'on'.
%  'fcsave_format'  - ['mean_snips'|'all_snips']  Option to save mean over snippets 
%                     (shape: 101,68,68) or all snippets (shape: n_snips,101,68,68). Default is 'mean_snips.'
%  'freqresolution' - [integer] Desired frequency resolution (in number of frequencies). 
%                     If specified, the signal is zero padded accordingly. Default is 0 (means no padding).
%  'fcomb'          - [struct] Frequency combination for which PAC is computed (in Hz). Must have fields 'low' and 
%                     'high' with fcomb.low < fcomb.high. For example, fcomb.low = 10 and fcomb.high = 50 if single 
%                     frequencies are used. fcomb.low = [4 8] and fcomb.high = [48 50] if frequency bands are used 
%                     (might take a long time to compute so use with caution). Default is {} (this will cause an error when PAC is selected).
%  'bs_outopts'     - [integer] Option which bispectral tensors should be stored in EEG.roi.PAC. Default is 1.
%                          1 - store all tensors: b_orig, b_anti, b_orig_norm, b_anti_norm
%                          2 - only store: b_orig, b_anti
%                          3 - only store: b_orig_norm, b_anti_norm
%                          4 - only store: b_orig, b_orig_norm
%                          5 - only store: b_anti, b_anti_norm
%  'roi_selection'  - [cell array of integers] Cell array of ROI indices {1, 2, 3, ...} indicating for which regions (ROIs) connectivity should be computed. 
%                     Default is empty (in this case, connectivity will be computed for all ROIs).
%  'conn_stats'     - ['on'|'off'] Run statistics on connectivity metrics. Default is 'off'.
%  'nshuf'          - [integer] number of shuffles for statistical significance testing. The first shuffle is the true value. Default is 1001. 
%  'freqrange'      - [min max] frequency range in Hz. This is used to compute and plot p-values. Default is to plot broadband power.
%  'poolsize'       - [integer] Number of workers in the parallel pool (check parpool documentation) for parallel computing
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
    uigeom = { [1] [1.2 1] [1.2 1] [1.2 1] [1.2 1] [1.2 1] [1.2 1] [1.2 1] [1] [0.2 1 0.35 0.8] [0.2 1 0.35 0.8] };
    uilist = { { 'style' 'text' 'string' 'Select connectivity measures' 'fontweight' 'bold' } ...
               { 'style' 'checkbox' 'string' 'Cross-spectrum'               'tag' 'cs' 'value' 1  } {} ...
               {'style' 'checkbox' 'string' '(Complex-valued) Coherency'                    'tag' 'ccoh' 'value' 0 } ...
               { 'style' 'checkbox' 'string' 'Coherence'                    'tag' 'acoh' 'value' 0 } ...
               { 'style' 'checkbox' 'string' 'Imaginary Coherency'                    'tag' 'icoh' 'value' 0  } ... 
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
%     if out.coh,   methods = [ methods { 'COH' } ]; end
    if out.ccoh,   methods = [ methods { 'cCOH' } ]; end
    if out.acoh,   methods = [ methods { 'aCOH' } ]; end
    if out.icoh,   methods = [ methods { 'iCOH' } ]; end
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
      'firstsnippet'   'string'   { 'on', 'off' }               'off'; 
      'errornosnippet' 'string'   { 'on', 'off' }               'off'; 
      'nepochs'        'real'                {}               [];
      'snip_length'    'integer'  { }                           60; 
      'fcsave_format'  'string'   { 'mean_snips', 'all_snips'}  'mean_snips';
      'freqresolution' 'integer'  { }                           0; 
      'fcomb'          'struct'   { }                           struct; 
      'bs_outopts'     'integer'  { }                           1; 
      'roi_selection'  'cell'     { }                           { }; 
      'conn_stats'     'string'   { }                           'off'; ...
      'nshuf'          'integer'  { }                           1001; ...
      'poolsize'       'integer'  { }                           1}, 'pop_roi_connect');
if ischar(g), error(g); end

% process multiple datasets
% -------------------------
if length(EEG) > 1
    if nargin < 2
        [ EEG, com ] = eeg_eval( 'pop_roi_connect', EEG, 'warning', 'off', 'params', options );
    else
        [ EEG, com ] = eeg_eval( 'pop_roi_connect', EEG, 'params', options );
    end
    return
end

% compute connectivity over snippets
if strcmpi(g.snippet, 'on') && strcmpi(g.conn_stats, 'off')
    % n_conn_metrics = length(g.methods); 

    snippet_length = g.snip_length; % seconds
    trials = size(EEG.roi.source_roi_data,3);
    pnts   = size(EEG.roi.source_roi_data,2);
    snip_eps = snippet_length/(pnts/EEG.roi.srate); % snip length/epoch length (how many trials for each snippet)
    nsnips = floor(trials/snip_eps); 
    if nsnips < 1
        if strcmpi(g.errornosnippet, 'on')
            error('Snippet length cannot exceed data length.\n')
        else
            fprintf(2, 'Snippet length cannot exceed data length, using the whole data\n')
            nsnips = 1;
        end
    end
    diff = (trials * pnts/EEG.roi.srate) - (nsnips * pnts/EEG.roi.srate * snip_eps);
    if diff ~= 0
        warning(strcat(int2str(diff), ' seconds are thrown away.'));
    end

    if strcmpi(g.firstsnippet, 'on')
        nsnips = 1;
    end
    
    % check if Parallel Processing Toolbox is available and licensed
    if license('test', 'Distrib_Computing_Toolbox') && ~isempty(ver('parallel'))
        if isfield(g, 'poolsize') && isnumeric(g.poolsize) && g.poolsize > 0
            % check if there's already an existing parallel pool
            currentPool = gcp('nocreate');
            if isempty(currentPool)
                parpool(g.poolsize);
            end
        end
    else
        disp('Parallel Processing Toolbox is not installed or licensed.');
    end

    tmplist1 = setdiff(g.methods, {'PAC'}); % list of fc metrics without PAC
    tmplist2 = intersect(g.methods, {'PAC'});
    % store each connectivity metric for each snippet in separate structure
    fc_matrices_snips = cell(nsnips, length(tmplist1));
    if ~isempty(tmplist2) 
        switch g.bs_outopts % number of PAC metrics (check documentation)
            case 1
                bs_matrices_snips = cell(nsnips, 4); 
                fns = cell(nsnips, 4);
            otherwise
                bs_matrices_snips = cell(nsnips, 2);
                fns = cell(nsnips, 2);
        end
    end
    source_roi_data_save = EEG.roi.source_roi_data;
    parfor isnip = 1:nsnips
%     for isnip = 1:nsnips
        EEG1 = EEG;
        begSnip = (isnip-1)* snip_eps + 1;
        endSnip = min((isnip-1)* snip_eps + snip_eps, size(source_roi_data_save,3));
        roi_snip = source_roi_data_save(:,:, begSnip:endSnip ); % cut source data into snippets
        EEG1.roi.source_roi_data = single(roi_snip);
        EEG1 = roi_connect(EEG1, 'morder', g.morder, 'naccu', g.naccu, 'methods', g.methods,'freqresolution', g.freqresolution, 'roi_selection', g.roi_selection); % compute connectivity over one snippet
        if ~isempty(intersect(g.methods, {'PAC'})) 
            EEG1 = roi_pac(EEG1, g.fcomb, g.bs_outopts, g.roi_selection);
        end
        if ~isempty(tmplist1)
            tmp_fc_matrices = cell(1, length(tmplist1));
            for fc = 1:length(tmplist1) 
                fc_name = g.methods{fc};
                fc_matrix = EEG1.roi.(fc_name);
                tmp_fc_matrices{fc} = fc_matrix;
            end
            fc_matrices_snips(isnip, :) = tmp_fc_matrices; 
        end
        if ~isempty(tmplist2)
            tmp_fns = fieldnames(EEG1.roi.PAC);
            tmp_bs_matrices = cell(1, length(tmp_fns));
            for bs = 1:length(tmp_fns)
                bs_matrix = EEG1.roi.PAC.(tmp_fns{bs});
                tmp_bs_matrices{bs} = bs_matrix;
            end
            bs_matrices_snips(isnip, :) = tmp_bs_matrices; 
            fns(isnip, :) = tmp_fns;
        end
    end

    % shut down current parallel pool only if the toolbox is available
    if license('test', 'Distrib_Computing_Toolbox') && ~isempty(ver('parallel'))
        poolobj = gcp('nocreate');
        if ~isempty(poolobj)
            delete(poolobj);
        end
    end
    
    % compute mean over connectivity of each snippet
    if ~isempty(tmplist1)
        for fc = 1:length(tmplist1)
            fc_name = g.methods{fc};
    
            [first_dim, second_dim, third_dim] = size(fc_matrices_snips{1,fc});
            conn_cell = fc_matrices_snips(:, fc); % store all matrices of one metric in a cell
            mat = cell2mat(conn_cell);
            reshaped = reshape(mat, first_dim, nsnips, second_dim, third_dim);
            reshaped = squeeze(permute(reshaped, [2, 1, 3, 4]));
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
    end
    if ~isempty(tmplist2)
        fns = fns(1, :);
        for bs = 1:length(fns)
            [second_dim, third_dim] = size(bs_matrices_snips{1, bs});
            conn_cell = bs_matrices_snips(:, bs); % store all matrices of one metric in a cell
            mat = cell2mat(conn_cell);
            reshaped = reshape(mat, second_dim, nsnips, third_dim);
            reshaped = squeeze(permute(reshaped, [2, 1, 3]));
            if strcmpi(g.fcsave_format, 'all_snips')
                EEG.roi.PAC.(fns{bs}) = reshaped;
            else
                if nsnips > 1
                    mean_conn = squeeze(mean(reshaped, 1));
                else
                    mean_conn = reshaped;
                end
                EEG.roi.PAC.(fns{bs}) = mean_conn; % store mean connectivity in EEG struct
            end
        end
    end
end

% TO-DO: add snippet option for stats mode
if strcmpi(g.conn_stats, 'on')
    EEG = roi_connstats(EEG, 'methods', g.methods, 'nshuf', g.nshuf, 'roi_selection', g.roi_selection, 'freqresolution', g.freqresolution, 'poolsize', g.poolsize, 'fcomb', g.fcomb);
end
if strcmpi(g.snippet, 'off')
    EEG = roi_connect(EEG, 'morder', g.morder, 'naccu', g.naccu, 'methods', g.methods,'freqresolution', g.freqresolution, 'roi_selection', g.roi_selection);
    if strcmpi(g.snippet, 'off') && ~isempty(intersect(g.methods, {'PAC'}))
        EEG = roi_pac(EEG, g.fcomb, g.bs_outopts, g.roi_selection);
    end
end
if nargout > 1
    com = sprintf( 'EEG = pop_roi_connect(EEG, %s);', vararg2str( options ));
end

