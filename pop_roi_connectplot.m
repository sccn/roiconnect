% pop_roi_connectplot - plot results of connectivity analysis computed
%                       by roi_connect.
% Usage:
%  pop_roi_connectplot(EEG, 'key', 'val', ...);
%
% Inputs:
%  EEG - EEGLAB dataset
%
% Required inputs:
%  'headmodel'   - [string] head model file in MNI space
%  'sourcemodel' - [string] source model file
% 
% Optional inputs:
%  'measure'    - ['psd'|'roipsd'|'trgc'|'crossspecimag'|'crossspecpow']
%                   'psd'   : Source power spectrum
%                   'psdroi': ROI based power spectrum
%                   'trgc'  : Time-reversed granger causality
%                   'crossspecimag': Imaginary part of coherence from cross-spectrum
%                   'crossspecpow' : Average cross-spectrum power for each ROI
%  'freqrange'  - [min max] frequency range in Hz. Default is to plot
%                 broadband power.
%  'smooth'     - [float] smoothing factor for cortex surface plotting
%  'plotcortex' - ['on'|'off'] plot results on smooth cortex. Default is 'on'
%  'plotmatrix' - ['on'|'off'] plot results on smooth cortex. Default is 'off'
%  'plotpsd'    - ['on'|'off'] plot PSD (for 'crossspecpow' only). Default is 'off'
%
% Author: Stefan Haufe and Arnaud Delorme, 2019
%
% Example:
%   % Requires prior call to pop_roi_connect
%   EEG = pop_roi_connectplot(EEG, 'measure', 'psd');

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

function [matrix, com] = pop_roi_connectplot(EEG, varargin)

matrix = [];
com = '';
if nargin < 1
    help pop_roi_connectplot;
    return
end

if ~isfield(EEG, 'roi')
    error('Compute connectivity first');
end

% if ~isfield(EEG.dipfit, 'hdmfile')
%     error('You need to select a head model file using DIPFIT settings first');
% end
% 
% if ~isequal(EEG.dipfit.coordformat, 'MNI')
%     error('You can only use this function with MNI coordinates - change head model');
% end

cortexFlag = isfield(EEG.roi.cortex, 'Faces');

splot = [];
splot(end+1).label  = 'Source power spectrum';
splot(end  ).acronym  = 'PSD';
splot(end  ).unit   = '?'; % not used yet
splot(end  ).cortex = cortexFlag;
splot(end  ).matrix = -1;
splot(end  ).psd    = -1;

splot(end+1).label  = 'ROI based power spectrum';
splot(end  ).acronym  = 'ROIPSD';
splot(end  ).unit   = '?'; % not used yet
splot(end  ).cortex = cortexFlag;
splot(end  ).matrix = -1;
splot(end  ).psd    = -1;

splot(end+1).label    = 'ROI to ROI cross-spectrum';
splot(end  ).labelshort = 'Cross-spectrum';
splot(end  ).acronym  = 'crossspecpow';
splot(end  ).unit   = 'Power (dB)';
splot(end  ).cortex = cortexFlag;
splot(end  ).matrix = -1;
splot(end  ).psd    = 0;

splot(end+1).label    = 'ROI to ROI imaginary part of cross-spectrum';
splot(end  ).labelshort = 'Img. part of cross-spectrum';
splot(end  ).acronym  = 'crossspecimag';
splot(end  ).unit   = 'net |iCOH|';
splot(end  ).cortex = cortexFlag;
splot(end  ).matrix = 1;
splot(end  ).psd    = -1;

splot(end+1).label    = 'ROI to ROI coherence';
splot(end  ).labelshort = 'Coherence';
splot(end  ).acronym  = 'Coh';
splot(end  ).unit   = '?';
splot(end  ).cortex = cortexFlag;
splot(end  ).matrix = 1;
splot(end  ).psd    = -1;

splot(end+1).label  = 'ROI to ROI granger causality';
splot(end  ).labelshort = 'Granger Causality';
splot(end  ).acronym  = 'GC';
splot(end  ).cortex = cortexFlag;
splot(end  ).matrix = 1;
splot(end  ).psd    = -1;

splot(end+1).label  = 'ROI to ROI time-reversed granger causality';
splot(end  ).labelshort = 'Time-rev. Granger Causality';
splot(end  ).acronym  = 'TRGC';
splot(end  ).unit   = '?'; % not used yet
splot(end  ).cortex = cortexFlag;
splot(end  ).matrix = 1;
splot(end  ).psd    = -1;

splot(end+1).label    = 'ROI to ROI mutual information C';
splot(end  ).labelshort = 'Mutual information C';
splot(end  ).acronym  = 'MIC';
splot(end  ).unit   = '?'; % not used yet
splot(end  ).cortex = cortexFlag;
splot(end  ).matrix = -1;
splot(end  ).psd    = 0;

splot(end+1).label    = 'ROI to ROI mutual information M';
splot(end  ).labelshort = 'Mutual information I';
splot(end  ).acronym  = 'MIM';
splot(end  ).unit   = '?'; % not used yet
splot(end  ).cortex = cortexFlag;
splot(end  ).matrix = -1;
splot(end  ).psd    = 0;

if nargin < 2

    cb_select = [ 'usrdat = get(gcf, ''userdata'');' ...
                  'usrdat = usrdat(get(findobj(gcf, ''tag'', ''selection''), ''value''));' ...
                  'fieldTmp = { ''cortex'' ''matrix'' ''psd'' };' ...
                  'for iField = 1:length(fieldTmp),' ...
                  '   if usrdat.(fieldTmp{iField}) == 1,' ...
                  '       set(findobj(gcf, ''tag'', fieldTmp{iField}), ''enable'', ''on'', ''value'', 1);' ...
                  '   elseif usrdat.(fieldTmp{iField}) == 0,' ...
                  '       set(findobj(gcf, ''tag'', fieldTmp{iField}), ''enable'', ''on'', ''value'', 0);' ...
                  '   else,' ...
                  '       set(findobj(gcf, ''tag'', fieldTmp{iField}), ''enable'', ''off'', ''value'', 0);' ...
                  '   end;' ...
                  'end;' ...
                  'clear iField fieldTmp usrdat;' ];

    uigeom = { [1 1] [1 1] 1 [0.3 1.1 1 1] };
    uilist = {{ 'style' 'text' 'string' 'Select a measure to plot' 'fontweight' 'bold'} ...
              { 'style' 'popupmenu' 'string' {splot.label} 'callback' cb_select 'value' 4 'tag' 'selection' } ...
              { 'style' 'text' 'string' 'Frequency range in Hz [min max]:'} ...
              { 'style' 'edit' 'string' ''} ...
              {} ...
              {} ...
              { 'style' 'checkbox' 'string' 'Plot on cortex' 'tag' 'cortex' 'value' 1 } ...
              { 'style' 'checkbox' 'string' 'Plot PSD' 'tag' 'psd'  'enable' 'off'  } ...
              { 'style' 'checkbox' 'string' 'Plot in matrix' 'tag' 'matrix' 'enable' 'off' } ...
              };
              
    [result,~,~,outs] = inputgui('geometry', uigeom, 'uilist', uilist, 'helpcom', 'pophelp(''pop_loadbv'')', ...
        'title', 'ROI connectivity', 'userdata', splot, 'eval', cb_select);
    if isempty(result), return, end
    
    options = {};
    options = { options{:} 'measure'   splot(result{1}).acronym };
    options = { options{:} 'freqrange' eval( [ '[' result{2} ']' ] ) };
    options = { options{:} 'plotcortex' fastif(outs.cortex, 'on', 'off') };
    options = { options{:} 'plotmatrix' fastif(outs.matrix, 'on', 'off') };
    options = { options{:} 'plotpsd'    fastif(outs.psd   , 'on', 'off') };
else
    options = varargin;
end

% decode input parameters
% -----------------------
g = finputcheck(options,  { 'measure'    'string'  {splot.acronym}  '';
                            'freqrange'  'real'    { }              [];
                            'smooth'     'real'    { }              0.35;
                            'plotcortex' 'string'  { 'on' 'off' }   'on';
                            'plotmatrix' 'string'  { 'on' 'off' }   'off';
                            'plotpsd'    'string'  { 'on' 'off' }   'off' }, 'pop_roi_connectplot');
if ischar(g), error(g); end
S = EEG.roi;

% colormap
load cm17;

% frequency range
if ~isempty(g.freqrange)
    frq_inds = find(S.freqs >= g.freqrange(1) & S.freqs < g.freqrange(2));
    titleStr = sprintf('%1.1f-%1.1f Hz frequency band', g.freqrange(1), g.freqrange(2));
else
    frq_inds = 1:length(S.freqs);
    titleStr = 'broadband';
end

% plotting options
allMeasures = { splot.acronym };
pos = strmatch( g.measure, allMeasures, 'exact');
plotOpt = splot(pos);

switch lower(g.measure)
    case { 'psd' 'roipsd' }
        if ~isfield(S, 'source_voxel_data')
            error('Cannot plot spectrum without source data');
        end
        
        if ~isempty(g.freqrange)
            % filter range of interest
            freqRatio = g.freqrange/max(S.freqs);
            [B,A] = butter(5,freqRatio); % sampling rate is 1 so divide by 2 as in S.freqs
            disp('Filtering...');
            X = filtfilt(B, A, double(S.source_voxel_data));
            
            disp('Applying hilbert transform...');
            Xhilbert = hilbert(X);
            P = reshape(sum(sum(Xhilbert.*conj(Xhilbert), 1), 3), [], 1);
        else
            P = reshape(sum(sum(S.source_voxel_data.^2, 1), 3), [], 1);
        end
        P_dB = 10*log10( P );
        
        if strcmpi(g.plotcortex, 'on')
            if strcmpi(lower(g.measure), 'roipsd')
                source_roi_power = zeros(1,S.nROI);
                for iROI = 1:S.nROI
                    source_roi_power(iROI) = mean(P(S.atlas.Scouts(iROI).Vertices));
                end
                source_roi_power_norm_dB = 10*log10( source_roi_power );
                allplots_cortex_BS(S.cortex, source_roi_power_norm_dB, [min(source_roi_power_norm_dB) max(source_roi_power_norm_dB)], cm17a, 'power [dB]', g.smooth);
                h = textsc([ 'ROI source power (' titleStr ')' ], 'title');
                set(h, 'fontsize', 20);
            else
                allplots_cortex_BS(S.cortex, P_dB, [min(P_dB) max(P_dB)], cm17a, 'power [dB]', g.smooth);
                h = textsc([ 'Source power (' titleStr ')' ], 'title');
                set(h, 'fontsize', 20);
            end
        end
        
    case { 'trgc' 'gc' }
        % calculation of net TRGC scores (i->j minus j->i), recommended procedure
        % TRGCnet = TRGC_(:, 1:2:end)-TRGC_(:, 2:2:end);
        % new way to compute net scores
        if strcmpi(g.measure, 'GC')
            TRGCnet = S.GC(:, :, 1) - S.GC(:, :, 2);
        else
            TRGCnet = S.TRGC(:, :, 1) - S.TRGC(:, :, 2);
        end
        TRGC = get_connect_mat( TRGCnet, S.nROI, -1);
        
        if strcmpi(g.plotmatrix, 'on')
            matrix = squeeze(mean(TRGC(frq_inds, :, :)));
            figure; imagesc(matrix); colorbar
            xlabel('ROI index (see Atlas for more info)');
            h = title([ 'ROI to ROI ' upper(g.measure) ' (' titleStr ')' ]);
            set(h, 'fontsize', 16);
        end        
        
        if strcmpi(g.plotcortex, 'on') 
            atrgc = mean(squeeze(mean(TRGC(frq_inds, :, :))), 2);
            allplots_cortex_BS(S.cortex, atrgc, [-max(abs(atrgc)) max(abs(atrgc))], cm17, upper(g.measure), g.smooth);
            h = textsc([ upper(g.measure) ' (' titleStr '); Red = net sender; Blue = net receiver' ], 'title');
            set(h, 'fontsize', 20);
        end
        
    case { 'mim' 'mic' }
        if strcmpi(g.measure, 'MIC')
            MI = S.MIC(:, :);
        else
            MI = S.MIM(:, :);
        end
        MI = get_connect_mat( MI, S.nROI, +1);

        if strcmpi(g.plotmatrix, 'on')
            matrix = squeeze(mean(MI(frq_inds, :, :)));
            figure; imagesc(matrix); colorbar
            xlabel('ROI index (see Atlas for more info)');
            h = title(['ROI to ROI imag. part of coherence (' titleStr ')']);
            set(h, 'fontsize', 16);
        end        
        
        if strcmpi(g.plotcortex, 'on')
            ami = mean(squeeze(mean(MI(frq_inds, :, :))), 2);
            allplots_cortex_BS(S.cortex, ami, [min(ami) max(ami)], cm17, upper(g.measure), g.smooth);
            h = textsc([ upper(g.measure) ' (' titleStr '); Red = net sender; Blue = net receiver' ], 'title');
            set(h, 'fontsize', 20);
        end
        
    case { 'crossspecpow' 'coh' 'crossspecimag' }
        if strcmpi(g.measure, 'coh')
            PS = S.COH; % do not know what to do here
            PS = squeeze(mean(mean(reshape(PS, S.srate+1, 3, S.nROI, 3, S.nROI), 2), 4));     
            PSmean = mean(squeeze(mean(PS(frq_inds, :, :))), 2);
        elseif strcmpi(g.measure, 'crossspecimag')
            PS = abs(imag(cs2coh(S.CS)));
            PS = squeeze(mean(mean(reshape(PS, S.srate+1, 3, S.nROI, 3, S.nROI), 2), 4));     
            PSmean = mean(squeeze(mean(PS(frq_inds, :, :))), 2);
        else
            PS = cs2psd(S.CS);
            apow = squeeze(sum(sum(reshape(PS(frq_inds, :), [], S.nPCA, S.nROI), 1), 2)).*S.source_roi_power_norm';
            PSmean = 10*log10(apow);
        end
        
        if strcmpi(g.plotmatrix, 'on')
            matrix = squeeze(mean(PS(frq_inds, :, :)));
            figure; imagesc(matrix); colorbar
            xlabel('ROI index (see Atlas for more info)');
            h = title([ plotOpt.label ' (' titleStr ')']);
            set(h, 'fontsize', 16);
        end
        
        if strcmpi(g.plotcortex, 'on')
            allplots_cortex_BS(S.cortex, PSmean, [min(PSmean) max(PSmean)], cm17a, plotOpt.unit, g.smooth);
            h = textsc([ plotOpt.labelshort ' (' titleStr ')' ], 'title');
            set(h, 'fontsize', 20);
        end
        
        if strcmpi(g.plotpsd, 'on') 
            figure; semilogy(S.freqs(frq_inds), PSmean(frq_inds, :)); grid on
            h = textsc(plotOptS.label, 'title');
            set(h, 'fontsize', 20);
        end
        
end

if nargin < 2
    com = sprintf('pop_roi_connectplot(EEG, %s);', vararg2str( options ));
end

function measure = get_connect_mat( measureOri, nROI, signVal)
        % create a ROI x ROI connectivity matrix, if needed
        % TRGCmat(f, ii, jj) is net TRGC from jj to ii
        measure = [];
        iinds = 0;
        for iroi = 1:nROI
            for jroi = (iroi+1):nROI
                iinds = iinds + 1;
                measure(:, iroi, jroi) = signVal * measureOri(:, iinds);
                measure(:, jroi, iroi) = measureOri(:, iinds);
            end
        end

