% pop_roi_connectivity_plot - plot results of connectivity analysis computed
%                             by roi_connectivity_process.
% Usage:
%  pop_roi_connectivity_plot(EEG, 'key', 'val', ...);
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
%   % Requires prior call to pop_roi_connectivity_process
%   EEG = pop_roi_connectivity_plot(EEG, 'measure', 'psd');

function com = pop_roi_connectivity_plot(EEG, varargin)

com = '';
if nargin < 1
    help pop_roi_connectivity_plot;
    return
end

if ~isfield(EEG, 'roiconnect')
    error('Compute connectivity first');
end

% if ~isfield(EEG.dipfit, 'hdmfile')
%     error('You need to select a head model file using DIPFIT settings first');
% end
% 
% if ~isequal(EEG.dipfit.coordformat, 'MNI')
%     error('You can only use this function with MNI coordinates - change head model');
% end

splot(1).label  = 'Source power spectrum';
splot(1).acronym  = 'PSD';
splot(1).cortex = 1;
splot(1).matrix = -1;
splot(1).psd    = -1;

splot(2).label  = 'ROI based power spectrum';
splot(2).acronym  = 'ROIPSD';
splot(2).cortex = 1;
splot(2).matrix = -1;
splot(2).psd    = -1;

splot(3).label  = 'ROI to ROI time-reversed granger causality';
splot(3).acronym  = 'TRGC';
splot(3).cortex = 1;
splot(3).matrix = 1;
splot(3).psd    = -1;

splot(4).label    = 'ROI to ROI imaginary part of cross-spectrum';
splot(4).acronym  = 'crossspecimag';
splot(4).cortex = 1;
splot(4).matrix = 1;
splot(4).psd    = -1;

splot(5).label    = 'ROI cross-spectrum';
splot(5).acronym  = 'crossspecpow';
splot(5).cortex = 1;
splot(5).matrix = -1;
splot(5).psd    = 0;

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
                            'plotpsd'    'string'  { 'on' 'off' }   'off' }, 'pop_roi_connectivity');
if ischar(g), error(g); end
S = EEG.roiconnect;

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

switch lower(g.measure)
    case { 'psd' 'roipsd' }
        if ~isempty(g.freqrange)
            % filter range of interest
            freqRatio = g.freqrange/max(S.freqs);
            [B,A] = butter(5,freqRatio); % sampling rate is 1 so divide by 2 as in S.freqs
            disp('Filtering...');
            X = filtfilt(B, A, S.source_voxel_data);
            
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
        
    case 'trgc'
        if strcmpi(g.plotmatrix, 'on') 
            figure; imagesc(squeeze(mean(S.TRGCmat(frq_inds, :, :)))); colorbar
            xlabel('ROI index (see Atlas for more info)');
            h = title([ 'ROI to ROI TRGC (' titleStr ')' ]);
            set(h, 'fontsize', 16);
        end        
        
        if strcmpi(g.plotcortex, 'on') 
            atrgc = squeeze(mean(mean(S.TRGCmat(frq_inds, :, :), 1), 2));
            allplots_cortex_BS(S.cortex, atrgc, [-max(abs(atrgc)) max(abs(atrgc))], cm17, 'TRGC', g.smooth);
            h = textsc([ 'TRGC (' titleStr '); Red = net sender; Blue = net receiver' ], 'title');
            set(h, 'fontsize', 20);
        end

    case 'crossspecimag'
        absiCOH = abs(imag(cs2coh(S.CS)));
        absiCOH = squeeze(mean(mean(reshape(absiCOH, S.srate+1, 3, S.nROI, 3, S.nROI), 2), 4));

        if strcmpi(g.plotmatrix, 'on')
            figure; imagesc(squeeze(mean(absiCOH(frq_inds, :, :)))); colorbar
            xlabel('ROI index (see Atlas for more info)');
            h = title(['ROI to ROI imag. part of coherence (' titleStr ')']);
            set(h, 'fontsize', 16);
        end        
        
        if strcmpi(g.plotcortex, 'on')
            netabsiCOH = mean(squeeze(mean(absiCOH(frq_inds, :, :))), 2);
            
            % red = highly interacting ROI
            allplots_cortex_BS(S.cortex, netabsiCOH, [min(netabsiCOH) max(netabsiCOH)], cm17a, 'net |iCOH|', g.smooth);
            h = textsc([ 'Imag. part of coherence (' titleStr '); red = highly interacting ROI'], 'title');
            set(h, 'fontsize', 20);
        end
        
    case 'crossspecpow'
        PS = cs2psd(S.CS);
        apow = squeeze(sum(sum(reshape(PS(frq_inds, :), [], S.nPCA, S.nROI), 1), 2)).*S.source_roi_power_norm';
        apow_dB = 10*log10(apow);
        
        if strcmpi(g.plotcortex, 'on') 
            allplots_cortex_BS(S.cortex, apow_dB, [min(apow_dB) max(apow_dB)], cm17a, 'power [dB]', g.smooth);
            h = textsc([ 'Cross-spectrum (' titleStr ')' ], 'title');
            set(h, 'fontsize', 20);
        end
        
        if strcmpi(g.plotpsd, 'on') 
            figure; semilogy(S.freqs(frq_inds), PS(frq_inds, :)); grid on
            h = textsc('Broadband ROI cross-spectrum average', 'title');
            set(h, 'fontsize', 20);
        end
        
end

if nargin < 2
    com = sprintf('pop_roi_connectivity_plot(EEG, %s);', vararg2str( options ));
end
