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
%  'measure'              - ['psd'|'roipsd'|'trgc'|'crossspecimag'|'crossspecpow'|'mic'|'mim']
%                           'psd'   : Source power spectrum
%                           'psdroi': ROI based power spectrum
%                           'TRGC'  : Time-reversed granger causality
%                           'GC'    : Granger causality
%                           'crossspecimag': Imaginary part of coherence from cross-spectrum
%                           'crossspecpow' : Average cross-spectrum power for each ROI
%                           'aCOH': Coherence 
%                           'iCOH': Absolute value of the imaginary part of Coherency
%                           'MIC' : Maximized Imaginary Coherency for each ROI
%                           'MIM' : Multivariate Interaction Measure for each ROI
%                           'PAC' : Phase-amplitude coupling for a certain frequency (band) combination based on bicoherence
%                           'PAC_anti': Phase-amplitude coupling for a certain frequency (band) combination based on the antisymmetrized bicoherence
%                           'TDE': Bispectral time-delay estimation
%                           'TDE_anti': Antisymmetrized bispectral time-delay estimation
%  'freqrange'            - [min max] frequency range or [integer] single frequency in Hz. Default is to plot broadband power.
%  'smooth'               - [float] smoothing factor for cortex surface plotting
%  'plotcortex'           - ['on'|'off'] plot results on smooth cortex. Default is 'on'
%  'plotcortexparams'     - [cell] ...
%  'plotcortexseedregion' - [string] plot seed voxel on cortex. Takes name of seed region as input.
%  'plot3d'               - ['on'|'off'] ... Default is 'off'
%  'plot3dparams'         - [cell] optional parameters for the generation of brain movies. Check the related documentation in roi_plotbrainmovie.m
%  'plotmatrix'           - ['on'|'off'] plot results as ROI to ROI matrix. Default is 'off'
%  'plotbutterfly'        - ['on'|'off'] plot results as a frequency x connectivity plot (butterfly plot). Default is 'off'
%  'plotbarplot'          - ['on'|'off'] plot ROI based power spectrum as barplot. Default is 'off'
%  'hemisphere'           - ['all'|'left'|'right'] hemisphere options for ROI to ROI matrix. Default is 'all'
%  'grouphemispheres'      - ['on'|'off'] group ROIs by hemispheres (left hemisphere, then right hemisphere). Default is 'off'
%  'region'               - ['all'|'cingulate'|'prefrontal'|'frontal'|'temporal'|'parietal'|'central'|'occipital'] region selection for ROI to ROI matrix. Default is 'all'
%  'largeplot'            - ['on'|'off'] plot MIM, TRGC and Power in a single large plot. Default is 'off'
%  'plotpsd'              - ['on'|'off'] plot PSD (for 'crossspecpow' only). Default is 'off'
%  'noplot'               - ['on'|'off'] when 'on', disable all plotting. Default is 'off'
%
% Output
%  matnet - connectivity matrix (nROI x nROI)
%
% Author: Stefan Haufe and Arnaud Delorme, 2019
%
% Example:
%   % Requires prior call to pop_roi_connect
%   matnet = pop_roi_connectplot(EEG, 'measure', 'ROIPSD');
%
% See also: STD_ROI_CONNECTPLOT, ROI_ACTIVITY, POP_ROI_ACTIVITY, POP_LEADFIELD

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

    if ~exist('roi_plotbrainmovie')
        fprintf(2, 'To plot connectivity in 3-D, install the brainmovie plugin\n');
        plot3dFlag = -1;
    else
        plot3dFlag = 0;
    end


    % if ~isfield(EEG.dipfit, 'hdmfile')
    %     error('You need to select a head model file using DIPFIT settings first');
    % end
    %
    % if ~isequal(EEG.dipfit.coordformat, 'MNI')
    %     error('You can only use this function with MNI coordinates - change head model');
    % end

    cortexFlag = isfield(EEG.roi.cortex, 'Faces')*2-1;
    areaNames  = [ 'All' { EEG.roi.atlas(1).Scouts.Label } ];

    splot = [];
    % splot(end+1).label  = 'Source power spectrum'; % we do not save that information anymore
    % splot(end  ).acronym  = 'PSD';
    % splot(end  ).unit   = '?'; % not used yet
    % splot(end  ).cortex = cortexFlag;
    % splot(end  ).matrix = -1;
    % splot(end  ).psd    = -1;

    if isfield(EEG.roi, 'source_roi_power')
        splot(end+1).label  = 'ROI based power spectrum';
        splot(end  ).labelshort = 'ROI PSD';
        splot(end  ).acronym  = 'ROIPSD';
        splot(end  ).unit   = 'Power (dB)'; % not used yet
        splot(end  ).cortex = cortexFlag;
        splot(end  ).matrix = -1;
        splot(end  ).psd    = -1;
        splot(end  ).plot3d = -1;
    end

    if isfield(EEG.roi, 'CS')
        splot(end+1).label    = 'ROI to ROI cross-spectrum';
        splot(end  ).labelshort = 'Cross-spectrum';
        splot(end  ).acronym  = 'crossspecpow';
        splot(end  ).unit   = 'Power (dB)';
        splot(end  ).cortex = cortexFlag;
        splot(end  ).matrix = -1;
        splot(end  ).psd    = 0;
        splot(end  ).plot3d = plot3dFlag;
    end

    if isfield(EEG.roi, 'CS')
        splot(end+1).label    = 'ROI to ROI imaginary part of cross-spectrum';
        splot(end  ).labelshort = 'Img. part of cross-spectrum';
        splot(end  ).acronym  = 'crossspecimag';
        splot(end  ).unit   = 'net |iCOH|';
        splot(end  ).cortex = cortexFlag;
        splot(end  ).matrix = 1;
        splot(end  ).psd    = -1;
        splot(end  ).plot3d = plot3dFlag;
    end

    if isfield(EEG.roi, 'aCOH')
        splot(end+1).label    = 'ROI to ROI coherence';
        splot(end  ).labelshort = 'Coherence';
        splot(end  ).acronym  = 'aCOH';
        splot(end  ).unit   = 'aCOH';
        splot(end  ).cortex = cortexFlag;
        splot(end  ).matrix = 1;
        splot(end  ).psd    = 0;
        splot(end  ).plot3d = plot3dFlag;
    end

    if isfield(EEG.roi, 'cCOH')
        splot(end+1).label    = 'ROI to ROI coherency';
        splot(end  ).labelshort = 'Coherency';
        splot(end  ).acronym  = 'cCOH';
        splot(end  ).unit   = 'cCOH';
        splot(end  ).cortex = cortexFlag;
        splot(end  ).matrix = 1;
        splot(end  ).psd    = 0;
        splot(end  ).plot3d = plot3dFlag;
        end

    if isfield(EEG.roi, 'iCOH')
        splot(end+1).label    = 'ROI to ROI absolute value of the imaginary part of coherency';
        splot(end  ).labelshort = 'Img. part of Coherency';
        splot(end  ).acronym  = 'iCOH';
        splot(end  ).unit   = 'iCOH';
        splot(end  ).cortex = cortexFlag;
        splot(end  ).matrix = 1;
        splot(end  ).psd    = 0;
        splot(end  ).plot3d = plot3dFlag;
    end

    if isfield(EEG.roi, 'GC')
        splot(end+1).label  = 'ROI to ROI granger causality';
        splot(end  ).labelshort = 'Granger Causality';
        splot(end  ).acronym  = 'GC';
        splot(end  ).unit   = 'GC';
        splot(end  ).cortex = cortexFlag;
        splot(end  ).matrix = 1;
        splot(end  ).psd    = 0;
        splot(end  ).plot3d = plot3dFlag;
    end

    if isfield(EEG.roi, 'TRGC')
        splot(end+1).label  = 'ROI to ROI time-reversed granger causality';
        splot(end  ).labelshort = 'Time-rev. Granger Causality';
        splot(end  ).acronym  = 'TRGC';
        splot(end  ).unit   = 'TRGC'; % not used yet
        splot(end  ).cortex = cortexFlag;
        splot(end  ).matrix = 1;
        splot(end  ).psd    = 0;
        splot(end  ).plot3d = plot3dFlag;
    end

    if isfield(EEG.roi, 'MIC')
        splot(end+1).label    = 'ROI to ROI Maximized Imag. Coh.';
        splot(end  ).labelshort = 'Maximized Imag. Coh.';
        splot(end  ).acronym  = 'MIC';
        splot(end  ).unit   = 'MIC'; % not used yet
        splot(end  ).cortex = cortexFlag;
        splot(end  ).matrix = -1;
        splot(end  ).psd    = 0;
        splot(end  ).plot3d = plot3dFlag;
    end

    if isfield(EEG.roi, 'MIM')
        splot(end+1).label    = 'ROI to ROI Multivariate Interaction Measure';
        splot(end  ).labelshort = 'Multivariate Interaction Measure';
        splot(end  ).acronym  = 'MIM';
        splot(end  ).unit   = 'MIM'; % not used yet
        splot(end  ).cortex = cortexFlag;
        splot(end  ).matrix = 1;
        splot(end  ).psd    = 0;
        splot(end  ).plot3d = plot3dFlag;
    end

    if isfield(EEG.roi, 'PAC')
        splot(end+1).label    = 'ROI to ROI Phase-amplitude coupling';
        splot(end  ).labelshort = 'Phase-amplitude coupling';
        splot(end  ).acronym  = 'PAC'; % PAC based on bicoherence
        splot(end  ).unit   = 'PAC'; % not used yet
        splot(end  ).cortex = cortexFlag;
        splot(end  ).matrix = 1;
        splot(end  ).psd    = 0;
        splot(end  ).plot3d = plot3dFlag;
    end

    if isfield(EEG.roi, 'PAC')
        splot(end+1).label    = 'ROI to ROI Phase-amplitude coupling';
        splot(end  ).labelshort = 'Phase-amplitude coupling';
        splot(end  ).acronym  = 'PAC_anti'; % PAC based on antisymmetrized bicoherence
        splot(end  ).unit   = 'PAC'; % not used yet
        splot(end  ).cortex = cortexFlag;
        splot(end  ).matrix = 1;
        splot(end  ).psd    = 0;
        splot(end  ).plot3d = plot3dFlag;
    end

    if isfield(EEG.roi, 'TDE')
        splot(end+1).label    = 'ROI to ROI Time-delay estimation';
        splot(end  ).labelshort = 'Time-delay estimation';
        splot(end  ).acronym  = 'TDE'; % TDE based on bispectrum
        splot(end  ).unit   = 'TDE'; % not used yet
        splot(end  ).cortex = cortexFlag;
        splot(end  ).matrix = 1;
        splot(end  ).psd    = -1;
        splot(end  ).plot3d = plot3dFlag;
    end

    if isfield(EEG.roi, 'TDE')
        splot(end+1).label    = 'ROI to ROI Time-delay estimation';
        splot(end  ).labelshort = 'Time-delay estimation';
        splot(end  ).acronym  = 'TDE_anti'; % TDE based on antisymmetrized bispectrum
        splot(end  ).unit   = 'TDE'; % not used yet
        splot(end  ).cortex = cortexFlag;
        splot(end  ).matrix = 1;
        splot(end  ).psd    = -1;
        splot(end  ).plot3d = plot3dFlag;
    end
   
    if nargin < 2

        cb_select = [ 'usrdat = get(gcf, ''userdata'');' ...
            'usrdat = usrdat(get(findobj(gcf, ''tag'', ''selection''), ''value''));' ...
            'fieldTmp = { ''cortex'' ''matrix'' ''psd'' ''plot3d'' };' ...
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
        
        fcregions = {'all', 'cingulate', 'prefrontal', 'frontal', 'temporal', 'parietal', 'central', 'occipital'};
        plotrow = [1 1];
        uigeom = { [1 1] [1 1] [1 0.5 0.2] 1 [1 1] plotrow [1] [3.3 1.5 1.1 1.2] plotrow };
        uilist = {{ 'style' 'text' 'string' 'Select a measure to plot' 'fontweight' 'bold'} ...
            { 'style' 'popupmenu' 'string' {splot.label} 'callback' cb_select 'value' 1 'tag' 'selection' } ...
            { 'style' 'text' 'string' 'Only from/to region' } ...
            { 'style' 'popupmenu' 'string' areaNames 'tag' 'seed_region'} ...
            { 'style' 'text' 'string' 'Frequency range in Hz [min max]:'} ...
            { 'style' 'edit' 'string' '' 'tag' 'freqs'} {} ...
            {} ...
            { 'style' 'text' 'string' 'Measure to plot' 'fontweight' 'bold' } ...
            { 'style' 'text' 'string' 'Measure parameters' 'fontweight' 'bold' } ...
            ...
            { 'style' 'checkbox' 'string' '3d static brainmovie' 'tag' 'plot3d' 'value' 1 } ...
            { 'style' 'edit'     'string' '''thresholdper'', 0.2' 'tag' 'plot3dparams' } ...
            ...
            { 'style' 'checkbox' 'string' 'Connectivity on cortex (requires surface atlas)' 'tag' 'cortex' 'value' 1 } ...
            ...
            { 'style' 'checkbox' 'string' 'Matrix representation' 'tag' 'matrix' 'enable' 'off' } ...
            { 'style' 'popupmenu' 'string' fcregions 'callback' cb_select 'value' 1 'tag' 'region' } ....
            { 'style' 'checkbox' 'string' 'left' 'tag' 'hemisphere_left' 'value' 1 } ...
            { 'style' 'checkbox' 'string' 'right' 'tag' 'hemisphere_right' 'value' 1 } ...
            ...
            { 'style' 'checkbox' 'string' 'Power spectral density' 'tag' 'psd'  'enable' 'off'  } {} ...
            };

        [result,~,~,outs] = inputgui('geometry', uigeom, 'uilist', uilist, 'helpcom', 'pophelp(''pop_loadbv'')', ...
            'title', 'ROI connectivity', 'userdata', splot, 'eval', cb_select);
        if isempty(result), return, end

        options = {};
        options = { options{:} 'measure'   splot(result{1}).acronym };
        options = { options{:} 'freqrange' eval( [ '[' outs.freqs ']' ] ) };
        options = { options{:} 'plotcortex' fastif(outs.cortex, 'on', 'off') };
        options = { options{:} 'plotcortexparams' {} };
        options = { options{:} 'plotcortexseedregion' outs.seed_region-1 };
        options = { options{:} 'plotmatrix' fastif(outs.matrix, 'on', 'off') };
        options = { options{:} 'plotpsd'    fastif(outs.psd   , 'on', 'off') };
        options = { options{:} 'plot3d'     fastif(outs.plot3d, 'on', 'off') };
        options = { options{:} 'plot3dparams' eval( [ '{' outs.plot3dparams '}' ] ) };
        options = { options{:} 'region' fcregions{outs.region} };
        % choose which hemisphere to plot
        if outs.hemisphere_left == 1 && outs.hemisphere_right == 0
            options = { options{:} 'hemisphere' 'left' };
        elseif outs.hemisphere_left == 0 && outs.hemisphere_right == 1
            options = { options{:} 'hemisphere' 'right' };
        else
            options = { options{:} 'hemisphere' 'all' };
        end
    else
        options = varargin;
    end

    % decode input parameters
    % -----------------------
    g = finputcheck(options,  { 'measure'    'string'  {splot.acronym}  '';
        'freqrange'             'real'     { }                     [];
        'smooth'                'real'     { }                     0.35;
        'plotcortex'            'string'   { 'on' 'off' }          'on';
        'plotcortexparams'      'cell'     { }                     {};
        'plotcortexseedregion'  'integer'  { }                     [];
        'plot3d'                'string'   { 'on' 'off' }          'off';
        'plot3dparams'          'cell'     { }                     {};
        'plotmatrix'            'string'   { 'on' 'off' }          'off';
        'plotbutterfly'         'string'   { 'on' 'off' }          'off';
        'noplot'                'string'   { 'on' 'off' }          'off';
        'plotbarplot'           'string'   { 'on' 'off'}           'off';
        'hemisphere'            'string'   {'all' 'left' 'right'}  'all';
        'grouphemispheres'      'string'   { 'on' 'off'}           'off';
        'region'                'string'   { 'all', 'cingulate', 'prefrontal', 'frontal', 'temporal', 'parietal', 'central', 'occipital' }  'all';
        'largeplot',            'string'   { 'on'  'off'  }        'off';
        'plotpsd',              'string'   { 'on' 'off' }          'off' }, 'pop_roi_connectplot');
    if ischar(g), error(g); end
    if isequal(g.plotcortexseedregion, 0)
        g.plotcortexseedregion = [];
    end
    S = EEG.roi;

    if isempty(g.measure)
        error('You must define a measure to plot');
    end
    
    if isfield(S, 'roi_selection')
        if ~isempty(S.roi_selection)
            warning("Plotting options ('region', 'hemisphere', 'grouphemispheres') are disabled when ROIs were explicitely selected.")
            g.region = 'all';
            g.hemisphere = 'all';
            g.grouphemispheres = 'all';
        end
    end
    
    % colormap
    load cm17;
    load cm18;
    
    % replace low-resolution with high-resolution cortex
    load cortex;
    
    % extract frequency indices
    if ~isempty(g.freqrange)
        if length(g.freqrange) == 1
            frq_inds = find(S.freqs == g.freqrange(1)); 
            titleStr = sprintf('%1.1f Hz', g.freqrange(1));
        else
            frq_inds = find(S.freqs >= g.freqrange(1) & S.freqs <= g.freqrange(2));
            titleStr = sprintf('%1.1f-%1.1f Hz frequency band', g.freqrange(1), g.freqrange(2));
        end
    else
        frq_inds = 1:length(S.freqs);
        titleStr = 'broadband';
    end

    % plotting options
    allMeasures = { splot.acronym };
    pos = strmatch( g.measure, allMeasures, 'exact');
    plotOpt = splot(pos);
    
    % either plot large plot with MIM, TRGC and power or only individual plots
    if strcmpi(g.largeplot, 'on')
        source_roi_power_norm_dB = 10*log10( mean(EEG.roi.source_roi_power(frq_inds,:)) ); % roipsd
        
%         TRGCnet = S.TRGC;
%         TRGCnet = TRGCnet - permute(TRGCnet, [1 3 2]);
%         TRGCnet = TRGCnet(:,:);
        
        TRGC_matrix = squeeze(mean(S.TRGC(frq_inds, :, :)));
        MIM_matrix = squeeze(mean(S.MIM(frq_inds, :, :)));
        
        roi_largeplot(EEG, MIM_matrix, TRGC_matrix, source_roi_power_norm_dB, titleStr)
    else     
        matrix = [];
        switch lower(g.measure)
            case { 'psd' 'roipsd' }
                if strcmpi(g.measure, 'psd')
                    % plot poower of individual voxels
                    % we would need to save the power in roi_activity. The function below can plot power
                    % allplots_cortex_BS(S.cortex, P_dB, [min(P_dB) max(P_dB)], cm17a, 'power [dB]', g.smooth);
                    error('This option is obsolete');
                end

                if strcmpi(g.plotcortex, 'on') && strcmpi(g.measure, 'roipsd')
                    cortexPlot = 10*log10( mean(EEG.roi.source_roi_power(frq_inds,:)) );
                end

                if strcmpi(g.plotbarplot, 'on') && ~strcmpi(g.noplot, 'on')
                    source_roi_power_norm_dB = 10*log10( mean(EEG.roi.source_roi_power(frq_inds,:)) );
                    roi_plotpower(EEG, source_roi_power_norm_dB, titleStr);
                end

            case { 'trgc' 'gc' }
                % calculation of net TRGC scores (i->j minus j->i), recommended procedure
                % new way to compute net scores
                if strcmpi(g.measure, 'GC')
%                     TRGCnet = S.GC(:, :, 1) - S.GC(:, :, 2);
                    TRGC = S.GC;
                else
%                     TRGCnet = S.TRGC(:, :, 1) - S.TRGC(:, :, 2);
                    TRGC = S.TRGC;
                end
                butterflyplot = squeeze(mean(TRGC, 2));
                matrix = squeeze(mean(TRGC(frq_inds, :, :), 1));
                cortexPlot  = mean(matrix, 2);
                cortexTitle = [ g.measure ' (' titleStr '); Red = net sender; Blue = net receiver' ];

            case { 'mim' 'mic' }
                if strcmpi(g.measure, 'MIC')
%                     MI = S.MIC(:, :);
                    MI = S.MIC;
                else
%                     MI = S.MIM(:, :);
                    MI = S.MIM;
                end
                butterflyplot = squeeze(mean(MI, 2));
                matrix = squeeze(mean(MI(frq_inds, :, :), 1));
                cortexPlot = mean(matrix, 2);
            
            case { 'acoh' 'ccoh' 'icoh'}
                if strcmpi(g.measure, 'aCOH')
                    butterflyplot = squeeze(mean(S.aCOH, 2));
                    matrix = squeeze(mean(S.aCOH(frq_inds, :, :), 1));
                elseif strcmpi(g.measure, 'cCOH')
                    error(['Complex values are not supported. To plot the absolute values, compute "aCOH", ' ...
                        'to plot the imaginary part, compute "iCOH".'])
                else
                    butterflyplot = squeeze(mean(S.iCOH, 2));
                    matrix = squeeze(mean(S.iCOH(frq_inds, :, :), 1));
                end
                cortexPlot = mean(matrix, 2);

            case { 'crossspecpow' 'coh' 'crossspecimag' }
                if strcmpi(g.measure, 'coh')
                    PS = abs(S.COH); % do not know what to do here
                    butterflyplot = squeeze(mean(PS, 2));
                    PSarea2area = squeeze(mean(PS(frq_inds, :, :)));
                    cortexPlot = mean(PSarea2area, 2);
                elseif strcmpi(g.measure, 'crossspecimag')
                    PS = abs(imag(cs2coh(S.CS)));
                    butterflyplot = squeeze(mean(PS, 2));
                    PSarea2area = squeeze(mean(PS(frq_inds, :, :)));
                    cortexPlot = mean(PSarea2area, 2);
                else
                    PS = cs2psd(S.CS);
                    butterflyplot = squeeze(mean(PS, 2));
                    apow = squeeze(sum(sum(reshape(PS(frq_inds, :), [], S.nROI), 1), 2)).*S.source_roi_power_norm';
                    cortexPlot = 10*log10(apow);
                    PSarea2area = [];
                end
                plotPSDFreq = S.freqs(frq_inds);
                plotPSD     = PS(frq_inds, :);
                matrix      = PSarea2area;

            case {'pac'}
                if isfield(S.PAC, 'b_orig_norm')
                    matrix = S.PAC.b_orig_norm;
                elseif isfield(S.PAC, 'b_orig')
                    matrix = S.PAC.b_orig;
                else
                    error('PAC (original bicoherence) cannot be plotted, field is missing.')
                end
                cortexPlot = mean(matrix, 2);

            case {'pac_anti'}
                if isfield(S.PAC, 'b_anti_norm')
                    matrix = S.PAC.b_anti_norm;
                elseif isfield(S.PAC, 'b_anti')
                    matrix = S.PAC.b_anti;
                else
                    error('PAC (antisymmetrized bicoherence) cannot be plotted, field is missing.')
                end
                cortexPlot = mean(matrix, 2);
                
            case {'tde' 'tde_anti'}
                if strcmpi(g.measure, 'tde')
                    t_xy = S.TDE.T_XY;
                    t_yx = S.TDE.T_YX;
                else
                    t_xy = S.TDE.aT_XY;
                    t_yx = S.TDE.aT_YX;
                end
                plot_tde(t_xy, S.TDE.shift, S.TDE.region_X, S.TDE.region_Y, S.TDE.method)
                plot_tde(t_yx, S.TDE.shift, S.TDE.region_Y, S.TDE.region_X, S.TDE.method)

                % disable cortex plot
                g.plotcortex = 'off';
        end

        % get seed
        if ~isempty(g.plotcortexseedregion)
            [coordinate, seed_idx] = get_seedregion_coordinate(EEG.roi.atlas.Scouts, g.plotcortexseedregion, EEG.roi.cortex.Vertices);
            seedMask = zeros(size(matrix));
            seedMask(seed_idx, :) = 1;
            seedMask(:,seed_idx) = 1;
        else
            seedMask = ones(size(matrix));
        end

        % frequency plot
        if strcmpi(g.plotpsd, 'on')
            figure; semilogy(plotPSDFreq, plotPSD); grid on
            h = textsc(plotOpt.label, 'title');
            set(h, 'fontsize', 20);
        end

        % plot on matrix
        if strcmpi(g.plotmatrix, 'on') && ~isempty(matrix)
            matrix = matrix.*seedMask;
            roi_plotcoloredlobes(EEG, matrix, titleStr, g.measure, g.hemisphere, g.grouphemispheres, g.region);
%             try
%                 roi_plotcoloredlobes(EEG, matrix, titleStr, g.measure, g.hemisphere, g.grouphemispheres, g.region);
%             catch
%                warning('Functionalities only available for the Desikan-Killiany atlas (68 ROIs).')
%                figure; imagesc(matrix);
%             end
        end

        % plot on cortical surface
        if strcmpi(g.plotcortex, 'on') && cortexFlag ~= -1
            if ~strcmpi(g.measure, 'gc') || ~strcmpi(g.measure, 'trgc')
                cortexTitle = [ plotOpt.labelshort ' (' titleStr ')' ];
            end
            if isempty(g.plotcortexseedregion)
                allplots_cortex_BS(cortex_highres, cortexPlot, [min(cortexPlot) max(cortexPlot)], cm17a, g.measure, g.smooth);
%                 allplots_cortex_BS(cortex_highres, cortexPlot, [min(cortexPlot) max(cortexPlot)], cm17a, upper(splot.unit), g.smooth);
%                 allplots_cortex_BS(S.cortex, cortexPlot, [min(cortexPlot) max(cortexPlot)], cm17a, plotOpt.unit, g.smooth);
            else
                cortexTitle = [ cortexTitle ' for area ' int2str(seed_idx)];
                cortexPlot = squeeze(matrix(seed_idx,:));
                allplots_cortex_BS(S.cortex, cortexPlot, [min(cortexPlot) max(cortexPlot)], cm17a, g.measure, g.smooth, [], {coordinate});
%                 allplots_cortex_BS(cortex_highres, cortexPlot, [min(cortexPlot) max(cortexPlot)], cm17a, upper(g.measure), g.smooth, [], {coordinate});
%                 allplots_cortex_BS(S.cortex, cortexPlot, [min(cortexPlot) max(cortexPlot)], cm17a, plotOpt.unit, g.smooth, [], {coordinate});
            end
            h = textsc(cortexTitle, 'title');
            set(h, 'fontsize', 20);
        elseif strcmpi(g.plotcortex, 'on') && cortexFlag == -1
            warning('EEG.roi.cortex does not contain the field "Faces" required to plot surface topographies.')
        end

        % make butterfly plot
        if strcmpi(g.plotbutterfly, 'on') && ~isempty(matrix)
            if strcmpi(g.measure, 'pac') || strcmpi(g.measure, 'pac_anti')
                warning('Butterfly plots (frequency x connectivity) cannot be computed for PAC because frequencies have already been specified for the computation.')
            else
                figure; plot(EEG.roi.freqs, butterflyplot, 'LineWidth', 1)
                h = title([ 'ROI to ROI ' replace_underscores(g.measure) ' (' titleStr ')' ]);
                set(h, 'fontsize', 16);
                xlabel('Frequency (Hz)')
                ylabel([replace_underscores(g.measure) ' (a.u.)'])
                grid on
            end

        end

        % plot 3D
        if strcmpi(g.plot3d, 'on') && ~isempty(matrix)
            PSarea2area = matrix.*seedMask;
            roi_plotbrainmovie(matrix, 'cortex', EEG.roi.cortex, 'atlas', EEG.roi.atlas, g.plot3dparams{:});
        end
    end

    if nargin < 2
        com = sprintf('pop_roi_connectplot(EEG, %s);', vararg2str( options ));
    end
end

function [coordinate, seed_idx] = get_seedregion_coordinate(scouts, seed_idx, vc)
    % determine voxel of selected seed region, if needed
    % assign region index to selected seed region (passed as index)
    if ~isempty(seed_idx)
        % ball not visible for these regions when plotting the mean voxel
        manual_region_idxs = [2, 16, 18, 25, 26, 31, 32, 45, 49, 50, 55, 56, 59, 60, 61, 64]; 
        pos_idx = scouts(seed_idx).Vertices;
        pos = vc(pos_idx,:);
        if seed_idx == 1
            coordinate = vc(736,:);
        elseif ismember(seed_idx, manual_region_idxs)
            pos_sorted = sortrows(pos, 3, 'descend'); % sort by descending Z-coordinate
            coordinate = pos_sorted(1,:);
        else
            mid_point = mean(pos,1);
            [~,closest_pos_idx] = min(eucl(mid_point, pos)); % determine mean voxel
            coordinate = pos(closest_pos_idx,:);
        end
    else
        error('Selected region not in cortex')
    end
end

% function labels = get_labels(EEG)
%     % retrieve labels from atlas
%     labels = strings(1,length(EEG.roi.atlas.Scouts));
%     for i = 1:length(labels)
%         scout = struct2cell(EEG.roi.atlas.Scouts(i));
%         labels(i) = char(scout(1));
%     end
%     labels = cellstr(labels);
% 
%     % remove region labels that were not selected
%     if isfield(EEG.roi, 'roi_selection')
%         if ~isempty(EEG.roi.roi_selection)
%             labels = labels(cell2mat(EEG.roi.roi_selection));
%         end
%     end
% end
% 
% function new_labels = replace_underscores(labels)
%     % remove underscores in label names to avoid bug
%     new_labels = strrep(labels, '_', ' ');
% end

function [colors, color_idxx, roi_idxx, labels_sorted, roi_loc] = get_colored_labels(EEG)
    labels = get_labels(EEG);

    colors = {[0,0,0]/255, [163, 107, 64]/255, [171, 163, 71]/255, [217, 37, 88]/255, [113, 15, 82]/255,[35, 103, 81]/255,[2, 45, 126]/255,};
    % assign labels to colors
    roi_loc ={'LT';'RT';'LL';'RL';'LF';'RF';'LO';'RO';'LT';'RT';'LPF';'RPF';'LT';'RT';'LP';'RP';'LT';'RT';'LT';'RT';'LL';'RL';'LO';'RO';'LPF';'RPF';'LO';'RO';'LPF';'RPF';'LT';'RT';'LC';'RC';'LT';'RT';'LF';'RF';'LPF';'RPF';'LF';'RF';'LO';'RO';'LC';'RC';'LL';'RL';'LC';'RC';'LP';'RP';'LL';'RL';'LF';'RF';'LF';'RF';'LP';'RP';'LT';'RT';'LP';'RP';'LT';'RT';'LT';'RT'};
    roi_loc = string(roi_loc);
    roi_loc = strrep(roi_loc, 'PF', '2');
    roi_loc = strrep(roi_loc, 'F', '3');
    roi_loc = strrep(roi_loc, 'T', '4');
    roi_loc = strrep(roi_loc, 'P', '5');
    roi_loc = strrep(roi_loc, 'C', '6');
    roi_loc = strrep(roi_loc, 'O', '7');
    roi_loc = strrep(roi_loc, 'LL', 'L1');
    roi_loc = strrep(roi_loc, 'RL', 'R1');
    roi_loc = strrep(roi_loc, 'L', '');
    roi_loc = strrep(roi_loc, 'R', '');

    % remove regions that were not selected
    if isfield(EEG.roi, 'roi_selection')
        if ~isempty(EEG.roi.roi_selection)
            roi_loc = roi_loc(cell2mat(EEG.roi.roi_selection));
        end
    end

    try
        [color_idxx,roi_idxx] = sort(str2double(roi_loc));
        labels_sorted = labels(roi_idxx);
    catch
        roi_idxx = 1:length(labels);
        color_idxx = mod(roi_idxx, length(colors))+1;
        labels_sorted = labels;
    end
end

function roi_plotpower(EEG, source_roi_power_norm_dB, titleStr)
    [colors, color_idxx, roi_idxx, labels_sorted, ~] = get_colored_labels(EEG);
    n_roi_labels = size(labels_sorted, 2);

    barh(source_roi_power_norm_dB(roi_idxx));
    set(gca, 'YDir', 'reverse');
    set(gca,'ytick',[1:n_roi_labels],'yticklabel',labels_sorted(1:end), 'fontweight','bold','fontsize', 9, 'TickLength',[0.015, 0.02], 'LineWidth',0.7);
    h = title([ 'ROI source power' ' (' titleStr ')' ]);
    set(h, 'fontsize', 16);
    ylabel('power [dB]')
    ax = gca;
    for i=1:numel(roi_idxx)
        ax.YTickLabel{ceil(i)} = sprintf('\\color[rgb]{%f,%f,%f}%s', colors{color_idxx(i)}, ax.YTickLabel{ceil(i)});
    end
    pos = get(gcf, 'Position');
    set(gcf, 'Position', [pos(1) pos(2) pos(3) pos(4)*1.8])
    movegui(gcf, 'south') % remove after
end
        
function roi_plotcoloredlobes( EEG, matrix, titleStr, measure, hemisphere, grouphems, region)
    % check if Desikan-Killiany atlas is used
    if EEG.roi.nROI == 68
        isDKatlas = true;
    else
        isDKatlas = false;
    end
    
    if ~strcmpi(region, 'all') && isDKatlas == 0
        error('Region plotting is only supported for the Desikan-Killiany atlas.');
    end
    
    % plot matrix with colored labels    
    load cm18
    switch lower(measure)
        case {'mim', 'mic', 'coh'}
            cmap = cm18a;
        otherwise
            cmap = cm18;
    end
    clim_min = min(matrix, [], 'all');
    clim_max = max(matrix, [], 'all');

    % hemisphere parameters to determine which labels to use 
    last_char = EEG.roi.atlas.Scouts(1).Label(end); 
    if strcmpi(hemisphere, 'left')
        if strcmpi(last_char, 'R')
            hem_idx = {2 2 2};  % use labels 2:2:end (first two values), only use 1/2 of the labels (3rd value)
        else 
            hem_idx = {1 2 2};  % use labels 1:2:end (first two values), only use 1/2 of the labels (3rd value)
        end
    elseif strcmpi(hemisphere, 'right')
        if strcmpi(last_char, 'L')
            hem_idx = {1 2 2};  
        else
            hem_idx = {2 2 2};  
        end
    else
        hem_idx = {1 1 1};
    end
    
    % sort matrix according to color scheme
    % reduce matrix to only keep components corresponding to selected region
    if isDKatlas == 1
        [colors, color_idxx, roi_idxx, labels_sorted, roi_loc] = get_colored_labels(EEG);
        labels = labels_sorted;

        % assign region input to an index
        [GC, ~] = groupcounts(roi_loc);
        switch lower(region)
            case 'cingulate'
                region_idx = 1;
            case 'prefrontal'
                region_idx = 2;
            case 'frontal'
                region_idx = 3;
            case 'temporal'
                region_idx = 4;
            case 'parietal'
                region_idx = 5;
            case 'central'
                region_idx = 6;
            case 'occipital'
                region_idx = 7;
            otherwise
                region_idx = 99;
        end

        matrix = matrix(roi_idxx, roi_idxx); 
        if not(region_idx == 99)
            if region_idx == 1
                start_idx = 1;
            else
                start_idx = 1 + sum(GC(1:region_idx-1));
            end
            end_idx = start_idx + GC(region_idx) - 1;
            matrix = matrix(start_idx:end_idx, start_idx:end_idx);  
            labels = labels(start_idx:end_idx);
            color_idxx = color_idxx(start_idx:end_idx);
        end
    else
        labels = get_labels(EEG);
    end
    n_roi_labels = size(matrix, 1);

    % remove underscores in labels to avoid plotting bug
    labels = replace_underscores(labels);

    % create dummy plot and add custom legend
    f = figure();
    %f.WindowState = 'maximized';
    hold on
    n_dummy_labels = 7;
    x = 1:10;

    % labels on dummy plot for positioning
    xlim([0 n_roi_labels])
    ylim([0 n_roi_labels])
    set(gca,'xtick',1:n_roi_labels,'xticklabel',labels(hem_idx{1}:hem_idx{2}:n_roi_labels));%, 'TickLabelInterpreter','none');
    ax = gca;
    if isDKatlas == 1
        for k=1:n_dummy_labels
            plot(x, x*k, '-', 'LineWidth', 9, 'Color', colors{k});
        end
        for i=hem_idx{1}:hem_idx{2}:n_roi_labels   
            ax.XTickLabel{ceil(i/hem_idx{3})} = sprintf('\\color[rgb]{%f,%f,%f}%s', colors{color_idxx(i)}, ax.XTickLabel{ceil(i/2)});
        end
        legend('Cingulate', 'Prefrontal', 'Frontal', 'Temporal', 'Parietal', 'Central', 'Occipital', 'Location', 'southeastoutside'); % modify legend position
    end
    xtickangle(90)
    pos = get(gca, 'Position');
    set(gca, 'Position', pos, 'DataAspectRatio',[1 1 1], 'visible', 'off')
    axes('pos', [pos(1) pos(2) pos(3) pos(4)]) % plot matrix over the dummy plot and keep the legend
    
    % group by hemispheres (left first, then right)
    if strcmp(grouphems, 'on')
        % sort matrix
        mat_left_row = matrix(1:2:end,:);
        mat_right_row = matrix(2:2:end,:);
        matrix = [mat_left_row; mat_right_row]; % sort rows
        mat_left_col = matrix(:,1:2:end);
        mat_right_col = matrix(:,2:2:end);
        matrix = horzcat(mat_left_col, mat_right_col); % sort columns
        
        % sort labels
        try % if color can be assigned
            lc = {labels; transpose(color_idxx)};
            for i = 1:length(lc)
                left = lc{i}(1:2:end);
                right = lc{i}(2:2:end);
                lc{i} = [left right];
            end
            labels = lc{1};
            color_idxx = transpose(lc{2});
        catch
            left = labels(1:2:end);
            right = labels(2:2:end);
            labels = [left right];
        end
    end

    % reduce matrix to keep only one hemisphere
    if strcmp(hemisphere, 'left')
        matrix = matrix(hem_idx{1}:hem_idx{2}:n_roi_labels,:);
        matrix(:,2:2:end) = []; 
        imagesc(matrix); colormap(cmap);  
    elseif strcmp(hemisphere, 'right')
        matrix = matrix(hem_idx{1}:hem_idx{2}:n_roi_labels,:);
        matrix(:,1:2:end) = [];
        imagesc(matrix); colormap(cmap);
    else
        imagesc(matrix); colormap(cmap);
    end
    cb = colorbar;
%     tf = isMATLABReleaseOlderThan("R2022a");
%     if tf
%         caxis([clim_min clim_max])
%     else
%         clim([clim_min clim_max])
%     end
    try 
        caxis([clim_min clim_max])
    catch
        clim([clim_min clim_max])
    end
    if isDKatlas == 1
        set(cb, 'Location', 'southoutside')
    else
        set(cb, 'Location', 'eastoutside')
    end
    set(gca, 'Position', pos, 'DataAspectRatio',[1 1 1], 'visible', 'on')

    % add colored labels with display option
    ax = gca;
    set(gca,'xtick',1:n_roi_labels,'xticklabel',labels(hem_idx{1}:hem_idx{2}:n_roi_labels));
    if isDKatlas == 1
        set(gca,'ytick',1:n_roi_labels,'yticklabel',labels(hem_idx{1}:hem_idx{2}:n_roi_labels), 'fontsize', 9, 'TickLength',[0.015, 0.02], 'LineWidth',0.75);
        for i=hem_idx{1}:hem_idx{2}:n_roi_labels  
            ax.XTickLabel{ceil(i/hem_idx{3})} = sprintf('\\color[rgb]{%f,%f,%f}%s', colors{color_idxx(i)}, ax.XTickLabel{ceil(i/hem_idx{3})});
            ax.YTickLabel{ceil(i/hem_idx{3})} = sprintf('\\color[rgb]{%f,%f,%f}%s', colors{color_idxx(i)}, ax.YTickLabel{ceil(i/hem_idx{3})});
        end
    else
        set(gca,'ytick',1:n_roi_labels,'yticklabel',labels(hem_idx{1}:hem_idx{2}:n_roi_labels), 'fontsize', 7, 'TickLength',[0.015, 0.02], 'LineWidth',0.75);
    end
    set(gca, 'YDir','normal')
    h = title([ 'ROI to ROI ' replace_underscores(measure) ' (' titleStr ')' ]);
    set(h, 'fontsize', 16);
    xtickangle(90)
end 

function roi_largeplot(EEG, mim, trgc, roipsd, titleStr)
    % plot MIM, TRGC and power (barplot) in a single large figure
    load cm18
    [colors, color_idxx, roi_idxx, labels_sorted, ~] = get_colored_labels(EEG);
    n_roi_labels = size(labels_sorted, 2);
    
    f = figure();
    f.WindowState = 'maximized';
    fc_matrices = cell(1,2);
    fc_matrices{1,1} = mim;
    fc_matrices{1,2} = trgc;
    fc_names = ["MIM", "TRGC"];
    
    for k = 1:2
        plt(k) = subplot(1,3,k);
        fc = fc_matrices{k};
        img = squeeze(fc)';
        img_sorted = img(roi_idxx, roi_idxx);
        imagesc(img_sorted)

        set(gca,'ytick',[1:n_roi_labels],'yticklabel',labels_sorted(1:end), 'fontsize', 5, 'TickLength',[0.015, 0.02], 'LineWidth',0.75);
        set(gca,'xtick',[1:n_roi_labels],'xticklabel',labels_sorted(1:end));
        h = title([ 'ROI to ROI ' fc_names{k} ' (' titleStr ')' ]);
        set(h, 'fontsize', 16);
        hcb = colorbar;
        hcb.Label.FontSize = 10;
        set(gca,'DataAspectRatio',[1 1 1])
        xtickangle(90)
        ax = gca;
        for i=1:numel(roi_idxx)    
            ax.XTickLabel{ceil(i)} = sprintf('\\color[rgb]{%f,%f,%f}%s', colors{color_idxx(i)}, ax.XTickLabel{ceil(i)});
            ax.YTickLabel{ceil(i)} = sprintf('\\color[rgb]{%f,%f,%f}%s', colors{color_idxx(i)}, ax.YTickLabel{ceil(i)});
        end
    end
    colormap(plt(1), cm18a)
    colormap(plt(2), cm18)
    
    % power
    subplot(1,3,3);
    barh(roipsd(roi_idxx));

    set(gca, 'YDir', 'reverse');
    set(gca,'ytick',[1:n_roi_labels],'yticklabel',labels_sorted(1:end), 'fontsize', 9, 'TickLength',[0.015, 0.02], 'LineWidth',0.7);
    h = title([ 'ROI source power' ' (' titleStr ')' ]);
    set(h, 'fontsize', 16);
    ylabel('power [dB]')
    ax = gca;
    for i=1:numel(roi_idxx)
        ax.YTickLabel{ceil(i)} = sprintf('\\color[rgb]{%f,%f,%f}%s', colors{color_idxx(i)}, ax.YTickLabel{ceil(i)});
    end
    
    % invisible plot to add legend
    hold on
    n_labels = 7;
    h = zeros(n_labels, 1);
    for k=1:numel(h)
        h(k) = plot(NaN, NaN, '-', 'LineWidth', 8, 'Color', colors{k});
    end
    if n_roi_labels == 68
        lgd = legend(h, 'Cingulate', 'Prefrontal', 'Frontal', 'Temporal', 'Parietal', 'Central', 'Occipital');
        lgd.FontSize = 10;
        set(lgd, 'Position', [0.44 0.06 0.25 0.25]);
    end
end

function plot_tde(T, shift, region_X, region_Y, method)
    [~, peak_idx] = max(T); 
    est_delay = shift(peak_idx); % in Hz
    
    figure; plot(shift, T, 'LineWidth', 1)
    xline(est_delay, '--r')
    xlabel('Time (s)')
    ylabel('a.u.')
    h = title(sprintf('%s -> %s TDE | Method %d', region_X, region_Y, method));
    set(h, 'fontsize', 16);
    subtitle("\tau = " + num2str(est_delay) + " (s)")
    grid on
end
