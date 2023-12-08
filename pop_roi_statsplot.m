% pop_roi_statsplot() - Generate p-values from FC null distributions and plots them. Based on Franziska Pellegrini's script
% fp_plot_fc_shuffletest.
%
% Inputs:
%   EEG        -  EEGLAB dataset with ROI activity computed.
%
% Optional inputs:
%   'measure'  - [cell] Cell of strings corresponding to methods.
%                       'CS'    : Cross spectrum
%                       'aCOH'  : Coherence
%                       'cCOH'  : (complex-valued) coherency
%                       'iCOH'  : absolute value of the imaginary part of coherency
%                       'wPLI'  : Weighted Phase Lag Index
%                       'PDC'   : Partial directed coherence
%                       'TRPDC' : Time-reversed partial directed coherence
%                       'DTF'   : Directed transfer entropy
%                       'TRDTF' : Time-reversed directed transfer entropy
%                       'MIM'   : Multivariate Interaction Measure for each ROI
%                       'MIC'   : Maximized Imaginary Coherency for each ROI
%  'freqrange' - [min max] frequency range or [integer] single frequency in Hz. Default is to plot broadband power.
%  'alpha'     - [integer] Significance level. Default is 0.05.
%
%   Author: Tien Dung Nguyen, tien-dung.nguyen@charite.de

function EEG = pop_roi_statsplot(EEG, varargin)

    if nargin < 2
        help roi_connstats;
        return
    end

    if ~isfield(EEG, 'roi') || ~isfield(EEG.roi, 'source_roi_data')
        error('Cannot find ROI data - compute ROI data first');
    end

    % decode input parameters
    % -----------------------
    g = finputcheck(varargin,  { 
        'measure'        'string'   { }    '';
        'freqrange'      'real'     { }    []; ...
        'alpha'          'integer'  { }    0.05}, 'pop_roi_statsplot');
    if ischar(g), error(g); end
    S = EEG.roi;

    if isempty(g.measure)
        error('You must define a measure to plot');
    end
    
    % extract frequency indices
    if ~isempty(g.freqrange)
        if length(g.freqrange) == 1
            frq_inds = find(S.freqs == g.freqrange(1)); 
            title = sprintf('%1.1f Hz', g.freqrange(1));
        else
            frq_inds = find(S.freqs >= g.freqrange(1) & S.freqs < g.freqrange(2));
            title = sprintf('%1.1f-%1.1f Hz frequency band', g.freqrange(1), g.freqrange(2));
        end
    else
        frq_inds = 1:length(S.freqs);
        title = 'broadband';
    end
    
    % select frequency or frequency band
    if length(frq_inds) > 1
        matrix = squeeze(mean(S.(g.measure)(frq_inds, :, :, :)));
    else
        matrix = squeeze(S.(g.measure)(frq_inds, :, :, :));
    end
    
    % generate p-values by comparing the true FC (first shuffle) to null distribution
    netFC = squeeze(mean(matrix, 2));
    FC_pn = sum(netFC(:, 1) < netFC(:, 2:end), 2)./(size(matrix, 3) - 1);

    % use FDR-correction for multiple comparison's correction
    [p_fdr, ~] = fdr(FC_pn, g.alpha);
    FC_pn(FC_pn > p_fdr) = 1;

    % plot 
    load cm17;
    load cortex; 
    FC_pn(FC_pn==0) = 1 / (size(netFC, 2) - 1);  % 1 / nshuf
    data = -log10(FC_pn);
    try
        allplots_cortex_BS(cortex_highres, data, [min(data) max(data)], cm17a ,'-log(p)', 0.3);
        h = textsc(title, 'title');
        set(h, 'fontsize', 20);
    catch
        warning('There are no "significant" p-values to be plotted.')
    end
end