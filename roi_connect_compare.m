% roi_connect_compare - Compare ROI in collection of datasets.
%
% Usage:
%  roi_connect_compare(ALLEEG, datCond, 'key', 'val', ...);
%
% Inputs:
%  ALLEEG - EEGLAB datasets
%  datCond - [cell] list of datasets to compare { [1:10] [11:20 } to
%            compare datasets 1-10 with 11:20
%
% Optional inputs:
%
% Output:
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

function [EEG] = roi_connect_compare(EEG, datCond, varargin)

if nargin < 2
    help roi_activity;
    return
end

% decode input parameters
% -----------------------
[g,addopts] = finputcheck(varargin,  { 'measure'    'string'  {}  'MIM';
    'freqrange'             {'real' 'cell' }     []                    [];
    'quantile'              'real'     { }                     [];
    'plotsubjmat'           'string'   { 'on' 'off' }          'off'; ...
    'threshold'             'real'     { }                     1.1;
    'outputdir'             'string'   { }                     '.';
    'network'               'string'   { 'on' 'off' }          'off'; ...
    'addrois'               'string'   { 'on' 'off' }          'off'; ...
    'pmask'                 'string'   { 'on' 'off' }          'off' }, 'pop_roi_connectplot', 'ignore');
if ischar(g), error(g); end
if isempty(g.freqrange)
    error('You must define freqrange');
end
if ~iscell(g.freqrange)
    g.freqrange = { g.freqrange };
end

% select data
% -----------
dataCell = cell(1, length(datCond));
freqs    = EEG(1).roi.freqs;
for iFreq = 1:length(g.freqrange)
    for iCond = 1:length(datCond)
        for iSet = 1:length(datCond{iCond}) % can only index datasets
            setNum = datCond{iCond}(iSet);
            dataCell{iFreq, iCond}(:,:,iSet) = selectconnect(EEG(setNum).roi.(g.measure), freqs, g.freqrange{iFreq}, g.quantile);
        end
    end
end

% compute stats
% -------------
freqTextSig = cell(1, length(g.freqrange));
freqTextSig(:) = { '' };
if strcmpi(g.pmask, 'on')
    freqTextSig = {};
    for iFreq = 1:length(g.freqrange)
        pcond = std_stat(dataCell(iFreq,:)', 'condstats', 'on', 'method', 'permutation', 'naccu', 50, 'mcorrect', 'fdr');
        %pcond = std_stat(dataCell', 'condstats', 'on', 'method', 'permutation', 'naccu', 200); %, 'mcorrect', 'fdr'); pcond{1} = pcond{1} < 0.5;
        pcond = pcond{1};
        
        for iDiag = 1:size(pcond,1), pcond(iDiag, iDiag) = 1; end
        if all(pcond(:) == 1)
            fprintf('Frequency %1.1f to %1.1f -> Nothing significant\n', g.freqrange{iFreq}(1), g.freqrange{iFreq}(2));
            freqTextSig{iFreq} = 'ns';
        else
            fprintf(2, 'Frequency %1.1f to %1.1f -> Something significant\n', g.freqrange{iFreq}(1), g.freqrange{iFreq}(2));
            freqTextSig{iFreq} = 'significant';
            for iCond = 1:length(datCond(:))
                dataCell{iFreq, iCond} = bsxfun(@times,dataCell{iFreq, iCond}, pcond<1 );
            end
        end
    end
end

% plot matrices
% -------------
if strcmpi(g.plotsubjmat, 'on')
    for iFreq = 1:length(g.freqrange)
        figure('position', [275  1075 1725 262]);
        numDat = size(dataCell{1},3);
        numRows = length(dataCell(:));
        cl = [0.03 0.1];
        %cl = [0.0 1];
        for iCond = 1:numRows
            for iDat = 1:numDat
                subplot(numRows, numDat, iDat+numDat*(iCond-1));
                imagesc(dataCell{iCond}(:,:,iDat));
                clim(cl);
                axis off;
            end
        end
        textsc('title', sprintf('Connectivity %g-%g Hz %s', g.freqrange{iFreq}(1), g.freqrange{iFreq}(2)), freqTextSig{iFreq});
    end
end

% plot surfaces
% -------------
if isfield(EEG(1).roi.cortex, 'Faces')
    for iFreq = 1:length(g.freqrange)
        mn = inf;
        mx = -inf;
        for iCond = 1:length(datCond)
            tmpData = mean(dataCell{iFreq, iCond},3);
            mn = min(mn, min(tmpData(:)));
            mx = max(mx, max(tmpData(:)));
        end
        load -mat cm17;
        for iCond = 1:length(datCond)
            tmpData = mean(dataCell{iFreq, iCond},3);
            figure; allplots_cortex_BS(EEG(1).roi.cortex, tmpData, [mn mx*0.5], cm17a, '', 0.35);
            textsc('title', sprintf('Connectivity %g-%g Hz condition %d', freqrange{iFreq}(1), freqrange{iFreq}(2)), iCond);
        end
    end
end

% define connectivity regions
% ---------------------------
if strcmpi(g.network, 'on')
    tab1 = readtable('NGNetworkROIs_v4.txt','delimiter', char(9));
    for iCond = 1:length(dataCell(:))
        if strcmpi(g.addrois, 'on')
            tab2 = readtable('NGNetworkROIs_area_definition_v2.txt','delimiter', char(9));
            [EEG(1),net,dataCell{iCond}] = roi_definenetwork(EEG(1), tab1, 'addrois', tab2, 'connectmat', dataCell{iCond}, 'ignoremissing', 'on'); % adding missing ROIs
        else
            [EEG(1),net,dataCell{iCond}] = roi_definenetwork(EEG(1), tab1, 'connectmat', dataCell{iCond}, 'ignoremissing', 'on'); % adding missing ROIs
        end
    end
    if size(dataCell,2) > 1
        for iNetwork = 1:length(net)
            netTmp = net([iNetwork iNetwork iNetwork]);
            netTmp(1).name = [ netTmp(1).name 'Condition 1' ];
            netTmp(2).name = [ netTmp(2).name 'Condition 2' ]; 
            netTmp(2).name = [ netTmp(2).name 'Difference' ];
            for iFreq = 1:length(g.freqrange)
    
                tit = sprintf('MIM %g_%g Hz', g.freqrange{iFreq}(1), g.freqrange{iFreq}(2));
    
                medthresh = nanmedian(dataCell{iFreq,1}(:))*g.threshold;
                roi_networkplot(EEG(1), netTmp, [ dataCell(iFreq,:) {dataCell{iFreq,2}-dataCell{iFreq,1}}], 'threshold', medthresh, 'subplots', 'on', 'title', tit, 'columns', 3, addopts{:}); %, 'limits'); %, [0.05 0.08]);
                tit(tit == ' ') = '_';
                print('-djpeg', [g.outputdir filesep 'Connectivity_maps_MIM_' net(iNetwork).name '_' tit '.jpg']);
                close;
            end
        end
    else
        for iFreq = 1:length(g.freqrange)
            tit = sprintf('MIM %g_%g Hz', g.freqrange{iFreq}(1), g.freqrange{iFreq}(2));
            roi_networkplot(EEG(1), net, dataCell{iFreq}, 'threshold', 0, 'subplots', 'on', 'title', tit, 'limits', [0.1 0.2]);
            tit(tit == ' ') = '_';
            print('-djpeg', [g.outputdir filesep 'Connectivity_maps_' tit '_fixed.jpg']);
            close;
        end
    end
end

% select connection based on frequencies and freq range
% -----------------------------------------------------
function MIMtmp = selectconnect(CON, freqs, range, quant)

[~,fBeg] = min(abs(freqs-range(1)));
[~,fEnd] = min(abs(freqs-range(2)));
MIMtmp = mean(CON(fBeg:fEnd,:,:),1);

if ~isempty(quant)
    TMPCON = MIMtmp(:);
    TMPCON(TMPCON == 0) = [];
    
    if quant < 0.5
        quant = 1-quant;
    end
    mn = quantile(TMPCON, quant);
    mx = quantile(TMPCON, 1-quant);

    MIMtmp = (MIMtmp-mn)/(mx-mn);
end

