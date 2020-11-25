% roi_network() - Compute ROI network activity
%
% Usage:
% results = roi_network( roiData, roiInds, 'key', val)
%
% Inputs:
%   ROIdata - [roi x times x trials] ROI activity. Assumes 1 row per ROI
%   ROIinds - indices of ROI in the network
%
% Optional spectral inputs:
%   'window'  - [integer] lenght of window in samples
%   'overlap' - [integer] number of overlaping samples
%   'nfft'    - [integer] number of sample for FFT 
%
% Optional inputs:
%   'srate'      - [real] sampling rate. Default is 128.
%   'freqranges' - [cell] cell array of requency ranges. For example
%                  for theta and alpha { [3 5] [9 11] }

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

function results = roi_network( spatiallyFilteredData, ROI_inds, varargin)

if nargin < 2
    help roi_network;
    return;
end

g = finputcheck(varargin, { ...
    'nfft'        'real'      { }             [];
    'overlap'     'real'      { }             0;
    'window'      'real'      { }             128;
    'srate'       'real'      { }             128;
    'postprocess' 'struct'    {}              struct([]);
    'freqranges'  'cell'      {}              {};
    }, 'roi_network');
if isstr(g)
    error(g);
end

networkData = spatiallyFilteredData(ROI_inds,:);
[S,freqs] = cpsd_welch(networkData,g.window,g.overlap,g.nfft,g.srate);
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
connectSpecSelect = zeros(size(absiCOH,1), size(absiCOH,2), length(g.freqranges));
for iSpec = 1:length(g.freqranges)
    freqRangeTmp = intersect( find(freqs >= g.freqranges{iSpec}(1)), find(freqs <= g.freqranges{iSpec}(end)) );
    connectSpecSelect(:,:,iSpec) = mean(absiCOH(:,:,freqRangeTmp),3); % mean power in frequency range
end

if ~isempty(g.postprocess)
    connectprocess = g.postprocess;
    connectprocessFields = fieldnames(connectprocess);
    for iProcess = 1:length(connectprocessFields)
        results.([ connectprocessFields{iProcess} ]) = feval(connectprocess.(connectprocessFields{iProcess}), connectSpecSelect);
    end
else
    results = connectSpecSelect;
end

% -----------------------------

function [S,freqs] = cpsd_welch(X,window,noverlap, nfft, srate)

if isempty(nfft) || nfft < window, nfft = window; end
h = nfft/2+1;
n = size(X,1);
S = complex(zeros(n,n,h));
for i = 1:n
    [S(i,i,:), freqs] = pwelch(X(i,:),window,noverlap,nfft, srate);          % auto-spectra
    for j = i+1:n % so we don't compute cross-spectra twice
        S(i,j,:) = cpsd(X(i,:),X(j,:),window,noverlap,nfft,srate); % cross-spectra
    end
end
S = S/pi; % the 'pi' is for compatibility with 'autocov_to_cpsd' routine

