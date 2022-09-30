% roi_definenetwork() - define network between brain areas. Requires to
%                       have an Atlas loaded in EEG.roi.atlas
% Usage:
%  [EEG, networks] = roi_definenetwork(EEG, netTable, 'key', 'val'); 
%
% Inputs:
%  EEG - EEGLAB dataset with ROI activity computed and Atlas loaded
%  netTable      - [string|table] Define network based on existing
%                  ROIs. If a string is provided, the file is loaded as a
%                  table. First row contains the names of the new ROI.
%                  other rows contain the name of the areas to group (see
%                  example). Networks are defined as groups of ROIs.
%
% Optional input:
%   'addrois'       - [string|table] Define additional ROIs based on existing
%                     ROIs (see example).
%   'ignoremissing' - ['on'|'off'] ignore missing names 'on' or issue
%                     an error 'off'. Default is 'off'.
%
% Output:
%   EEG - EEG structure with EEG.roi.atlas.Scout and EEG.roi.atlas.networks
%         field updated and now containing new ROI or network.
%   networks - Same as EEG.roi.atlas.networks
%
% Example:
%   DNM = [1 2 3 4 5]';
%   EEG = roi_definenetwork(EEG, table(DNM)); % define network DNM comprising ROI 1, 2, 3, 4 and 5
%
% Example:
%   DNM = { 'Brodmann area 10L' 'Brodmann area 10R' }';
%   EEG = roi_definenetwork(EEG, table(DNM)); % define network DNM comprising ROI name 24L and 24R
%
% Example:
%   A = { 'Brodmann area 10L' 'Brodmann area 10R' }';
%   A = { 'Brodmann area 31L' 'Brodmann area 31R' }';
%   DNM = { 'A' 'B' }';
%   EEG = roi_definenetwork(EEG, table(DNM), 'addrois', table(A,B)); % define network DNM comprising ROI A and B
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

function [EEG,networks] = roi_definenetwork(EEG, roiTable, varargin)

if nargin < 2
    help roi_definenetwork;
    return
end

g = finputcheck(varargin, { ...
    'ignoremissing'    'string'    {'on' 'off'}    'off';
    'addrois'          ''          {}              [];
    }, 'roi_definenetwork');
if isstr(g)
    error(g);
end

if ischar(roiTable)
    roiTable = readtable(roiTable);
end
if ischar(g.addrois) && ~isempty(g.addrois)
    g.addrois = readtable(g.addrois);
end

try
    allLabels = { EEG.roi.atlas.Scouts.Label };
catch
    error('Atlas not found. Use pop_leadfield to choose a source model which contains an Atlas.');
end

% add new ROIs
if ~isempty(g.addrois)
    colNames = fieldnames(g.addrois);
    for iCol = 1:size(g.addrois,2) % scan columns
        
        EEG.roi.atlas.Scouts(end+1).Label = colNames{iCol};
        
        inds = [];
        if isnumeric(g.addrois(1,iCol))
            inds = g.addrois(:,iCol);
        else
            for iRow = 1:size(g.addrois,1)
                val = g.addrois{iRow, iCol}{1};
                if ~isempty(val)
                    indTmp = strmatch(val, allLabels, 'exact');
                    if length(indTmp) ~= 1
                        if strcmpi(g.ignoremissing, 'off')
                            error('Area %s not found or duplicate', val);
                        else
                            fprintf('Area %s not found or duplicate, using the first one\n', val);
                            indTmp = [];
                        end
                    else
                        inds = [ inds;indTmp ];
                    end
                end
            end
        end
        EEG.roi.atlas.Scouts(end).Vertices = vertcat(EEG.roi.atlas.Scouts(inds).Vertices);
    end
end
            
% define networks
networks = [];
colNames = fieldnames(roiTable);
allLabels = { EEG.roi.atlas.Scouts.Label };
for iCol = 1:size(roiTable,2) % scan columns
    
    networks(end+1).name = colNames{iCol};

    inds = [];
    if isnumeric(roiTable{1,iCol})
        inds = [ roiTable{:,iCol} ]';
    else
        for iRow = 1:size(roiTable,1)
            val = roiTable{iRow, iCol}{1};
            if ~isempty(val)
                indTmp = strmatch(val, allLabels, 'exact');
                if length(indTmp) ~= 1
                    if strcmpi(g.ignoremissing, 'off')
                        error('Area %s not found or duplicate', val);
                    else
                        fprintf('Area %s not found or duplicate, using the first one\n', val);
                        if length(indTmp) > 1 indTmp = indTmp(1); end
                    end
                end
                inds = [ inds;indTmp ];
            end
        end
    end
    networks(end).ROI_inds = inds';
    
end

EEG.roi.atlas.networks = networks;
