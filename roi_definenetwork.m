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
%   'connectmat'    - [array] connectivity matrix. When provided return
%                     the new connectivity with added ROI ('addrois' input)
%
% Output:
%   EEG - EEG structure with EEG.roi.atlas.Scout and EEG.roi.atlas.networks
%         field updated and now containing new ROI or network.
%   networks   - Same as EEG.roi.atlas.networks
%   connectmat - Updated connectivity matrix (when provided as input)
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
% Example:
%  [EEG, net] = roi_definenetwork(EEG, 'NGNetworkROIs_v4.txt', 'addrois', 'NGNetworkROIs_area_definition_v2.txt', 'ignoremissing', 'on'); 
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

function [EEG,networks,connectmat] = roi_definenetwork(EEG, roiTable, varargin)

if nargin < 2
    help roi_definenetwork;
    return
end

g = finputcheck(varargin, { ...
    'ignoremissing'    'string'    {'on' 'off'}    'off';
    'addrois'          ''          {}              [];
    'connectmat'       ''          {}              [];
    }, 'roi_definenetwork');
if isstr(g)
    error(g);
end

if ischar(roiTable)
    if ~exist(roiTable, 'file')
        p = fileparts(which('roi_definenetwork'));
        roiTable2 = fullfile(p, roiTable);
        if ~exist(roiTable2, 'file')
            error('File not found %s', roiTable);
        end
        roiTable = roiTable2;
    end
    roiTable = readtable(roiTable,'Delimiter', char(9));
end
if ischar(g.addrois) && ~isempty(g.addrois)
    if ~exist(g.addrois, 'file')
        p = fileparts(which('roi_definenetwork'));
        tmpTable = fullfile(p, g.addrois);
        if ~exist(tmpTable, 'file')
            error('File not found %s', g.addrois);
        end
        g.addrois = tmpTable;
    end
    g.addrois = readtable(g.addrois,'delimiter', char(9));
end

try
    allLabels = { EEG.roi.atlas.Scouts.Label };
catch
    error('Atlas not found. Use pop_leadfield to choose a source model which contains an Atlas.');
end

% add new ROIs
connectmat = g.connectmat;
if ~isempty(g.addrois)
    colNames = fieldnames(g.addrois);
    ROIinds = cell(1, size(g.addrois,2));

    for iCol = 1:size(g.addrois,2) % scan columns
        EEG.roi.atlas.Scouts(end+1).Label = colNames{iCol};
        
        inds = [];
        if isnumeric(g.addrois(1,iCol))
            inds = g.addrois(:,iCol);
        else
            for iRow = 1:size(g.addrois,1)
                val = g.addrois{iRow, iCol}{1};
                if ~isempty(val)
                    indTmp1 = strmatch(val, allLabels, 'exact');
                    indTmp2 = strmatch([ 'Brodmann area ' val], allLabels, 'exact');
                    indTmp = [ indTmp1 indTmp2 ];
                    if length(indTmp) == 0
                        if strcmpi(g.ignoremissing, 'off')
                            error('Area %s not found ', val);
                        else
                            fprintf('Area %s not found, ignoring it\n', val);
                            indTmp = [];
                        end
                    elseif length(indTmp) > 1
                        if strcmpi(g.ignoremissing, 'off')
                            error('Area %s duplicate', val);
                        else
                            fprintf('Area %s duplicate, ignoring\n', val);
                            indTmp = [];
                        end
                    else
                        inds = [ inds;indTmp ];
                    end
                end
            end
        end
        ROIinds{iCol} = inds;
        EEG.roi.atlas.Scouts(end).Vertices = vertcat(EEG.roi.atlas.Scouts(inds).Vertices);
    end
    
    % add ROIs to connectivity matrix by combining info from origin ROIs (not ideal, better compute it directly)
    if ~isempty(connectmat)
        if iscell(connectmat)
            for iMat = 1:length(connectmat)
                connectmat{iMat} = augmentConnectivity(connectmat{iMat}, ROIinds);
            end
        else
            connectmat = augmentConnectivity(connectmat, ROIinds);
        end
    end
end

% only add ROIs - return
if isempty(roiTable) 
    networks = [];
    return;
end

% define networks
networks = [];
colNames = fieldnames(roiTable);
allLabels = lower({ EEG.roi.atlas.Scouts.Label });
for iCol = 1:size(roiTable,2) % scan columns
    
    networks(end+1).name = colNames{iCol};

    inds = [];
    if isnumeric(roiTable{1,iCol})
        inds = [ roiTable{:,iCol} ]';
    else
        for iRow = 1:size(roiTable,1)
            val = roiTable{iRow, iCol}{1};
            if ~isempty(val)
                indTmp1 = strmatch(lower(val), allLabels, 'exact');
                indTmp2 = strmatch(lower([ 'Brodmann area ' val]), allLabels, 'exact');
                indTmp = [ indTmp1 indTmp2 ];
                if isempty(indTmp)
                    if strcmpi(g.ignoremissing, 'off')
                        error('Area %s not found', val);
                    else
                        fprintf('Area %s not found\n', val);
                    end
                elseif length(indTmp) > 1
                    fprintf('Area %s duplicate, using the first one\n', val);
                    indTmp = indTmp(1);
                end
                inds = [ inds;indTmp ];
            end
        end
    end
    networks(end).ROI_inds = inds';
    
end

EEG.roi.atlas.networks = networks;

% % augment connectivity rows/cols A, B, C, D
% % new areas (A,B) and (C,D)
% % The connectivity between (A,B) and (C,D) is ( A->C + A->D + B->C + B->D )/4 
% % equals (4 -2 +1 -1)/4 = 0.5 in the example below
% connectmat = [ 0 1 4 -2; 1 0 1 -1; 4 1 0 1; -2 -1 1 0];
% ROIinds = { [1 2] [ 3 4] }
% newconnectmat = augmentConnectivity(connectmat,ROIinds)

function newconnectmat = augmentConnectivity(connectmat,ROIinds)

sz = [size(connectmat) 1 1];
permFlag = false;
if sz(2) == sz(3) && sz(1) ~= sz(2)
    connectmat = permute(connectmat, [2 3 1]);
    permFlag = true;
end

% accross subjects if any
nVals = size(connectmat,1);
if size(connectmat,3) > 1
    newconnectmat = zeros(nVals+length(ROIinds), nVals+length(ROIinds), size(connectmat,3));
    for iSubject = 1:size(connectmat,3)
        newconnectmat(:,:,iSubject) = augmentConnectivity(connectmat(:,:,iSubject),ROIinds);
    end
    if permFlag
        newconnectmat = permute(newconnectmat, [3 1 2]);
    end
    return;
end

newconnectmat = zeros(nVals+length(ROIinds), nVals+length(ROIinds));
newconnectmat(1:nVals,1:nVals) = connectmat;

for iCol = 1:length(ROIinds)
    newconnectmat(nVals+iCol,1:nVals) =  mean(connectmat(ROIinds{iCol},:),1);
    newconnectmat(1:nVals,nVals+iCol) =  mean(connectmat(:,ROIinds{iCol}),2);
end
for iCol1 = 1:length(ROIinds)
    for iCol2 = 1:length(ROIinds)
        if iCol1 ~= iCol2
            newconnectmat(nVals+iCol1,nVals+iCol2) =  mean(newconnectmat(nVals+iCol1,ROIinds{iCol2}));
            if 0
                % brute force check (but slower)
                tot = 0;
                for iRoi1 = ROIinds{iCol1}(:)'
                    for iRoi2 = ROIinds{iCol2}(:)'
                        tot = tot + connectmat(iRoi1, iRoi2);
                    end
                end
                tot = tot/length(ROIinds{iCol1})/length(ROIinds{iCol2});
                if abs(tot - newconnectmat(nVals+iCol1,nVals+iCol2)) > 1e-15
                    error('Non equal value');
                end
            end
        end
    end
end

