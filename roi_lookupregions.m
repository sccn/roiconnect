% roi_lookupregions() - look up regions from network name
%
% Usage:
%     regions = roi_lookupregions(networkname, EEG);
%     regions = roi_lookupregions(networkname, networkdefs);
%
% Input:
%     networkname - [string] name of network.
%     EEG         - [EEGLAB structure] use EEG dataset to look up network
% Or
%     networkdefs - [string|struct] network definition file (same input as
%                   roi_definenetwork()). Alternatively network structure
%                   (second output of roi_definenetwork())
%                   You may use also EEG.roi.atlas.networks
%
% Outputs:
%     regions - [cell array of string] name of the regions making up the
%               network
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

function regions = roi_lookupregions(networkname, networkdefs)

regions = {};
if ischar(networkdefs)
    if ~exist(networkdefs, 'file')
        p = fileparts(which('roi_definenetwork'));
        roiTable2 = fullfile(p, networkdefs);
        if ~exist(roiTable2, 'file')
            error('File not found %s', networkdefs);
        end
        roiTable = roiTable2;
    end
    roiTable = readtable(networkdefs,'Delimiter', char(9));
    allnetworknames = fieldnames(roiTable);
elseif ~isstruct(networkdefs)
    roiTable = networkdefs;
    allnetworknames = fieldnames(roiTable);
else
    EEG = networkdefs;
    if isfield(EEG, 'roi') && isfield(EEG.roi, 'atlas')
        atlas = EEG.roi.atlas;
    else
        atlas = EEG;
    end
    if ~isfield(atlas, 'networks')
        error('Input structure must contain a field named "networks"')
    end

    allnetworknames = { atlas.networks.name };
end

indNet = [];
len    = 0;
for iNet = 1:length(allnetworknames)
    if contains(networkname, allnetworknames{iNet})
        if length(  allnetworknames{iNet} ) > len
            indNet = iNet;
            len    = length( allnetworknames{iNet} );
        end
    end
end

if isempty(indNet)
    fprintf('Network %s not found\n', networkname);
elseif ~isstruct(networkdefs)
    regions = roiTable.(allnetworknames{iNet});
else
    regions = { atlas.Scouts(atlas.networks(indNet).ROI_inds).Label };
end
