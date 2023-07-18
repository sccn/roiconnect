% roi_sourceplot - plot activity in source model (even when it contains
%                  holes)
%
% Usage:
%  EEG = roi_sourceplot(freqs, sourceact, sourcemodel, 'key', val);
%
% Inputs:
%  freqs       - [real] array of frequencies
%  sourceact   - [voxels x freq] source activities
%  sourcemodel - [string] file name of source model or source model
%                structure. Must contain a field Vertices of [voxels x 3]
%                MNI locations.
%
% Required inputs:
%  'freqselect' - [real] frequency of interest or frequency range of interest.
%                 Defaut is all frequencies.
%
% Example:
%  % run ROI activity using ROI connect to compute activity in ROI (and voxels)
%  roi_sourceplot(EEG.roi.freqs, EEG.roi.source_voxel_power, EEG.roi.cortex );
%
% Author:  Arnaud Delorme

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

function alldata = roi_sourceplot(freqs, sourceact, sourcemodel, varargin)

if nargin < 3
    help roi_sourceplot;
    return;
end

g = finputcheck(varargin, { ...
    'freqrange'       {'cell' 'real'}  { [] [] }       {};
    'saveasfile'      'string'         { }             '';
    'mask'            'string'         { }             '';
    'noplot'          'string'         { 'on' 'off' }  'off';
    'precomputed'     'struct'         { }             struct([]);
    'limits'          ''               []              [];
    'slice'           'integer'        []              []}, 'roi_sourceplot');
if ischar(g)
    error(g);
end
if isempty(g.freqrange)
    g.freqrange = [freqs(1) freqs(end)];
end
if ~iscell(g.freqrange)
    g.freqrange = { g.freqrange };
end

if strcmpi(g.noplot, 'off')
    figure('paperpositionmode', 'auto', 'position', [1440  200  814 1138]);
end
alldata = [];
for iFreq = 1:length(g.freqrange)
    
   
    % transform to volume
    if ischar(sourcemodel)
        sourceProjtmp = load('-mat', sourcemodel);
    else
        sourceProjtmp = sourcemodel;
    end
    if ~isfield(sourceProjtmp, 'Vertices')
        sourceProjtmp.Vertices = sourceProjtmp.pos;
    end

    stride = unique(abs(diff(sourceProjtmp.Vertices)));
    if iscell(sourceact)
        % select precomputed activity matrix
        volMat = loreta2volume(sourceact{iFreq}, sourceProjtmp.Vertices, [stride(2) stride(2) stride(2)]);
    else
        % select frequency range in activity matrix
        g.freqselect = g.freqrange{iFreq};
        if isempty(g.freqselect)
            indFreq = 1:length(freqs);
        elseif length(g.freqselect) == 1
            [~,indFreq] = min(abs(freqs-g.freqselect));
            indFreq(2) = indFreq;
        elseif length(g.freqselect) == 2
            [~,indFreq1] = min(abs(freqs-g.freqselect(1)));
            [~,indFreq2] = min(abs(freqs-g.freqselect(2)));
            indFreq = indFreq1:indFreq2;
        else
            error('Frequency selection must be an array of 1 or 2 elements');
        end
        if isempty(indFreq)
            error('No frequency found');
        end
        volMat = loreta2volume(mean(sourceact(:,indFreq), 2), sourceProjtmp.Vertices, [stride(2) stride(2) stride(2)]);
    end
    if stride(2) ~= floor(stride(2))
        error('Non-integer stride')
    end

%     xl = [min(sourceProjtmp.Vertices(:,1)) max(sourceProjtmp.Vertices(:,1)) ];
%     yl = [min(sourceProjtmp.Vertices(:,2)) max(sourceProjtmp.Vertices(:,2)) ];
%     zl = [min(sourceProjtmp.Vertices(:,3)) max(sourceProjtmp.Vertices(:,3)) ];
%     downscale = 10;
%     volMat = zeros(ceil(diff(xl)/downscale+1), ceil(diff(yl)/downscale+1), ceil(diff(zl)/downscale+1));
%     clear sourcemodelout;
%     sourcemodelout.dim = size(volMat);
%     [r,c,v] = ind2sub(size(volMat),find(volMat == 0));
%     sourcemodelout.pos = [r,c,v];
%     sourcemodelout.inside    = zeros(size(sourcemodelout.pos,1),1);
%     sourcemodelout.unit      = 'mm';
%     %sourcemodelout.transform = traditionaldipfit([xl(1)-5 yl(1)-5 zl(1)-5 0 0 0 5 5 5]);
%     sourcemodelout.transform = [5 0 0 -75;0 5 0 -105;0 0 5 -45;0 0 0 1]; % 5mm grid
%     allInds = zeros(size(sourceProjtmp.Vertices));
%     allIndVolume = zeros(length(sourceProjtmp.Vertices),1);
%     for iVert = 1:length(sourceProjtmp.Vertices)
%         xVal = round((sourceProjtmp.Vertices(iVert,1)-xl(1))/downscale)+1;
%         yVal = round((sourceProjtmp.Vertices(iVert,2)-yl(1))/downscale)+1;
%         zVal = round((sourceProjtmp.Vertices(iVert,3)-zl(1))/downscale)+1;
%         ind = sub2ind(size(volMat), xVal, yVal, zVal);
%         volMat(xVal, yVal, zVal) = mean(sourceact(iVert,indFreq), 2);
%         allIndVolume(iVert) = ind;
%         allInds(iVert,:) = [xVal yVal zVal];
%         sourcemodelout.inside(ind) = true;
%     end
    
    % put precomputed data in VolMat
    if ~isempty(g.slice)
        tmpSlice = g.slice;
    else
        tmpSlice = ceil(linspace(3, size(volMat,3)-2, 5));
    end
    for iSlice = 1:length(tmpSlice)
        fieldVal = sprintf('loreta%1.0fto%1.0fHz_slice%d', g.freqrange{iFreq}(1), g.freqrange{iFreq}(2), tmpSlice(iSlice));
        if isfield(g.precomputed, fieldVal)
            res = g.precomputed.(fieldVal);
            res(isnan(res)) = 0;
            res(:,:) = res(end:-1:1,:);
            res = res';
            volMat(:,:,tmpSlice(iSlice)) = res;
        end
    end
    
    % plot
    res = squeeze(volMat(:,:,tmpSlice));
    if isempty(g.limits)
        mi = min(res(:),[],'omitnan');
        mx = max(res(:),[],'omitnan');
        fprintf('Loreta limits from data: %1.2f to %1.2f\n', mi, mx);
    else
        mi = 0;
        mx = g.limits(iFreq);
    end
    cmap = colormap('turbo');
    for iSlice = 1:length(tmpSlice)
        res = squeeze(volMat(:,:,tmpSlice(iSlice)));
        res = res';
        res(:,:) = res(end:-1:1,:);
        
        % save and retreive data
        fieldVal = sprintf('loreta%1.0fto%1.0fHz_slice%d', g.freqrange{iFreq}(1), g.freqrange{iFreq}(2), tmpSlice(iSlice));
        alldata.(fieldVal) = res;
        
        if strcmpi(g.noplot, 'off')
            resrgb = ones([size(res) 3]);
            for iPix1 = 1:size(res,1)
                for iPix2 = 1:size(res,2)
                    if res(iPix1,iPix2) ~= 0
                        if isnan(res(iPix1,iPix2))
                            resrgb(iPix1,iPix2,:) = [0.9 0.9 0.9];
                        else
                            ind = ceil((res(iPix1,iPix2)-mi)/(mx-mi)*(size(cmap,1)-1))+1;
                            ind = max(1, ind);
                            ind = min(size(cmap,1), ind);
                            resrgb(iPix1,iPix2,:) = cmap(ind,:);
                        end
                    end
                end
            end
            
            subplot(length(tmpSlice), length(g.freqrange), iFreq + length(g.freqrange)*(iSlice-1));
            imagesc(resrgb); axis equal; axis off;
            if iSlice == 1
                if length(g.freqrange{iFreq}) == 2
                    h = title(sprintf('%1.1f-%1.1f Hz', g.freqrange{iFreq}(1), g.freqrange{iFreq}(2)));
                else
                    h = title(sprintf('%1.1f Hz', g.freqrange{iFreq}));
                end
                set(h, 'fontsize', 12);
            end
        end
    end
end
if ~isempty(g.saveasfile) && strcmpi(g.noplot, 'off')
    print('-djpeg', g.saveasfile);
    close
end
return

%figure; imagesc(squeeze(volMat(:,:,10))); axis equal; axis off;

% mimick a source to plot
% sourceProj =
%
%   struct with fields:
%
%        dim: [29 35 22]
%        pos: [22330×3 double]
%       time: [1×129 double]
%        mom: {22330×1 cell}
%     inside: [22330×1 logical]
%        cfg: [1×1 struct]
%
sourceProj = [];
sourceProj.cfg = [];
sourceProj.time   = freqs;
sourceProj.dim    = sourcemodelout.dim;
sourceProj.pos    = sourcemodelout.pos;
sourceProj.inside = sourcemodelout.inside;
sourceProj.mom    = cell(length(sourcemodelout.inside),1);
insideInds = find(sourcemodelout.inside);
for iVert = 1:length(allIndVolume)
    sourceProj.mom{allIndVolume(iVert)} = sourceact(iVert,:);
end

cfg              = [];
cfg.method       = 'otho';
cfg.funparameter = 'mom';
cfg.location     = [26 8 10];
cfg.latency      = 10;
cfg.slicepos     = [8 9];
ft_sourceplot(cfg, sourceProj);

% convert coordinates to volume
% -----------------------------
function vol = loreta2volume(act, pos, stride)

minx = min(pos(:,1)); maxx = max(pos(:,1));
miny = min(pos(:,2)); maxy = max(pos(:,2));
minz = min(pos(:,3)); maxz = max(pos(:,3));

vol = zeros((maxx-minx)/stride(1), (maxy-miny)/stride(2), (maxz-minz)/stride(3));

for iPos = 1:size(pos,1)
    vol( (pos(iPos,1)-minx)/stride(1)+1, ...
         (pos(iPos,2)-miny)/stride(2)+1, ...
         (pos(iPos,3)-minz)/stride(3)+1) = act(iPos);
    % in case using the sourcemodel with shrunk coordinates
    %     vol( floor((pos(iPos,1)-minx)/stride(1)+1), ...
    %          floor((pos(iPos,2)-miny)/stride(2)+1), ...
    %          floor((pos(iPos,3)-minz)/stride(3)+1)) = act(iPos);
end
