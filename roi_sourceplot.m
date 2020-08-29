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
% Author:  Arnaud Delorme

function roi_sourceplot(freqs, sourceact, sourcemodel, varargin)

if nargin < 3
    help roi_sourceplot;
    return;
end

[g] = finputcheck(varargin, { ...
    'freqselect'      'real'  { }             [] }, 'roi_sourceplot');
if isstr(g)
    error(g);
end

if isempty(g.freqselect)
    indFreq = 1:length(freqs);
elseif length(g.freqselect) == 1
    [~,indFreq] = min(abs(freqs-g.freqselect));
    indFreq = [ indFreq indFreq ];
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

% transform to volume
if ischar(sourcemodel)
    sourceProjtmp = load('-mat', sourcemodel);
end
xl = [min(sourceProjtmp.Vertices(:,1)) max(sourceProjtmp.Vertices(:,1)) ];
yl = [min(sourceProjtmp.Vertices(:,2)) max(sourceProjtmp.Vertices(:,2)) ];
zl = [min(sourceProjtmp.Vertices(:,3)) max(sourceProjtmp.Vertices(:,3)) ];
volMat = zeros(diff(xl)/5+1, diff(yl)/5+1, diff(zl)/5+1);
clear sourcemodelout;
sourcemodelout.dim = size(volMat);
[r,c,v] = ind2sub(size(volMat),find(volMat == 0));
sourcemodelout.pos = [r,c,v];
sourcemodelout.inside    = zeros(size(sourcemodelout.pos,1),1);
sourcemodelout.unit      = 'mm';
%sourcemodelout.transform = traditionaldipfit([xl(1)-5 yl(1)-5 zl(1)-5 0 0 0 5 5 5]);
sourcemodelout.transform = [5 0 0 -75;0 5 0 -105;0 0 5 -45;0 0 0 1]; % 5mm grid
allInds = zeros(size(sourceProjtmp.Vertices));
allIndVolume = zeros(length(sourceProjtmp.Vertices),1);
for iVert = 1:length(sourceProjtmp.Vertices)
    xVal = (sourceProjtmp.Vertices(iVert,1)-xl(1))/5+1;
    yVal = (sourceProjtmp.Vertices(iVert,2)-yl(1))/5+1;
    zVal = (sourceProjtmp.Vertices(iVert,3)-zl(1))/5+1;
    ind = sub2ind(size(volMat), xVal, yVal, zVal);
    volMat(xVal, yVal, zVal) = mean(sourceact(iVert,indFreq), 2);
    allIndVolume(iVert) = ind;
    allInds(iVert,:) = [xVal yVal zVal];
    sourcemodelout.inside(ind) = true;
end

res = squeeze(volMat(:,:,10));
res = res';
res(:,:) = res(end:-1:1,:);
cmap = colormap('jet');
mi = min(res(:));
mx = max(res(:));
resrgb = ones([size(res) 3]);
for iPix1 = 1:size(res,1)
    for iPix2 = 1:size(res,2)
        if res(iPix1,iPix2) ~= 0
            ind = ceil((res(iPix1,iPix2)-mi)/(mx-mi)*(size(cmap,1)-1))+1;
            resrgb(iPix1,iPix2,:) = cmap(ind,:);
        end
    end
end

imagesc(resrgb); axis equal; axis off;
%figure; imagesc(squeeze(volMat(:,:,10))); axis equal; axis off;
return

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
% sourceProj2 = sourceProj;
% sourceProj2.time = freqs;
% for iVert = 1:length(allIndVolume)
%     sourceProj2.mom{allIndVolume(iVert)} = brainDx(iVert,:);
% end

cfg              = [];
cfg.method       = 'otho';
cfg.funparameter = 'mom';
cfg.location     = [26 8 10];
cfg.latency      = 10;
cfg.slicepos     = [8 9];
ft_sourceplot(cfg, sourceProj);