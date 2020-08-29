% surface only
function [sourcemodelout, transform] = transform_move_inward(sourcemodel, headmodel, transform)

if ischar(headmodel)
    headmodel = load('-mat', headmodel);
    if isfield(headmodel, 'vol')
        headmodel = headmodel.vol;
        headmodel.unit = 'mm';
    end
end
if ischar(sourcemodel)
    try
        sourcemodel = load('-mat', sourcemodel);
    catch
        % Likely a volume atlas
        sourcemodelout = sourcemodel;
        return
    end
    if isfield(sourcemodel, 'cortex')
        sourcemodel = sourcemodel.cortex;
    end
end
if isfield(sourcemodel, 'inside')
    pos = sourcemodel.transform * [sourcemodel.pos(logical(sourcemodel.inside),:) ones(sum(sourcemodel.inside),1) ]';
    sourcemodel = [];
    sourcemodel.pos = pos(1:3,:)';
end
    
newsourcemodel = [];
if isfield(sourcemodel, 'Vertices') && isfield(sourcemodel, 'Faces')
    newsourcemodel.pos = sourcemodel.Vertices;
    newsourcemodel.tri = sourcemodel.Faces;
elseif isfield(sourcemodel, 'Vertices')
    newsourcemodel.pos = sourcemodel.Vertices;
    newsourcemodel.tri = [];
elseif isfield(sourcemodel, 'vertices')
    newsourcemodel.pos = sourcemodel.vertices;
    newsourcemodel.tri = sourcemodel.faces;
else
    newsourcemodel.pos = sourcemodel.pos;
    if isfield(newsourcemodel, 'tri')
        newsourcemodel.tri = sourcemodel.tri;
    else
        newsourcemodel.tri = [];
    end
end

cfg = [];
pos = [newsourcemodel.pos ones(size(newsourcemodel.pos,1),1) ];
if ~isempty(transform)
    pos = traditionaldipfit(transform)*pos';
else
    pos = pos';
end
pos(4,:) = [];
cfg.sourcemodel.pos = pos';
cfg.sourcemodel.tri = newsourcemodel.tri;
cfg.sourcemodel.unit = headmodel.unit;
cfg.moveinward = 1;
cfg.headmodel = headmodel;
disp('moving source model inward if necessary');
sourcemodelout = ft_prepare_sourcemodel(cfg);
transform = [];

% plot3dmeshalign(headmodel);
% 
% hold on;
% plot3dmeshalign(tmp2, [], [1 0 0])
