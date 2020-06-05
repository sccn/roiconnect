function sourcemodel = convertleadfieldtofieldtrip(lf, varargin)

if nargin < 2
    help convertleadfieldtofieldtrip;
    return; 
end

g = finputcheck(varargin, { 'unit'   'string' {} '';
                            'inside' 'string' {} '';
                            'pos'    'real'   {} [];
                            'label'  'cell'   {} {};
                            
                            'unit'   'string' {} 'mm' }, 'convertleadfieldtofieldtrip');
if ischar(g), error(g); end

if isfield(lf, 'Vertices')
    sourcemodel.pos = lf.Vertices;
else
    sourcemodel.pos = g.pos(:,1:3);
end
if isempty(sourcemodel.pos)
    error('Position of vertices must be present either in the Leadfield Matrix structure or provided as additional argument');
end
% leadfield
if isfield(lf, 'LFM')
    sourcemodel.pos = sourcemodel.pos(1:end/3,:); % duplicated positions
    
    sourcemodel.leadfield = cell(size(sourcemodel.pos,1),1);
    nChans = size(lf.LFM,1);
    nChans = 32;
    
    lf.LFM = reshape(lf.LFM, size(lf.LFM,1), size(lf.LFM,2)/3, 3);
    for iPos = 1:length(sourcemodel.leadfield)
        sourcemodel.leadfield{iPos} = squeeze(lf.LFM(1:nChans,iPos,:));
    end
end


% units
if isempty(g.unit)
    if max(sourcemodel.pos(:,1)) > 50
        sourcemodel.unit = 'mm';
    elseif max(sourcemodel.pos(:,1)) > 5
        sourcemodel.unit = 'cm';        
    else
        sourcemodel.unit = 'm';        
    end
    fprintf('Assuming unit in %s based on vertices coordinates\n', sourcemodel.unit);
end

% inside
if isempty(g.inside)
    sourcemodel.inside = ones(size(sourcemodel.pos,1),1) > 0;
    fprintf('Assuming all vertices are contained within the head model\n', sourcemodel.unit);
end

% electrode labels
if ~isempty(g.label)
    sourcemodel.label = g.label;
else
    for iChan = 1:nChans
        sourcemodel.label{iChan} = sprintf('E%d', iChan);
    end
end

% dim ord
sourcemodel.leadfielddimord = '{pos}_chan_ori';
return

% test
cfg = [];
dataPre = eeglab2fieldtrip(EEG, 'preprocessing', 'dipfit');
cfg.method      = 'eloreta';
cfg.sourcemodel = sourcemodel;
%cfg.headmodel   = vol.vol;
source          = ft_sourceanalysis(cfg, dataPre);  % compute the source


