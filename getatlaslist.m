function [atlasstr, atlas] = getatlaslist(fileName)

data = load('-mat', fileName);

atlasstr = {};
atlas    = {};
for iAtlas = 1:length(data.Atlas)
    if ~isempty(data.Atlas(iAtlas).Scouts) && ~strcmpi(data.Atlas(iAtlas).Name, 'Structures')
        atlasstr{end+1} = [data.Atlas(iAtlas).Name ' (' int2str(length(data.Atlas(iAtlas).Scouts)) ' ROIs)' ];
        atlas{   end+1} = data.Atlas(iAtlas).Name;
    end
end