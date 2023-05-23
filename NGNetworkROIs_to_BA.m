clear
allAreas = readtable('NGNetworkROIs_area_definition_v2.txt', 'delimiter', char(9));
networks = readtable('NGNetworkROIs_v4.txt', 'delimiter', char(9));
networksNew = networks;
networksNew(2:end,:) = [];

networkNames  = fieldnames(networks);
allAreasNames = fieldnames(allAreas);
allAreasNames = allAreasNames(1:end-3);
for iNet = 1:length(networkNames)-3

    count = 1;
    areas = networks.(networkNames{iNet});
    addAreaList = {};

    % get the list of areas
    for iArea = 1:length(areas)
        areaNameTmp = areas{iArea};

        if ~isempty(areaNameTmp)
            if areaNameTmp(end) ~= 'L' && areaNameTmp(end) ~= 'R'
                areaNameTmp1 = [ areaNameTmp '_R' ];
                areaNameTmp2 = [ areaNameTmp '_L' ];
                colPos1 = strmatch(areaNameTmp1, allAreasNames, 'exact');
                colPos2 = strmatch(areaNameTmp2, allAreasNames, 'exact');
                addAreaList = { addAreaList{:} allAreas{:,colPos1}{:} allAreas{:,colPos2}{:} };
            else
                colPos = strmatch(areaNameTmp, allAreasNames, 'exact');
                addAreaList = { addAreaList{:} allAreas{:,colPos}{:} };
            end
        end
    end

    % write the list of areas
    addAreaList(cellfun(@isempty, addAreaList)) = [];
    for iArea = 1:length(addAreaList)
        networksNew(iArea,iNet) = { addAreaList{iArea} };
    end
end
writetable(networksNew, 'BANetworkROIs_v4.txt', 'delimiter', char(9));
