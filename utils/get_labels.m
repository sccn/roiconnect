function labels = get_labels(EEG)
    % retrieve labels from atlas
    labels = strings(1,length(EEG.roi.atlas.Scouts));
    for i = 1:length(labels)
        scout = struct2cell(EEG.roi.atlas.Scouts(i));
        labels(i) = char(scout(1));
    end
    labels = cellstr(labels);

    % remove region labels that were not selected
    if isfield(EEG.roi, 'roi_selection')
        if ~isempty(EEG.roi.roi_selection)
            labels = labels(cell2mat(EEG.roi.roi_selection));
        end
    end
end
