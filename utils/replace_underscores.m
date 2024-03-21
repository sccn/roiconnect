function new_labels = replace_underscores(labels)
    % remove underscores in label names to avoid bug
    new_labels = strrep(labels, '_', ' ');
end