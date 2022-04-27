function [imgFileName,txtFileName] = plotconnectivitymultiple(networkfile, connectSpecSelect, varargin)

if nargin < 2
    help plotconnectivitymultiple;
    return
end

g = finputcheck(varargin, { ...
    'subplots'    'string'    {'on' 'off'}    'off';
    'exporttxt'   'string'    {'on' 'off'}    'on';
    'title'       'string'    {}              '';
    'filename'    'string'    {}              '';
    }, 'roi_network');
if isstr(g)
    error(g);
end
if ~isempty(g.filename) % remove file extension
    [filePath,g.filename] = fileparts(g.filename);
    g.filename = fullfile(filePath, g.filename);
end

if nargin < 5
    % create a single packet
    singlefig = true;
end

% decode network param
if ischar(networkfile)
    networkfile = load('-mat', networkfile);
end
loreta_Networks = networkfile.loreta_Networks;
loreta_ROIS     = networkfile.loreta_ROIS;

if strcmpi(g.subplots, 'on')
    figure('position', [100 100 1000 700], 'paperpositionmode', 'auto');
    ncol = ceil(sqrt(length(loreta_Networks)));
    nrow = ceil(length(loreta_Networks)/ncol);
end
imgFileName = {};
txtFileName = {};
for iNet = 1:length(loreta_Networks)
    if strcmpi(g.subplots, 'on')
        subplot(nrow, ncol, iNet)
    else
        figure('position', [100 100 400 700], 'paperpositionmode', 'auto');
    end
    
    labels = { loreta_ROIS(loreta_Networks(iNet).ROI_inds).Label };
    plotconnectivity(connectSpecSelect{iNet}(:,:), 'labels', labels, 'axis', gca, 'threshold', 0.1);
    tmpTitle = loreta_Networks(iNet).name;
    tmpTitle(tmpTitle == '_') = ' ';
    tmpTitle = [ tmpTitle ' ' g.title ];
    h = title(tmpTitle, 'interpreter', 'none');
    pos = get(h, 'position');
    set(h, 'position', pos + [0 0.1 0]);
    
    % save individual plots
    if strcmpi(g.exporttxt, 'on') && ~isempty(g.filename)
        tmptable = array2table(connectSpecSelect{iNet}(:,:), 'variablenames', labels, 'rownames', labels);
        tmpFileName = [ g.filename '_' loreta_Networks(iNet).name '.txt' ];
        txtFileName{end+1} = tmpFileName;
        writetable(tmptable, tmpFileName,'WriteRowNames', true );
    end
    
    % save text
    if ~strcmpi(g.subplots, 'on') && ~isempty(g.filename)
        set(h, 'fontsize', 16, 'fontweight', 'bold');
        tmpFileName = [ g.filename '_' loreta_Networks(iNet).name '.jpg' ];
        imgFileName{end+1} = tmpFileName;
        print('-djpeg', tmpFileName );
        close
    end
end

if strcmpi(g.subplots, 'on')
    h = textsc(titl, 'title');
    set(h, 'fontsize', 16, 'fontweight', 'bold');
    print('-djpeg', g.filename);
    close
end