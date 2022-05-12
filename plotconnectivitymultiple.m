function [imgFileName,txtFileName] = plotconnectivitymultiple(networkfile, connectSpecSelect, varargin)

if nargin < 2
    help plotconnectivitymultiple;
    return
end

g = finputcheck(varargin, { ...
    'subplots'    'string'    {'on' 'off'}    'off';
    'exporttxt'   'string'    {'on' 'off'}    'on';
    'title'       'string'    {}              '';
    'plotmode'    'string'    {'2D' '3D' 'both' }  '2D';
    'filename'    'string'    {}              '';
    }, 'roi_network');
if isstr(g)
    error(g);
end
if ~isempty(g.filename) % remove file extension
    [filePath,g.filename] = fileparts(g.filename);
    g.filename = fullfile(filePath, g.filename);
end
if strcmpi(g.subplots, 'on')
    if ~strcmpi(g.plotmode, '2D')
        error('When using subplots, you cannot use 3-D plots');
    end
end
if ~strcmpi(g.plotmode, '2D')
    if ~exist('roi_plotbrainmovie')
        error('You need to install the Brainmovie plugin to plot sources in 3-D');
    end
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
    tmpTitle = loreta_Networks(iNet).name;
    tmpTitle(tmpTitle == '_') = ' ';
    tmpTitle = [ tmpTitle ' ' g.title ];
    
    % 2-D plot
    if strcmpi(g.plotmode, '2D') || strcmpi(g.plotmode, 'both')
        plotconnectivity(connectSpecSelect{iNet}(:,:), 'labels', labels, 'axis', gca, 'threshold', 0.1);
        h = title(tmpTitle, 'interpreter', 'none');
        pos = get(h, 'position');
        set(h, 'position', pos + [0 0.1 0]);
    
        % save individual plots
        if ~strcmpi(g.subplots, 'on') && ~isempty(g.filename)
            set(h, 'fontsize', 16, 'fontweight', 'bold');
            tmpFileName = [ g.filename '_' loreta_Networks(iNet).name '.jpg' ];
            imgFileName{end+1} = tmpFileName;
            print('-djpeg', tmpFileName );
            close
        end
    end
    
    % 3-D plot
    if strcmpi(g.plotmode, '3D') || strcmpi(g.plotmode, 'both')
        options = {};
        if ~strcmpi(g.subplots, 'on') && ~isempty(g.filename)
            tmpFileName = [ g.filename '_' loreta_Networks(iNet).name '_3d' ];
            options = { 'filename'  tmpFileName };
            imgFileName{end+1} = [ tmpFileName '.xhtml' ];
        end
        roi_plotbrainmovie(connectSpecSelect{iNet}(:,:), 'labels', labels, 'threshold', 0.1, options{:});
    end
    
    % save text
    if strcmpi(g.exporttxt, 'on') && ~isempty(g.filename)
        tmptable = array2table(connectSpecSelect{iNet}(:,:), 'variablenames', labels, 'rownames', labels);
        tmpFileName = [ g.filename '_' loreta_Networks(iNet).name '.txt' ];
        txtFileName{end+1} = tmpFileName;
        writetable(tmptable, tmpFileName,'WriteRowNames', true );
    end
    
end

if strcmpi(g.subplots, 'on')
    h = textsc(titl, 'title');
    set(h, 'fontsize', 16, 'fontweight', 'bold');
    print('-djpeg', g.filename);
    close
end