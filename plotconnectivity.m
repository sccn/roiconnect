% plotconnectivity() - plot circle or brain image showing connectivity 
%                      between regions
%
% Usage:
%     plotconnectivity(array, 'key', val);
%
% Input:
%     array   - [n x n] square array indicating which cell are connected
%
% Optional inputs:
%    'labels'  - [cell] name for each row/column
%    'axis'    - [axis handle] axis to plot the figure (otherwise creates
%                a new figure)
%    'brainimg' - ['on'|'off'] plot results on a 2-D image of the brain
%    'threshold' - [real] only show connections above a given threshold
%
% Author: Arnaud Delorme

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

% Test
% plotconnectivity(rand(4,4), 'labels', { 'Dorso_lateral_prefrontal_cortex' 'Parietal_lobe' 'Thalamus' 'Visual_cortex' }, 'brainimg', 'off');
% plotconnectivity(rand(8,8), 'brainimg', 'bilateral', 'labels', { 'Dorso_lateral_prefrontal_cortex_L' 'Parietal_lobe_L' 'Thalamus_L' 'Visual_cortex_L' 'Dorso_lateral_prefrontal_cortex_R' 'Parietal_lobe_R' 'Thalamus_R' 'Visual_cortex_R' });
% plotconnectivity(rand(4,4), 'brainimg', 'bilateral', 'labels', { 'Dorso_lateral_prefrontal_cortex_L' 'Parietal_lobe_L' 'Thalamus_L' 'Visual_cortex_L' });

function limits = plotconnectivity(array, varargin)

if nargin < 2
    help plotconnectivity
    return
end

radius = 0.5;
linewidth = 1;

g = finputcheck(varargin, { ...
    'labels'      'cell'      { }             {};
    'labelsgroup'  'cell'      { }             {};
    'axis'        ''          {}              [];
    'colormap'    ''          {}              cool;
    'limits'      'real'      {}              [];
    'cbar'        'string'     {'on' 'off' }     'off';
    'brainimg'   'string'     {'on' 'off' 'bilateral'}     'bilateral';
    'threshold'   'real'      {}              0.25;
    }, 'roi_network');
if isstr(g)
    error(g);
end
limits = g.limits;

if g.threshold > 0
    array(abs(array) < g.threshold) = 0;
end

if size(array,1) ~= size(array,2)
    error('Input array must be square');
end
if isempty(g.labels)
    if strcmpi(g.brainimg, 'on')
        disp('Cannot plot on brain with area labels')
        g.brainimg = 'off';
    end
end

% colors for labelgroups
if ~isempty(g.labelsgroup)
    groups = unique(g.labelsgroup);
    %               green         yellow              purple    grey
    colors = { 'r' [0 0.8 0] 'b' [0.7 0.7 0] 'k' 'm' [0.6 0.6 1] [0.5 0.5 0.5]};

    % reorder by group
    allinds = [];
    for iGroup = 1:length(groups)
        inds = strmatch(groups{iGroup}, g.labelsgroup, 'exact');
        inds = sort(inds);
        allinds = [ allinds; inds];
    end
    g.labelsgroup = g.labelsgroup(allinds);
    g.labels      = g.labels(allinds);
    array = array(allinds, :);
    array = array(:,allinds);
end

if ~strcmpi(g.brainimg, 'off')
    p = fileparts(which('plotconnectivity.m'));
    if strcmpi(g.brainimg, 'bilateral')
        [img, map, alphachannel]  = imread(fullfile(p, 'brain2.png'));
        coords = loadtxt('brain_coords2.txt');
    else
        [img, map, alphachannel]  = imread(fullfile(p, 'brain.png'));
        coords = loadtxt('brain_coords.txt');
    end
    coords(:,1) = [];
    for indLab = 1:length(g.labels)
        indCoord = strmatch( lower(g.labels{indLab}), lower(coords(1,:)) );
        if isempty(indCoord)
            disp('Could not find brain areas, plotting on circle')
            g.brainimg = 'off';
            break;
            indLab = length(g.labels);
        else
            x(indLab) = coords{2,indCoord};
            y(indLab) = coords{3,indCoord};
        end
    end
end
if strcmpi(g.brainimg, 'off')
    plotImg = false;
    anglesInit = linspace(0,2*pi,size(array,1)+1) + pi/size(array,1);
    x = sin(anglesInit)*radius;
    y = cos(anglesInit)*radius;
end

% settings
% --------
% g = struct(varargin{:});
% try g.maxcoh;     catch, g.maxcoh   = max(max(abs(array))); end
% try g.colormap;   catch, g.colormap = 'bluered'; end
% try g.electrodes; catch, g.electrodes = 'off'; end
% try g.usethreshold; catch, g.usethreshold = 2; end
% if strcmpi(g.colormap, 'bluered'), cmap = redbluecmap; cmap = cmap(end:-1:1,:);
% else                               cmap = yellowredbluecmap;
% end

% make lines between pairs of electrodes
% --------------------------------------
if isempty(g.axis)
    if ~strcmpi(g.brainimg, 'off')
        figure('position', [0 0 400 700])
    else
        figure;
    end
else
    axes(g.axis); hold on;
end
if ~strcmpi(g.brainimg, 'off')
    imagesc(img, 'AlphaData', alphachannel); axis off;
    alpha(0.2)
    hold on;
    axis equal
    set(gca, 'ydir', 'reverse');
    pos = get(gca, 'position');
    axes('position', pos); axis off; hold on;
    set(gca, 'ydir', 'normal');
else
    hold on;
    axis equal
end
g.axis = gca;
axis equal;
axis off;

% plot dots
if ~strcmpi(g.brainimg, 'off')
    plot(x,y,'k.');
    plot(x,y,'.','markersize', 12);
else
    plot(x,y,'k-');
    x(end) = [];
    y(end) = []; % remove duplicate last point
    for iX = 1:length(x)
        if ~isempty(g.labelsgroup)
            ind = strmatch(g.labelsgroup{iX}, groups, 'exact');
            col = colors{ind}; %mod(ind-1, length(colors))+1};
        else
            col = 'r';
        end
        plot(x(iX),y(iX),'.','markersize', 12, 'color', col);
    end
end

% rename labels
% -------------
if ~strcmpi(g.brainimg, 'off') % CTC added .... to branch around code if not plotting
if isempty(g.labels)
    for iPnt = 1:length(anglesInit)
        g.labels{iPnt} = sprintf('  Area %d', iPnt);
    end
elseif length(x) ~= length(g.labels)
    error('Wrong number of labels');
else
    for iPnt = 1:length(g.labels)
        g.labels{iPnt} = strrep(g.labels{iPnt}, 'Brodmann area', 'BA');
        g.labels{iPnt} = [ ' ' g.labels{iPnt} ];
    end
end
end % CTC added .... to branch around above code if not plotting
warning off;
if isempty(limits)
    limits(2) = max(array(:));
    arrayTmp = array(:);
    arrayTmp(arrayTmp == 0) = Inf;
    limits(1) = min(arrayTmp(:));
end
for ind1 = 1:size(array,1)
    for ind2 = 1:size(array,2)
        if ind1 ~= ind2
            if array(ind1, ind2) ~= 0

                aa = [x(ind1) y(ind1)];
                bb = [x(ind2) y(ind2)];
                distance = sqrt(sum(abs(aa-bb).^2));

                center = distance*4*2*(aa+bb)/2;
                radius = sqrt(sum(abs(aa-center).^2));
                value0to1 = (array(ind1, ind2)-limits(1))/(limits(2)-limits(1));
                col = ceil( value0to1*( size(g.colormap,1)-1 ) )+1;
                col = max(col,1);
                col = min(col,size(g.colormap,1));
                if sum(abs(center)) < 1e-8 || ~strcmpi(g.brainimg, 'off')
                    patchline([aa(1) bb(1)],[aa(2) bb(2)],'edgecolor', g.colormap(col, :), 'edgealpha', max(0.1,col/256), 'linewidth', 2);
                    %plot([aa(1) bb(1)],[aa(2) bb(2)],'-', 'color', g.colormap(col, :), 'alpha', 0.5, 'linewidth', 2);
                else
                    angle1 = atan2(aa(1)-center(1), aa(2)-center(2));
                    angle2 = atan2(bb(1)-center(1), bb(2)-center(2));
                    angles = [ angle1 angle2 ];
                    if angle2 < angle1
                        angles = [ angle2 angle1 ];
                    end
                    if diff(angles) > pi
                        angles = [angles(2) 2*pi+angles(1)];
                    end

                    pnts = linspace(angles(1),angles(2),round(diff(angles)*10));
                    x2 = sin(pnts)*radius+center(1);
                    y2 = cos(pnts)*radius+center(2);
                    plot(x2,y2,'-', 'color', g.colormap(col, :));
                end
            end
        end
    end
    if ~strcmpi(g.brainimg, 'off')
        str = formatlabel(g.labels{ind1});
        xx = x(ind1)-0.1;
        if any(str == 10)   
            yy = y(ind1)+0.1;
        else
            yy = y(ind1)+0.05;
        end
        h = text( xx, yy, 0, str, 'interpreter', 'none', 'fontsize', 8);
    else
        if ~isempty(g.labelsgroup)
            ind = strmatch(g.labelsgroup{ind1}, groups, 'exact');
            col = colors{ind}; %mod(ind-1, length(colors))+1};
        else
            col = 'k';
        end
        h = text( x(ind1), y(ind1), 0, [ ' ' g.labels{ind1} ], 'interpreter', 'none', 'fontsize', 8, 'color', col);
        set(h, 'HorizontalAlignment','left', 'rotation', 90-anglesInit(ind1)/pi*180);
    end
end
if ~strcmpi(g.brainimg, 'off')
    if strcmpi(g.brainimg, 'bilateral')
        xlim([0 1])
        ylim([0 2])
    else
        xlim([0 1])
        ylim([0 1])
    end
else
    xlim([-0.7 0.7]);
    ylim([-0.7 0.7]);
end
if strcmpi(g.cbar, 'on')
    h = cbar;
    n = size(g.colormap,1);
    transparency = linspace(0,1,n)';
    transparency(transparency < 0.1) = 0.1;
    imagesc([0 1], linspace(limits(1), limits(2), n), reshape(cool,n,1,3), 'AlphaData', transparency);
    set(h, 'xticklabel', [], 'YAxisLocation', 'right')
    ylim(limits);
    axes(g.axis);
end

function str = formatlabel(str)

str = strip(str);
str(str == '_') = ' ';
if length(str) > 10
    sp = find(str == ' ');
    if length(sp) > 1 && length(str(sp(2)+1:end) ) > 1
        str = [ str(1:sp(2)-1) 10 str(sp(2)+1:end) ];
    end
end