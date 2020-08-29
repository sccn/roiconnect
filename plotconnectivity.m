% plotconnectivity() - plot coher map

function plotconnectivity(array, varargin)

radius = 0.5;
linewidth = 1;

g = finputcheck(varargin, { ...
    'labels'      'cell'      { }             {};
    'axis'        ''          {}              [];
    'threshold'   'real'      {}              0.25;
    }, 'roi_network');
if isstr(g)
    error(g);
end

if g.threshold > 0
    array(array < g.threshold) = 0;
end

if size(array,1) ~= size(array,2)
    error('Input array must be square');
end

anglesInit = linspace(0,2*pi,size(array,1)+1) + pi/size(array,1);
x = sin(anglesInit)*radius;
y = cos(anglesInit)*radius;

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

% find channel coordinates
% ------------------------
if isempty(g.labels)
    for iPnt = 1:length(anglesInit)
        g.labels{iPnt} = sprintf('  Area %d', iPnt);
    end
elseif length(anglesInit)-1 ~= length(g.labels)
    error('Wrong number of labels');
else
    for iPnt = 1:length(g.labels)
        g.labels{iPnt} = strrep(g.labels{iPnt}, 'Brodmann area', 'BA');
        g.labels{iPnt} = [ ' ' g.labels{iPnt} ];
    end
end

% make lines between pairs of electrodes
% --------------------------------------
if isempty(g.axis)
    figure; hold on;
else
    axes(g.axis); hold on;
end
axis equal;
axis off;

plot(x,y,'k');
plot(x,y,'.','markersize', 12);
x(end) = [];
y(end) = [];

warning off;
color = 'r';
for ind1 = 1:size(array,1)
    for ind2 = 1:size(array,2)
        if ind1 ~= ind2
            if array(ind1, ind2) > 0

                aa = [x(ind1) y(ind1)];
                bb = [x(ind2) y(ind2)];
                distance = sqrt(sum(abs(aa-bb).^2));

                center = distance*4*2*(aa+bb)/2;
                radius = sqrt(sum(abs(aa-center).^2));
                if sum(abs(center)) < 1e-8
                    plot([aa(1) bb(1)],[aa(2) bb(2)],'-');
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
                    plot(x2,y2,'-');
                    0;
                end
            end
        end
    end
    h = text( x(ind1), y(ind1), 0, g.labels{ind1});
    set(h, 'HorizontalAlignment','left', 'rotation', 90-anglesInit(ind1)/pi*180);
    0;
end
xlim([-0.7 0.7]);
ylim([-0.7 0.7]);

