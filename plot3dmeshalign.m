function plot3dmeshalign(filename, transform, color1)

if isempty(filename)
    return;
end

if ischar(filename)
    [~,~,ext] = fileparts(filename);
    if strcmpi(ext, '.nii')
        atlas = ft_read_atlas(filename);
        mri = sum(atlas.tissue(:,:,:,:),4) > 0;
        [r,c,v] = ind2sub(size(mri),find(mri));
        xyz = [r c v ones(length(r),1)];
        xyz = atlas.transform*xyz';
        if nargin > 1 && ~isempty(transform)
            xyz = traditionaldipfit(transform)*xyz;
        end
        plot3(xyz(1,:),xyz(2,:),xyz(3,:), '.');
        return
    elseif strcmpi(ext, '.head')
        afni = ft_read_atlas(filename);
        mri = sum(afni.brick0(:,:,:,:),4) > 0;
        mri = afni.brick1(:,:,:,1);
        [r,c,v] = ind2sub(size(mri),find(mri));
        xyz = [r c v ones(length(r),1)];
        xyz = afni.transform*xyz';
        if nargin > 1 && ~isempty(transform)
            xyz = traditionaldipfit(transform)*xyz;
        end
        plot3(xyz(1,:),xyz(2,:),xyz(3,:), '.');
        return
    elseif strcmpi(ext, '.mat')
        f = load('-mat', filename);
    else
        error('Unknown file format');
    end
else
    f = filename;
end
if isfield(f, 'cortex')
    f = f.cortex;
end
if isfield(f, 'SurfaceFile') % Brainstrom leadfield
    p = fileparts(fileparts(fileparts(filename)));
    try
        f = load('-mat', fullfile(p, 'anat', f.SurfaceFile));
    catch
        error('Cannot find Brainstorm mesh file')
    end
end
if isfield(f, 'pos')
    pos = f.pos;
    tri = f.tri;
elseif isfield(f, 'Vertices')
    pos = f.Vertices;
    tri = f.Faces;
elseif isfield(f, 'vertices')
    pos = f.vertices;
    tri = f.faces;
elseif isfield(f, 'vol')
    pos = f.vol.bnd(3).pnt;
    tri = f.vol.bnd(3).tri;
end

if size(pos,1) == 3
    pos = pos';
end

if nargin > 1 && ~isempty(transform)
    pos = traditionaldipfit(transform)*[pos ones(size(pos,1),1)]';
    pos(4,:) = [];
    pos = pos';
end
if nargin > 1
    xl = diff(xlim);
    p = max(max(pos));
    if p > xl*100 || xl > p*100
        disp('Warning: widely different scale, one of the mesh might not be visible');
    end
end
if nargin < 3
    color1  = [1 1 1]*0.5;
end

pairwiseDist = ones(size(pos,1),4);
colors = pairwiseDist(:,4)*color1; % some initial color (erased below)

% fig = figure; 
hold on
fig = gcf;
patch('Faces',tri,'Vertices',pos, 'FaceVertexCdata',colors,'facecolor','interp','edgecolor','none', 'facealpha', 0.5);
userDat.axisBrain = gca;
axis equal;
axis off;
hold on;

% change lighting
set(fig, 'renderer', 'opengl');
lighting(userDat.axisBrain, 'phong');
hlights = findobj(userDat.axisBrain,'type','light');
delete(hlights)
hlights = [];
camlight(0,0);
camlight(90,0);
camlight(180,0);
camlight(270,0);
camproj orthographic
axis vis3d
camzoom(1);