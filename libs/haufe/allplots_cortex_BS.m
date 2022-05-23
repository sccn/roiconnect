function allplots_cortex_BS(cortex, data_in, colorlimits, cm, unit, smooth, printfolder, varargin)
% (C) 2018 Stefan Haufe
%
% If you use this code in a publication, please cite
%
% Haufe, S., & Ewald, A. (2016). A simulation framework for benchmarking
% EEG-based brain connectivity estimation methodologies. Brain topography, 1-18.

if nargin < 7
    printfolder = '';
end

if length(data_in) == size(cortex.Vertices, 1)
    data = data_in;
else
    % find Atlas with the same number of ROIs
    for iatl = 1:length(cortex.Atlas)
        if length(data_in) == length(cortex.Atlas(iatl).Scouts)
            data = nan*ones(size(cortex.Vertices, 1), 1);
            for iroi = 1:length(cortex.Atlas(iatl).Scouts)
                data(cortex.Atlas(iatl).Scouts(iroi).Vertices) = data_in(iroi);
            end
            break;
        end
    end
end

cortex.Vertices = cortex.Vertices(:, [2 1 3]);
cortex.Vertices(:, 1) = -cortex.Vertices(:, 1);


for iatl = 1:length(cortex.Atlas)
    if isequal(cortex.Atlas(iatl).Name, 'Structures')
        break;
    end
end

cortex.Faces_left = cortex.Faces;
%cortex.Faces_left(min(ismember(cortex.Faces_left, cortex.Atlas(iatl).Scouts(1).Vertices), [], 2) == 0, :) = [];

cortex.Faces_right = cortex.Faces;
%cortex.Faces_right(min(ismember(cortex.Faces_right, cortex.Atlas(iatl).Scouts(2).Vertices), [], 2) == 0, :) = [];

set(0,'DefaultFigureColor',[1 1 1])
if ~isempty(printfolder)
    printfolder = [printfolder '/'];
    mkdir(printfolder)
end
res = '150';

if smooth
    SurfSmoothIterations = ceil(300 * smooth * length(cortex.Vertices) / 100000);
    vc = tess_smooth(cortex.Vertices, 1, SurfSmoothIterations, tess_vertconn(cortex.Vertices, cortex.Faces), 1);
    sm = '_smooth';
else
    vc = cortex.Vertices;
    sm = '';
end

surface_pars = struct('alpha_const', 1, 'colormap', cm, 'colorlimits', colorlimits, ...
    'showdirections', 0, 'colorbars', 0, 'dipnames', [], 'mymarkersize', 15, 'directions', [0 0 1 1 1 1], ...
    'printcbar', 1, 'userticks', []);

if length(varargin) > 0
    varargin1 = varargin{1};
else
    varargin1 = {};
end

if length(varargin) > 1
    input_pars = varargin{2};
    finames = fieldnames(input_pars);
    for ifi = 1:length(finames)
        surface_pars = setfield(surface_pars, finames{ifi}, getfield(input_pars, finames{ifi}));
    end
end

figure('position', [60 828 950 500]);
subplot(2,3,1);
surface_pars.myviewdir = [-1 0 0];
showsurface3(vc, cortex.Faces_left, surface_pars, data, varargin1{:});
%export_fig([printfolder 'cortex' sm '_left'], ['-r' num2str(res)], '-a2', '-transparent');

%subplot(2,3,2);
%showsurface3(vc, cortex.Faces_right, surface_pars, data, varargin1{:});
%export_fig([printfolder 'cortex' sm '_right_inner'], ['-r' num2str(res)], '-a2', '-transparent');

surface_pars.myviewdir = [1 0 0];
subplot(2,3,2);
showsurface3(vc, cortex.Faces_right, surface_pars, data, varargin1{:});
%export_fig([printfolder 'cortex' sm '_right'], ['-r' num2str(res)], '-a2', '-transparent');

%subplot(2,3,4);
%showsurface3(vc, cortex.Faces_left, surface_pars, data, varargin1{:});
%export_fig([printfolder 'cortex' sm '_left_inner'], ['-r' num2str(res)], '-a2', '-transparent');

surface_pars.myviewdir = [-1e-10 0 1];
surface_pars.directions = [1 0 1 1 0 0];

% figure; showsurface3(vc, sa.cortex75K.tri, surface_pars, data, varargin1{:});
% export_fig([printfolder 'cortex' sm '_top'], ['-r' num2str(res)], '-a2', '-transparent');

surface_pars.myviewdir = [0 0 1];

subplot(2,3,4);
showsurface3(vc, cortex.Faces, surface_pars, data, varargin1{:});
%export_fig([printfolder 'cortex' sm '_top_upright'], ['-r' num2str(res)], '-a2', '-transparent');


% surface_pars.myviewdir = [-1e-10 0 -1];
%
% figure; showsurface3(vc, sa.cortex75K.tri, surface_pars, data, varargin1{:});
% export_fig([printfolder 'cortex' sm '_bottom'], ['-r' num2str(res)], '-a2', '-transparent');

surface_pars.myviewdir = [0 1e-10 -1];

subplot(2,3,5);
showsurface3(vc, cortex.Faces, surface_pars, data, varargin1{:});
% export_fig([printfolder 'cortex' sm '_bottom_upright'], ['-r' num2str(res)], '-a2', '-transparent');

if isfield(surface_pars, 'printcbar') && surface_pars.printcbar
    h = subplot(2,3,3);
    hf = imagesc(randn(5)); colormap(cm)
    set(h, 'clim', colorlimits, 'visible', 'off'); %, 'position', [0.1 0.1 0.6 0.8]
    set(hf, 'visible', 'off')
    cb = colorbar;
    set(cb, 'fontsize', 30)
    if ~isempty(surface_pars.userticks)
        set(cb, 'xtick', sort([colorlimits, surface_pars.userticks]))
    end
    ylabel(cb, unit)
    %export_fig([printfolder 'cortex_cbar'], ['-r' num2str(res)], '-a2', '-transparent')
end
if ~isempty(printfolder)
    export_fig([printfolder 'cortex' sm ''], ['-r' num2str(res)], '-a2', '-transparent');
end

% set(0,'DefaultFigureColor','remove')
% figure;
% hf = imagesc(randn(5)); colormap(cm)
% set(gca, 'clim', colorlimits, 'position', [0.1 0.1 0.6 0.8], 'visible', 'off')
% set(hf, 'visible', 'off')
% cb = colorbar;
% set(cb, 'fontsize', 30)
% ylabel(cb, unit)
% export_fig([printfolder 'slices_cbar'], ['-r' num2str(res)], '-a2')

% set(0,'DefaultFigureColor',[1 1 1])
%
% figure; showmri_transp3(sa.mri, struct('orientation', 'axial', ...
%     'trcut', 0.5, 'colormaps', {{cm}}, 'colorbars', 0, 'colorlimits', colorlimits), [sa.cortex75K.vc(sa.cortex5K.in_from_cortex75K, :) data(sa.cortex5K.in_from_cortex75K)], varargin1{:})
% export_fig([printfolder 'slices_axial'], '-r300', '-a2');
%
% figure; showmri_transp3(sa.mri, struct('orientation', 'sagittal', ...
%     'trcut', 0.5, 'colormaps', {{cm}}, 'colorbars', 0, 'colorlimits', colorlimits), [sa.cortex75K.vc(sa.cortex5K.in_from_cortex75K, :) data(sa.cortex5K.in_from_cortex75K)], varargin1{:})
% export_fig([printfolder 'slices_sagittal'], '-r300', '-a2');
%
% figure; showmri_transp3(sa.mri, struct('orientation', 'coronal', ...
%     'trcut', 0.5, 'colormaps', {{cm}}, 'colorbars', 0, 'colorlimits', colorlimits), [sa.cortex75K.vc(sa.cortex5K.in_from_cortex75K, :) data(sa.cortex5K.in_from_cortex75K)], varargin1{:})
% export_fig([printfolder 'slices_coronal'], '-r300', '-a2');

% close all











