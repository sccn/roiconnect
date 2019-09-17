function allplots_cortex(cortex, data_in, colorlimits, cm, unit, smooth, printfolder, varargin)

if length(data_in) == size(cortex.vc, 1)
  data = data_in;
else
  if length(data_in) == length(unique(cortex.in_HO)) || length(data_in) == length(unique(cortex.in_HO))-1
    data = nan*ones(size(cortex.vc, 1), 1);
    for iroi = 1:length(data_in)
      data(find(cortex.in_HO == iroi)) = data_in(iroi);
    end
  end
end

set(0,'DefaultFigureColor',[1 1 1])
printfolder = [printfolder '/'];
mkdir(printfolder)

res = '150';

if smooth 
    vc = cortex.vc_smooth;
    sm = '_smooth';
else
    vc = cortex.vc;
    sm = '';
end

surface_pars = struct('alpha_const', 1, 'colormap', cm, 'colorlimits', colorlimits, 'colorbars', 0, 'dipnames', [], 'directions', [0 0 1 1 1 1], 'showdirections', 0);

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
    
surface_pars.myviewdir = [-1 0 0];
figure; showsurface3(vc, cortex.tri_left, surface_pars, data, varargin1{:});
export_fig([printfolder 'cortex' sm '_left'], ['-r' num2str(res)], '-a2'); 


figure; showsurface3(vc, cortex.tri_right, surface_pars, data, varargin1{:});
export_fig([printfolder 'cortex' sm '_right_inner'], ['-r' num2str(res)], '-a2'); 


surface_pars.myviewdir = [1 0 0];

figure; showsurface3(vc, cortex.tri_left, surface_pars, data, varargin1{:});
export_fig([printfolder 'cortex' sm '_left_inner'], ['-r' num2str(res)], '-a2'); 

figure; showsurface3(vc, cortex.tri_right, surface_pars, data, varargin1{:});
export_fig([printfolder 'cortex' sm '_right'], ['-r' num2str(res)], '-a2'); 


surface_pars.myviewdir = [-1e-10 0 1];
surface_pars.directions = [1 1 1 1 0 0];

figure; showsurface3(vc, cortex.tri, surface_pars, data, varargin1{:});
export_fig([printfolder 'cortex' sm '_top'], ['-r' num2str(res)], '-a2'); 


surface_pars.myviewdir = [0 0 1];

figure; showsurface3(vc, cortex.tri, surface_pars, data, varargin1{:});
export_fig([printfolder 'cortex' sm '_top_upright'], ['-r' num2str(res)], '-a2'); 


surface_pars.myviewdir = [-1e-10 0 -1];

figure; showsurface3(vc, cortex.tri, surface_pars, data, varargin1{:});
export_fig([printfolder 'cortex' sm '_bottom'], ['-r' num2str(res)], '-a2'); 


surface_pars.myviewdir = [0 1e-10 -1];

figure; showsurface3(vc, cortex.tri, surface_pars, data, varargin1{:});
export_fig([printfolder 'cortex' sm '_bottom_upright'], ['-r' num2str(res)], '-a2'); 


figure; 
hf = imagesc(randn(5)); colormap(cm)
set(gca, 'clim', colorlimits, 'position', [0.1 0.1 0.6 0.8], 'visible', 'off')
set(hf, 'visible', 'off')
cb = colorbar; 
set(cb, 'fontsize', 30)
ylabel(cb, unit)
export_fig([printfolder 'cortex_cbar'], ['-r' num2str(res)], '-a2')  


% set(0,'DefaultFigureColor','remove')
% figure; 
% hf = imagesc(randn(5)); colormap(cm)
% set(gca, 'clim', colorlimits, 'position', [0.1 0.1 0.6 0.8], 'visible', 'off')
% set(hf, 'visible', 'off')
% cb = colorbar; 
% set(cb, 'fontsize', 30)
% ylabel(cb, unit)
% export_fig([printfolder 'slices_cbar'], ['-r' num2str(res)], '-a2')  

set(0,'DefaultFigureColor',[1 1 1])
% 
% figure; showmri_transp3(sa.mri, struct('orientation', 'axial', ...
%     'trcut', 0.5, 'colormaps', {{cm}}, 'colorbars', 0, 'colorlimits', colorlimits), [cortex.vc(sa.cortex5K.in_from_cortex75K, :) data(sa.cortex5K.in_from_cortex75K)], varargin1{:})
% export_fig([printfolder 'slices_axial'], '-r300', '-a2'); 
% 
% figure; showmri_transp3(sa.mri, struct('orientation', 'sagittal', ...
%     'trcut', 0.5, 'colormaps', {{cm}}, 'colorbars', 0, 'colorlimits', colorlimits), [cortex.vc(sa.cortex5K.in_from_cortex75K, :) data(sa.cortex5K.in_from_cortex75K)], varargin1{:})
% export_fig([printfolder 'slices_sagittal'], '-r300', '-a2'); 
% 
% figure; showmri_transp3(sa.mri, struct('orientation', 'coronal', ...
%     'trcut', 0.5, 'colormaps', {{cm}}, 'colorbars', 0, 'colorlimits', colorlimits), [cortex.vc(sa.cortex5K.in_from_cortex75K, :) data(sa.cortex5K.in_from_cortex75K)], varargin1{:})
% export_fig([printfolder 'slices_coronal'], '-r300', '-a2'); 

% close all











