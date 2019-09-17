function han = showvc_prog_plain(loc,tri,view_dir,para);
% Guido Nolte, Stefan Haufe
%
% If you use this code in a publication, please ask Guido Nolte for the correct 
% citation. In addition, cite
%
% Haufe, S., & Ewald, A. (2016). A simulation framework for benchmarking 
% EEG-based brain connectivity estimation methodologies. Brain topography, 1-18.
%

cmap=1;
voxelkont=0;
colorbars=1;
alphaclimPolicy = 'none';
alpha_const = 1;
if nargin>3 
  if isfield(para, 'alpha_const')
    alpha_const = para.alpha_const;
  end
  
  if isfield(para,'voxelfield')
     voxelfield=para.voxelfield;
     voxelkont=1;
      if isfield(para,'colorlimits');
         colmin=para.colorlimits(1);
         colmax=para.colorlimits(2);
      end

    source_val = voxelfield;
    source_val_max=max(abs(source_val));
    source_val_min=min(abs(source_val));
    
  
    if isfield(para, 'climPolicy')
     climPolicy = para.climPolicy;
    else
     if ~isfield(para, 'colorlimits')
         if sign(max(voxelfield)) == sign(min(voxelfield)) || abs(min(voxelfield)) < eps 
             climPolicy = 'minmax';
         else
             climPolicy = 'sym';
         end
     else
         climPolicy = 'none';
     end
    end
     
    if isequal(climPolicy, 'minmax')
        colmax=max(source_val);
        colmin=min(source_val);
    elseif isequal(climPolicy, 'sym')
        colmax=source_val_max;
        colmin=-source_val_max;
    end 
    
    if ~isfield(para, 'colormap')
      if isequal(climPolicy, 'minmax') || sign(colmax) == sign(colmin) || abs(colmin) < eps
          load('cm8', 'cm8')
          mycolormap= cm8;
      else
          load('cm9', 'cm9')
          mycolormap= cm9;
      end
    else
      mycolormap = para.colormap;
    end
  end
  
  if isfield(para,'alphafield')
    alphafield=para.alphafield;
    alpha_val = alphafield;
    alpha_val_max=max(abs(alpha_val));
    alpha_val_min=min(abs(alpha_val));
    if isequal(alphaclimPolicy, 'minmax')
        alphamax=max(alpha_val);
        alphamin=min(alpha_val);
    elseif isequal(alphaclimPolicy, 'sym')
        alphamax=alpha_val_max;
        alphamin=-alpha_val_max;
    end
  end
 if isfield(para,'mycolormap')
     mycolormap=para.mycolormap;
 end
 if isfield(para,'colorbars')
     colorbars=para.colorbars;
 end
 if isfield(para, 'climPolicy')
    climPolicy = para.climPolicy;
 end

  if isfield(para, 'alphaclimPolicy')
    alphaclimPolicy = para.alphaclimPolicy;
  end
  
  if isfield(para, 'colormap')
    mycolormap = para.colormap;
  end
  
  
end

[ntri,ndum]=size(tri);


%colormap(newmap(30:150,:));
% brighten(-.6);


locm=mean(loc);
[nloc,ndum]=size(loc);
relloc=loc-repmat(locm,nloc,1);
dis=(sqrt(sum((relloc.^2)')))';thresh=3;dis(dis<thresh)=thresh;
%map=colormap('jet');newmap=colormap_interpol(map,3);size(newmap);colormap(
%newmap);
% cortexcolor=[255 213 119]/255;
cortexcolor=[.8 .8 .8];
% cortexcolor=[1 1 1];
% cortexcolor=[234 183 123]/255;
h=patch('vertices',loc,'faces',tri);
set(h,'FaceColor',cortexcolor);
view(view_dir);
% set(h,'edgecolor',[0 0 0], 'linewidth', 1,'facelighting','flat');
set(h,'edgecolor', 'none','facelighting','phong');
% set(h,'specularexponent', 100);
set(h,'specularstrength', 0);
set(h,'ambientstrength',.25);
set(h,'diffusestrength',.8)
% dis0=0*dis+0;
% set(h,'facevertexalphadata',dis0)
% set(h,'alphadatamapping','direct')
set(h,'facealpha', alpha_const)

material dull
camlight('headlight','infinite');
axis equal 



if voxelkont>0
   if length(voxelfield) == size(loc, 1)
       if isfield(para, 'alphafield')
           if isequal(alphaclimPolicy, 'none')
               adata = alphafield;               
           elseif isequal(alphaclimPolicy, 'minmax') || sign(alphamax) == sign(alphamin) || abs(alphamin) < eps
               adata = alphafield;
               adata = adata-min(adata);
               adata = adata/max(adata);
               adata = tanh(8*adata - 3);
               adata = adata-min(adata);
               adata = adata/max(adata);
               adata = adata*0.6;
               adata = adata + 0.4;
           elseif isequal(alphaclimPolicy, 'sym')
               adata = abs(alphafield);
               adata = adata/max(adata);
               adata = tanh(8*adata - 3);
               adata = adata-min(adata);
               adata = adata/max(adata);   
               adata = adata*0.6;
               adata = adata + 0.4;
           end           
       else
         if isfield(para, 'alpha_const')
           adata = para.alpha_const*ones(size(voxelfield));
         else
           if isequal(climPolicy, 'minmax') || sign(colmax) == sign(colmin) || abs(colmin) < eps
               adata = voxelfield;
               adata = adata-min(adata);
               adata = adata/max(adata);
               adata = tanh(8*adata - 3);
               adata = adata-min(adata);
               adata = adata/max(adata);
               adata = adata*0.7;
               adata = adata + 0.3;
           else
               adata = abs(voxelfield);
               adata = adata/max(adata);
               adata = tanh(8*adata - 3);
               adata = adata-min(adata);
               adata = adata/max(adata);
               adata = adata*0.7;
               adata = adata + 0.3;
           end
         end
       end
       
       voxelfield(voxelfield > colmax) = colmax;
       voxelfield(voxelfield < colmin) = colmin;
       
       map=colormap(mycolormap);
       newmap=colormap_interpol(map,3);
       
       col_inds = floor((voxelfield-colmin)/((colmax-colmin)/(size(newmap, 1)-1))+1);
       face_cols = zeros(size(col_inds, 1), 3);
       face_cols(~isnan(col_inds), :) = newmap(col_inds(~isnan(col_inds)), :);
       face_cols(isnan(col_inds), :) = repmat([0.8 0.8 0.8], size(find(isnan(col_inds)), 1), 1);
       
       face_cols(isnan(voxelfield), :) = 0.4;
       
       set(h, 'facealpha', 'interp', 'facevertexalphadata', adata, 'alphadatamapping', 'none', 'FaceVertexCData', face_cols, 'facecolor', 'interp')
%        'facealpha', 'interp', 'facevertexalphadata', adata, 'alphadatamapping', 'none', 
       alim([0 1]);
%        caxis([colmin colmax])
%        colormap(newmap);
%        keyboard
   else
       set(h, 'FaceVertexCData',voxelfield, 'facecolor', 'flat')
       map=colormap(mycolormap);
       newmap=colormap_interpol(map,3);size(newmap);colormap(newmap);
   end
   
%      vf=voxelfield(abs(voxelfield)>=0);
% load cm5
%      map=colormap(mycolormap);
%      newmap=colormap_interpol(map,3);size(newmap);colormap(newmap);
%      [nc,nx]=size(newmap);
%      vfint=ceil((nc-1)*((voxelfield-min(vf))/(max(vf)-min(vf))));
%      vfint(voxelfield==0)=1;
%      vftruecolor=newmap(vfint);
   %h=patch('vertices',loc,'faces',tri_new,'FaceVertexCData',voxelfield,...
   % 'facecolor','interp','edgecolor','none','facelighting','phong');
%     h=patch('vertices',loc,'faces',tri_new,'FaceVertexCData',vftruecolor,...
%     'facecolor','interp','edgecolor','none','facelighting','phong');
   
   if colorbars==1  
      pos=get(gca,'pos');
      set(gca,'pos',[pos(1) pos(2) 0.8*pos(3) pos(4)]);
      pos=get(gca,'pos');
      han.cb = colorbar('location','eastoutside', 'position', [pos(1)+pos(3)+0.1 pos(2)+0.05 0.05 pos(4)-0.1]);
   end
   
end

% set(h,'specularexponent',50000);
han.patch = h;

return;