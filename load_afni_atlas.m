function [newBrick, xyz, labels, labelsstr ] = load_afni_atlas(sourcemodel, headmodel, sourcemodel2mni)

    afni = ft_read_atlas(sourcemodel);
    mri = afni.brick1(:,:,:,1);
    labelsstr = afni.brick1label;
    [r,c,v] = ind2sub(size(mri),find(mri));
    
    sub = 4; % Subsample by 4
    newBrickCell = cell(ceil(size(mri)/sub));
    for iDip = 1:length(r)
        newBrickCell{round(r(iDip)/sub), round(c(iDip)/sub), round(v(iDip)/sub)} = [ newBrickCell{round(r(iDip)/sub), round(c(iDip)/sub), round(v(iDip)/sub)} mri(r(iDip),c(iDip),v(iDip)) ];
    end
    
    % get dominating label for voxel
    newBrick = zeros(ceil(size(mri)/sub));
    for iDip = 1:length(newBrickCell(:))
        if ~isempty(newBrickCell{iDip})
            uniq = unique(newBrickCell{iDip});
            if length(uniq) > 1
                count = histc(newBrickCell{iDip}, uniq);
                [~,indMax] = max(count);
                newBrick(iDip) = uniq(indMax);
            else
                newBrick(iDip) = uniq;
            end
        end
    end
    
    % recompute coordinates
    inds = find(newBrick);
    [r,c,v] = ind2sub(size(newBrick),find(newBrick));
    r = (r-1)*4 + 2.5;
    c = (c-1)*4 + 2.5;
    v = (v-1)*4 + 2.5;
    labels = newBrick(inds);
    
    % transform coordinates
    xyz = [r c v ones(length(r),1)];
    xyz = afni.transform*xyz';
    if nargin > 2 && ~isempty(sourcemodel2mni)
        xyz = traditionaldipfit(sourcemodel2mni)*xyz;
    end
    xyz(4,:) = [];
    xyz = xyz';
    
    headmodel = load('-mat', headmodel);
    try 
        inside = ft_inside_headmodel(xyz, headmodel);
    catch
        % ft_inside_headmodel not compatible with headmodel, using custom code
        p = pwd;
        cd('/data/matlab/eeglab/plugins/fieldtrip/forward/private')
        inside = bounding_mesh(xyz, headmodel.vol.bnd(end).pnt, headmodel.vol.bnd(end).tri);
        cd(p);
    end
    inside = inside > 0;

    if 0
        figure; plot3dmeshalign(headmodel);
        plot3(xyz(inside,1),xyz(inside,2),xyz(inside,3), 'b.');
        hold on; plot3(xyz(~inside,1),xyz(~inside,2),xyz(~inside,3), 'r.');
    end

    xyz(~inside,:) = [];
    labels(~inside) = [];
    