function [inds, PCA_inds] = fp_npcs2inds(npcs)

% Copyright (c) 2022 Franziska Pellegrini and Stefan Haufe

%%
beg_inds = cumsum([1 npcs(1:end-1)]);
end_inds = cumsum([npcs]);

for iroi = 1:numel(npcs)
    PCA_inds{iroi} = beg_inds(iroi):end_inds(iroi);
end

inds = {}; ninds = 0;
for iroi = 1:numel(npcs)
    for jroi = (iroi+1):numel(npcs)
        inds{ninds+1} = {PCA_inds{iroi}, PCA_inds{jroi}};
        ninds = ninds + 1;
    end
end