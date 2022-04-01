function plotconnectivitymultiple(networkfile, connectSpecSelect, titl, filename)

% decode network param
if ischar(networkfile)
    networkfile = load('-mat', networkfile);
end
loreta_Networks = networkfile.loreta_Networks;
loreta_ROIS     = networkfile.loreta_ROIS;

figure('position', [100 100 1000 700], 'paperpositionmode', 'auto');
ncol = ceil(sqrt(length(loreta_Networks)));
nrow = ceil(length(loreta_Networks)/ncol);

for iNet = 1:length(loreta_Networks)
    subplot(nrow, ncol, iNet)
    labels = { loreta_ROIS(loreta_Networks(iNet).ROI_inds).Label };
    plotconnectivity(connectSpecSelect{iNet}(:,:), 'labels', labels, 'axis', gca, 'threshold', 0.1);
    h = title(loreta_Networks(iNet).name, 'interpreter', 'none');
    pos = get(h, 'position');
    set(h, 'position', pos + [0 0.1 0]);
end

h = textsc(titl, 'title');
set(h, 'fontsize', 16, 'fontweight', 'bold');
print('-djpeg', filename);
close
