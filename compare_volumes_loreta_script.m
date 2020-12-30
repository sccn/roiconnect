filename1 = '/data/matlab/eeglab/plugins/dipfit/standard_BEM/standard_vol.mat';
p  = fileparts(which('ft_defaults.m'));
filename2 = fullfile(p, 'template','atlas','afni','TTatlas+tlrc.HEAD');

% PUT BREAKPOINT LINE 44 of plot3dmeshalign (which will load and align the
% AFNI model
plot3dmeshalign(filename1, filename2, [], [1 0 0], [0 0 1], {'Thalamus'});

% then copy and paste the code below
filename2 = 'LORETA-Talairach-BAs.mat';
tmp = load('-mat', filename2);

indRm = [];
for iVert = 1:length(tmp.Vertices)
    if min(sqrt(sum(bsxfun(@minus, pos2(:,:), tmp.Vertices(iVert,:)).^2,2))) > 20
        indRm = [indRm iVert];
    end
end
tmp.Vertices(indRm,:) = [];
figure; h = plot3(pos2(:,1),pos2(:,2),pos2(:,3), '.', 'color', color2);
hold on; plot3(tmp.Vertices(:,1),tmp.Vertices(:,2), tmp.Vertices(:,3), '.', 'color', color1);
axis equal;
xlabel('x')
ylabel('y')
zlabel('z')

pos3 = [-3 -4;
        -3 -11;
        -3 -18;
        -10 -4;
        -10 -11;
        -10 -18;
        -10 -25;
        -17 -11;
        -17 -18;
        -17 -25;
        -17 -32;
        -24 -25;
        -24 -32];
pos3 = [ pos3; -pos3(:,1)+1 pos3(:,2)];
pos3 = [ [pos3 ones(size(pos3,1),1)*1];[pos3 ones(size(pos3,1),1)*8];[pos3 ones(size(pos3,1),1)*15] ];
hold on;  h = plot3(pos3(:,1),pos3(:,2), pos3(:,3), '.', 'color', 'g');
set(h, 'MarkerSize', 20)

% Add to Loreta source model
cortex = load('-mat', 'LORETA-Talairach-BAs.mat');
cortex.Vertices = [ cortex.Vertices; pos3];
cortex.Atlas.Scouts(end+1).Label = 'ThalamusL';
cortex.Atlas.Scouts(end).Vertices = [2395:2395+size(pos3,1)/2-1];
cortex.Atlas.Scouts(end+1).Label = 'ThalamusR';
cortex.Atlas.Scouts(end).Vertices = [2395+size(pos3,1)/2:size(cortex.Vertices,1)];
save('-mat', 'LORETA-Talairach-BAs.mat', '-struct', 'cortex');
