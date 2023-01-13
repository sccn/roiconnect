% distanceMat = neighbormat( EEG(1).roi.cortex.Vertices, EEG(1).roi.cortex.Atlas.Scouts, 'distance', 50); figure; imagesc(distanceMat);
% distanceMat = neighbormat( EEG(1).roi.cortex.Vertices, EEG(1).roi.cortex.Atlas.Scouts, 'pairwise', 15); figure; imagesc(distanceMat);
function distanceMat = neighbormat(vertices, scouts, method, threshold)

if nargin < 3
    method = 'distance';
    threshold = 40;
end
distanceMat = ones(length(scouts), length(scouts))*1000;

if strcmpi(method, 'distance')
    % compute center
    centers = zeros(3, length(scouts));
    for iScout = 1:length(scouts)
        centers(:,iScout) = mean(vertices(scouts(iScout).Vertices,:),1);
    end
    
    for iScout1 = 1:length(scouts)
        for iScout2 = iScout1+1:length(scouts)
            distanceMat(iScout1, iScout2) = sqrt(sum((centers(:,iScout1) - centers(:,iScout2)).^2));
            distanceMat(iScout2, iScout1) = sqrt(sum((centers(:,iScout1) - centers(:,iScout2)).^2));
        end
    end
    distanceMat = distanceMat < threshold;
else
    verticesTmp = cell(1, length(scouts));
    for iScout = 1:length(scouts)
        verticesTmp{iScout} = vertices(scouts(iScout).Vertices,:);
    end
    
    for iScout1 = 1:length(scouts)
        for iScout2 = iScout1+1:length(scouts)
            minDist = zeros(1, length(verticesTmp{iScout1}));
            for iVertice = 1: length(verticesTmp{iScout1})
                diffsq = bsxfun(@minus, verticesTmp{iScout2}, verticesTmp{iScout1}(iVertice,:)).^2;
                minDist(iVertice) = min(sqrt(sum(diffsq,2)));
            end
            distanceMat(iScout1, iScout2) = min(minDist);
            distanceMat(iScout2, iScout1) = distanceMat(iScout1, iScout2);
        end
    end
    distanceMat = distanceMat < threshold;
end

