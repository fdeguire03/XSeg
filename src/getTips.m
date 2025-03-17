function sortedTips = getTips(skeletonPoints, comparison_point)
[y, x] = find(skeletonPoints);
candidateTips = cat(2, x, y);
distances = vecnorm(candidateTips - comparison_point, 2, 2);
[sortedDistances, fingerTipIdx] = sort(distances, 'descend');
sortedTips = candidateTips(fingerTipIdx, :);
sortedTips(:, 3) = sortedDistances;
centerHand = sortedTips(end-1:end, 1:2); % remember the last two because these will be the center of the hand

% filter out points that are too close to each other

D = squareform(pdist(sortedTips));  % Distance matrix
minDist = 75;  
closePairs = D < minDist & D > 0;  % Exclude self-distances
removeIdx = false(size(sortedTips, 1), 1);  % Initialize removal mask
for i = 1:size(candidateTips, 1)
    for j = i+1:size(candidateTips, 1)
        if closePairs(i, j)
            removeIdx(j) = true;
        end
    end
end
sortedTips = sortedTips(~removeIdx, :, :);

% filter out points that have a similar angle to the comparison point

% Compute angle of each branch point relative to the palm center
angles = atan2(sortedTips(:,2) - comparison_point(2), sortedTips(:,1) - comparison_point(1));

% Cluster branch points that have similar angles (one per finger)
[sortedAngles, sortIdx] = sort(angles);
filteredBranchPoints = [];

angleThreshold = 0.075; % Adjust based on finger spread

for i = 1:length(sortedAngles)
    if isempty(filteredBranchPoints) || abs(sortedAngles(i) - filteredBranchPoints(end,end)) > angleThreshold
        filteredBranchPoints = [filteredBranchPoints; sortedTips(sortIdx(i), :), sortedAngles(i)];
    else
        if sortedTips(sortIdx(i), 3) > filteredBranchPoints(end, 3)
            filteredBranchPoints(end, :) = [sortedTips(sortIdx(i), :), sortedAngles(i)];
        end
    end
end

[~, fingerTipIdx] = sort(filteredBranchPoints(:, 3), 'descend');
sortedTips = filteredBranchPoints(fingerTipIdx, 1:2); % Remove angles and distances from final result
sortedTips = cat(1, sortedTips, centerHand);

end