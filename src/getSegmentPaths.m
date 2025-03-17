function segments = getSegmentPaths(keyPoints, skeleton, handMask, patchSize)

    covered = zeros(size(handMask));
    testSpace = zeros(size(handMask));
    distMaps = zeros(size(keyPoints, 1), size(skeleton, 1), size(skeleton, 2));
    for i = 1:size(keyPoints, 1)
        % Compute geodesic distance from fingertip along skeleton
        distMap = bwdistgeodesic(skeleton, keyPoints(i,1), keyPoints(i,2), 'quasi-euclidean');
        if i > 5
            distMap = distMap * 1.5;  % penalize the center of the hand so we follow patches along the fingers for longer
        end
        distMaps(i, :, :) = distMap;
    end
    [~, Idx] = min(distMaps, [], 1);
    Idx = squeeze(Idx) .* skeleton;
    segments = cell(size(distMaps, 1), 1);

    num_patches = 1;
    for i = 1:max(Idx, [], 'all')
        [row, col] = find(Idx == i);
        bbox = zeros(4,1);
        bbox(1:2) = [min(col) min(row)];
        width = (max(col) - bbox(1))*1.5;
        height = (max(row) - bbox(2))*1.25;
        bbox(3) = width;
        bbox(4) = height;
        bbox(1) = floor(bbox(1) - width/6);
        bbox(2) = floor(bbox(2) - height/8);

        if i == size(distMaps, 1)
            patchSize = [patchSize(1)*1.5 patchSize(2)*1.5];
        end

        numPatchesHeight = ceil(height / patchSize(1));
        numPatchesWidth = ceil(width / patchSize(2));
        patchSizeHeight = ceil(height / numPatchesHeight);
        patchSizeWidth = ceil(width / numPatchesWidth);
        patches = {};
        count = 1;
        for j = 1:numPatchesHeight
            for k = 1:numPatchesWidth
                centerY = bbox(2) + (j-1)*patchSizeHeight + floor(patchSizeHeight/2);
                centerX = bbox(1) + (k-1)*patchSizeWidth + floor(patchSizeWidth/2);

                yRange = max(centerY - ceil(1.2*patchSizeHeight/2), 1) : min([centerY + ceil(1.2*patchSizeHeight/2), bbox(2)+bbox(4), size(handMask, 1)]);
                xRange = max(centerX - ceil(1.2*patchSizeWidth/2), 1) : min([centerX + ceil(1.2*patchSizeWidth/2), bbox(1)+bbox(3), size(handMask, 2)]);

                patchHeight = max(yRange) - min(yRange);
                patchWidth = max(xRange) - min(xRange);
                %testSpace(yRange, xRange) = 1;
                %testSkeletonOverlap = testSpace & skeleton;
                %sum(testSkeletonOverlap(:))
                if sum(max(handMask(yRange, xRange) - covered(yRange, xRange), 0), 'all') / (patchHeight * patchWidth) > 0.1%  && max(covered(yRange, xRange), [], 'all') < 3%&& sum(testSkeletonOverlap(:)) > 1% at least 10% of pixels should be new
                    patch = cell(3, 1);
                    patch{1} = handMask(yRange, xRange);
                    patch{2} = yRange;
                    patch{3} = xRange;
                    patches{count} = patch;
                    covered(yRange, xRange) = covered(yRange, xRange) + 1;
                    count = count + 1;
                    num_patches = num_patches + 1;
                end
                testSpace(yRange, xRange) = 0;
                
            end
        end
        segments{i} = patches;
    end
    %{
    imshow(covered)
    hold on;
    overlayColor = cat(3, ones(size(handMask)), ones(size(handMask)), zeros(size(handMask))); 
    h = imshow(overlayColor);
    set(h, 'AlphaData', 0.5 * double(handMask));
    overlayColor = cat(3, ones(size(handMask)), zeros(size(handMask)), zeros(size(handMask))); 
    h = imshow(overlayColor);
    set(h, 'AlphaData', 0.5 * double(skeleton));
    %}
end