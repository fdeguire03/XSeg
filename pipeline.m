clear;
clc;
d = readmatrix('boneage-training-dataset.csv');
p = 1;
ages = [];
genders = [];
bone_props = [];
centerOfMasses = [];
aspectRatios = [];
centerOfMassesHand = [];
aspectRatiosHand = [];
make_plot = 0;
sample = 0;
for index = 1000:16000
%while p < 5
    try
        %index = [6561 9765 10643 7760];
        %index = index(p);
        if sample
            index = randsample(d(:,1), 1);
            i_str = num2str(index);
            f = ['boneage-training-dataset/' i_str '.png'];
            attempts = 0;
            while attempts < 5 && ~isfile(f)
                index = index + 1;
                i_str = num2str(i);
                %while length(i_str) < 5
                %    i_str = ['0' i_str];
                %end
                f = ['boneage-training-dataset/' i_str '.png'];
                attempts = attempts + 1;
            end
        else
            i_str = num2str(index);
            f = ['boneage-training-dataset/' i_str '.png'];
        end
        f
        im = mat2gray(imread(f));
        info = d(d(:,1) == index, :);
        gender = info(end);
        if gender
            gender = 'M';
        else
            gender = 'F';
        end
        age = info(2);
        test = im;
        [rows, cols] = size(im);
        for m = 1:2
            if m == 1
                lower_bound = 25;
                alpha = 0.15;
                g_size = 20;
                mult = 0.2;
                size_thresh_mult = 0.45;
                I = imadjust(test, [prctile(test, lower_bound, 'all'), 1]);
            else
                lower_bound = 1;
                alpha = 0.2;
                g_size = 10;
                mult = 0.5;
                size_thresh_mult = 0.35;
                %I = adapthisteq(test);
                I = imadjust(test, [prctile(test(test>0), lower_bound, 'all'), 1]);
            end
            [hand_mask, area, tim, solidity, bbox, orientation] = getHandMask(I, g_size, mult, size_thresh_mult, alpha, 1);
            test = test .* imdilate(hand_mask, strel('disk', 20));
            if m == 1
                last_area = area;
                last_solidity = solidity;
                last_bbox = bbox;
                last_orientation = orientation;
                tim1 = tim;
                hand_mask1 = hand_mask;
            end
        end

        %{
        angle = orientation - 90;
        
        if abs(angle) > 15
            im = imrotate(im, angle, 'nearest', 'crop');
            hand_mask = imrotate(hand_mask, angle, 'nearest', 'crop');
            hand_mask1 = imrotate(hand_mask1, angle, 'nearest', 'crop');
            tim1 = imrotate(tim1, angle, 'nearest', 'crop');
        end
        %}


        [E, t] = edge(tim1, "Canny");
        [E2, t] = edge(hand_mask1, "Canny");
        [E3, t] = edge(hand_mask, "Canny");
        
        se = strel('disk', 3);
        thick_edges = imdilate(E, se);
        thick_edges2 = imdilate(E2, se);
        thick_edges3 = imdilate(E3, se);
        
        %I = hand_mask .* I
        shrinkage = area / last_area;
        isHand = checkHand(hand_mask, solidity, bbox, shrinkage);
        if ~isHand
            isHand = checkHand(hand_mask1, last_solidity, last_bbox, 1);
            if ~isHand
                disp('Hand mask not generated with high enough confidence, skipping...')
                continue
            end
            temp = thick_edges3;
            thick_edges3 = thick_edges2;
            thick_edges2 = temp;
            hand_mask = hand_mask1;
        end
        disp('Found hand mask')

        [rows, cols] = find(hand_mask);

        % Get the leftmost and rightmost edge at each row
        minCols = accumarray(rows, cols, [], @min);  % Left edge
        maxCols = accumarray(rows, cols, [], @max);  % Right edge
        width = maxCols - minCols;
        baseWidth = mean(width(end-15:end-5));
        handIdx = find(width > 1.1*baseWidth, 1, 'last');
        handIdx = round(0.8*handIdx + 0.2*max(rows));  % go 20% back towards the end of the rows
        noWristMask = hand_mask;
        noWristMask(handIdx:end, :) = 0;

        hand_mask = noWristMask;
        [rows, cols] = find(hand_mask);
        row_max = max(rows);
        bottom_center = [median(cols(rows == row_max)) row_max];
        skeleton = bwmorph(hand_mask, 'skel', Inf);
        
        branchPoints = bwmorph(skeleton, 'branchpoints');
        tipPoints = getTips(branchPoints, bottom_center);
        centerHand = tipPoints(end-1:end, :);
        if length(tipPoints) > 6
            fingerTips = tipPoints(1:5, :);
        else
            endPoints = bwmorph(skeleton, 'endpoints');
            tipPoints = getTips(endPoints, bottom_center);
            fingerTips = tipPoints(1:min(length(tipPoints), 5), :);
        end
        skeleton = imdilate(skeleton, strel('disk', 4));
        
        keyPoints = cat(1, fingerTips, centerHand);

        boneMask = zeros(size(hand_mask));
        pixelThresholds = zeros(size(hand_mask));
        thresholdCounts = zeros(size(hand_mask));
        segments = getSegmentPaths(keyPoints(1:6, :), skeleton, hand_mask, [250 175]);
        test = im .* hand_mask;
        %imshow(test);
        for i = length(segments):-1:1
            s = segments{i};
            for j = 1:length(s)
                patchMask = s{j}{1};
                yRange = s{j}{2};
                xRange = s{j}{3};
                region = test(yRange, xRange);
                if i > 8
                    rectangle('Position', [min(xRange) min(yRange) max(xRange)-min(xRange) max(yRange)-min(yRange)], 'EdgeColor', 'red')
                end
                if i < 6
                    tr = medfilt2(adapthisteq(region), [5 5]);
                    %tr = imgaussfilt(tr, 2.5);
                    %tr = medfilt2(region, [9 9]);
                    %tr = imgaussfilt(adapthisteq(region), 1.5);
                    mult = 0.6;
                    mult2 = 1;
                else
                    tr = imgaussfilt(adapthisteq(region), 2.5);
                    %tr = imgaussfilt(region, 2.5);
                    mult = 0.45;
                    mult2 = 1;
                end
                t = graythresh(tr(patchMask>0))*mult2;
                mask = tr > t;
                %mask = imclose(mask, strel('disk', 4));
                %mask = imfill(mask, 6, 'holes');
                attempts = 0;
                last_t = 0;
                while attempts < 5 && sum(mask(:)) > mult * sum(patchMask(:))
                    t = graythresh(tr(tr > last_t))*mult2;
                    if t - last_t > 0.01
                        t = 0.6*t + 0.4*last_t;  % prevent threshold from growing too quickly
                    end
                    
                    %t = graythresh(tr(tr > last_t));
                    mask = tr > t;
                    attempts = attempts + 1;
                    last_t = t;
                    %mask = imclose(mask, strel('disk', 4));
                    %mask = imfill(mask, 6, 'holes');
                end
                %mask = imclose(mask, strel('disk', 4));
                mask = imfill(mask, 6, 'holes');
                %imshow(region)
                %hold on;
                %overlayColor = cat(3, ones(size(mask)), ones(size(mask)), zeros(size(mask))); 
                %h = imshow(overlayColor);
                %set(h, 'AlphaData', 0.1 * double(mask));
                boneMask(yRange, xRange) = max(boneMask(yRange, xRange), mask);
                pixelThresholds(yRange, xRange) = pixelThresholds(yRange, xRange) + t;
                thresholdCounts(yRange, xRange) = thresholdCounts(yRange, xRange) + 1;
                
            end
            
        end
        pixelThresholds = pixelThresholds ./ max(1, thresholdCounts);
        pixelThresholds(pixelThresholds == 0) = 1;
        boneMask2 = adapthisteq(test) > pixelThresholds;
        boneMask2 = imclose(boneMask2, strel('disk', 3));
        boneMask2 = imfill(boneMask2, 6, 'holes');
        
        combinedMask = (boneMask + boneMask2) / 2;
        combinedMask = imfill(combinedMask, 'holes');

        edges = edge(adapthisteq(im), "Canny", [0.005 0.05], 2.5) & ~thick_edges3;%, [0.01 0.05], 5);
        %edges = bwareaopen(edges, 25);
        edges = imclose(edges, strel('disk', 4));
        edges = bwareaopen(edges, 200);
        edges = imfill(edges, 'holes');
        finalBoneMask = bwareaopen(edges .* combinedMask, 100);
        

        boneProportion = sum(finalBoneMask(:)) / sum(noWristMask(:));
        [xGrid, yGrid] = meshgrid(1:size(finalBoneMask, 2), 1:size(finalBoneMask, 1)); % Create coordinate grids

        % Compute weighted sums
        totalWeight = sum(finalBoneMask(:));
        xCenter = sum(finalBoneMask(:) .* xGrid(:)) / totalWeight;
        yCenter = sum(finalBoneMask(:) .* yGrid(:)) / totalWeight;
        [rows, cols] = find(finalBoneMask);
        height = (max(rows) - min(rows));
        width = (max(cols) - min(cols));
        yWeight = (yCenter - min(rows)) / height;
        xWeight = (xCenter - min(cols)) / width;
        centerOfMasses(p,:) = [xWeight yWeight];
        aspectRatios(p) = height / width;

        totalWeight = sum(hand_mask(:));
        xCenter = sum(hand_mask(:) .* xGrid(:)) / totalWeight;
        yCenter = sum(hand_mask(:) .* yGrid(:)) / totalWeight;
        [rows, cols] = find(hand_mask);
        height = (max(rows) - min(rows));
        width = (max(cols) - min(cols));
        yWeight = (yCenter - min(rows)) / height;
        xWeight = (xCenter - min(cols)) / width;
        centerOfMassesHand(p,:) = [xWeight yWeight];
        aspectRatiosHand(p) = height / width;

        if make_plot
            subplot(2,2,p)
            edge_overlay = imoverlay(im, finalBoneMask, 'green');
            outline =  bwboundaries(hand_mask);
            outline = outline{1};
            imshow(edge_overlay);
            hold on;
            plot(outline(:,2), outline(:,1), 'y');
            overlayColor = cat(3, ones(size(boneMask)), ones(size(boneMask)), zeros(size(boneMask))); 
            h = imshow(overlayColor);
            set(h, 'AlphaData', 0.25 * double(combinedMask));

            s = sprintf('ID %d (age %d, gender %s): bone prop %f', index, age, gender, boneProportion);
            title(s);
        end
        ages(p) = age;
        genders(p) = gender;
        bone_props(p) = boneProportion;
        p = p+1
        



    catch ME
        fprintf('Error: %s\n', ME.message);
    end
end

save('data.mat', 'ages', 'genders', 'bone_props', 'centerOfMasses', 'aspectRatios', 'centerOfMassesHand', 'aspectRatiosHand')

binSize = 18;
binEdges = min(ages)-1: binSize : max(ages)+1;  % bin every six months together
binCenters = binEdges(1:end-1) + binSize/2;
binnedData = struct();
for i = 1:length(binCenters)
    binnedData(i).ageRange = [binEdges(i), binEdges(i+1)];
    binIndices = (ages >= binEdges(i) & ages <= binEdges(i+1));
    binnedData(i).boneProportion = bone_props(binIndices);
    binnedData(i).xCenter = centerOfMasses(binIndices, 1);
    binnedData(i).yCenter = centerOfMasses(binIndices, 2);
    binnedData(i).aspectRatio = aspectRatios(binIndices);
    binnedData(i).xCenterHand = centerOfMassesHand(binIndices, 1);
    binnedData(i).yCenterHand = centerOfMassesHand(binIndices, 2);
    binnedData(i).aspectRatioHand = aspectRatiosHand(binIndices);
    binnedData(i).genders = genders(binIndices);
end

meanBoneProportion = zeros(length(binCenters), 1);
meanBoneProportionM = zeros(length(binCenters), 1);
meanBoneProportionF = zeros(length(binCenters), 1);
stdBoneProportion = zeros(length(binCenters), 1);
stdBoneProportionM = zeros(length(binCenters), 1);
stdBoneProportionF = zeros(length(binCenters), 1);
for i = 1:length(binCenters)
    binnedBoneProps = binnedData(i).boneProportion;
    meanBoneProportion(i) = mean(binnedBoneProps);
    meanBoneProportionM(i) = mean(binnedBoneProps(binnedData(i).genders == 'M'));
    meanBoneProportionF(i) = mean(binnedBoneProps(binnedData(i).genders == 'F'));
    stdBoneProportion(i) = std(binnedBoneProps);
    stdBoneProportionM(i) = std(binnedBoneProps(binnedData(i).genders == 'M'));
    stdBoneProportionF(i) = std(binnedBoneProps(binnedData(i).genders == 'F'));
end

figure;
scatter(binCenters, meanBoneProportion);
hold on;
scatter(binCenters, meanBoneProportionM, 'b', 'filled');
scatter(binCenters, meanBoneProportionF, 'r', 'filled');
%errorbar(binCenters, meanBoneProportionM, stdBoneProportionM)
%errorbar(binCenters, meanBoneProportionF, stdBoneProportionF)
xlabel('Age (months)');
ylabel('Average Bone Proportion');
title('Bone Proportion vs. Age');
grid on;