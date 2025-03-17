function [handMask, area, tim, selected_solidity, selected_bbox, selected_orientation] = getHandMask(I, g_size, init_mult, size_thresh_mult, alpha, choose_single)
    
    [rows, cols] = size(I);
    tukeyRow = tukeywin(rows, alpha);  % Row-wise window
    tukeyCol = tukeywin(cols, alpha)'; % Column-wise window (transpose)
    tukey2D = tukeyRow * tukeyCol; 
    
    % Apply the window to the image
    I = I .* tukey2D;
    gim = imgaussfilt(I, g_size);
    %gim = medfilt2(I);
    
    t = graythresh(gim(gim>0));
    tim = gim > (t * init_mult);
    size_thresh = numel(tim) * size_thresh_mult;
    attempts = 0;
    while attempts < 5 && (sum(tim, 'all') > size_thresh)
        init_mult = init_mult + 0.05;
        tim = gim > (t * init_mult);
        attempts = attempts + 1;
    end

    
    % zero out a row near the bottom to avoid connecting the entire hand to
    % the edge
    tim(end-floor(0.05*size(I,1)):end, :) = 0;
    

    if choose_single
        CC = bwconncomp(tim);
        stats = regionprops(CC, 'Area', 'BoundingBox', 'Solidity', 'Orientation');
        
        valid_idx = [];
        max_area = 0;
        selected_idx = 0;
        for i = 1:length(stats)
            bbox = stats(i).BoundingBox;
            solidity = stats(i).Solidity;
            aspect_ratio = bbox(3) / bbox(4);  % width/height
        
            if solidity > 0.2 && solidity < 0.9
                valid_idx = [valid_idx, i];  % Store valid component indices
        
                % Choose the largest valid component
                if stats(i).Area > max_area
                    selected_solidity = solidity;
                    selected_bbox = bbox;
                    selected_orientation = stats(i).Orientation;
                    max_area = stats(i).Area;
                    selected_idx = i;
                end
            end
        end
        
        handMask = false(size(I));
        handMask(CC.PixelIdxList{selected_idx}) = true;
        area = max_area;
    else
        area = 0;
        selected_solidity = 0;
        selected_bbox = 0;
        selected_orientation = 0;
        handMask = tim;
    end