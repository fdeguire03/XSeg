function isHand = checkHand(mask, solidity, bbox, shrinkage)

if shrinkage < 0.1
    isHand = false;
    return
end

if solidity < 0.5 || solidity > 0.9
    isHand = false;
    return
end

x = round(bbox(1));  % Left coordinate
y = round(bbox(2));  % Top coordinate
w_bbox = round(bbox(3));  % Width
h_bbox = round(bbox(4));  % Height

% Define a threshold for excessive foreground near the edges
edge_threshold = 0.05;  % 5% of the width/height

[rows, cols] = find(mask);
rows_start = min(rows);
h = max(rows) - rows_start;
cols_end = max(cols);
cols_start = min(cols);
w = cols_end - cols_start;

row_comparison = round(h * edge_threshold);
col_comparison = round(w * edge_threshold);

% Check for excessive mask presence in the top/right regions
top_region = sum(mask(rows_start:rows_start+row_comparison, :), 'all');  % Top 20%
right_region = sum(mask(:, cols_end - col_comparison:cols_end), 'all');  % Right 20%
left_region = sum(mask(:, cols_start:cols_start+col_comparison, :), 'all');  % Right 20%

top_fraction = top_region / (row_comparison * w);
right_fraction = right_region / (col_comparison * h);
left_fraction = left_region / (col_comparison * h);

% If a large fraction of pixels are in these regions, flag it as a false detection
if (top_fraction > 0.25) || (right_fraction > 0.25) || (left_fraction > 0.5)
    isHand = false;
    return
end

isHand = true;
return

row_disjoint = mean(rows) - median(rows)
col_disjoint = mean(cols) - median(cols)
disjointedness = max(abs(row_disjoint), abs(col_disjoint));
if disjointedness > 50
    isHand = false;
    return
end
isHand = true;
return

mask = bwconvhull(mask, 'objects'); % Ensure full convex shape
[B, L] = bwboundaries(mask, 'noholes');
if isempty(B)
    isHand = false;
    return;
end
% Count significant convex defects (finger gaps)
numFingers = length(B{1}) / 50 % Rough estimate (tune factor)

% Consider valid if 4-5 fingers are detected
isHand = (numFingers >= 3.5) && (numFingers <= 6.5);

end