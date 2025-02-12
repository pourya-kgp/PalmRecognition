% ----------------------------------------------------------------------------------------------------
% Author        : Pourya Khodagholipour (P.KH)
% Project Name  : Implementation of a Palmprint Recognition System based on Image Processing (M.S. Thesis)
% File Name     : Pretangent.m
% Description   : Processes a palmprint image and extracts the ROI using a tangent-based segmentation
%                 method
% Creation Date : 2010/01
% Revision Date : 2025/02/10
% ----------------------------------------------------------------------------------------------------

function roi = PreTangent(hand, database, figureFlag)

% Palmprint ROI Extraction Using an Enhanced Tangent-Based Preprocessing Method
%
% This function is based on the tangent-based preprocessing method introduced by 
% Dr. Zhang in the "Online Palmprint Identification" paper (2003), with significant 
% enhancements to improve compatibility across three different palmprint databases: 
% PolyU, Delhi, and CASIA.
%
% Enhancements and Modifications:
% 1. Hand Orientation Estimation:
%    - Developed a new method to determine the hand's orientation before ROI extraction.
% 2. Boundary Tracing & Gap Selection:
%    - Introduced a method to determine the optimal starting point for contour tracing.
%    - Enhanced finger gap selection criteria for improved tangent line computation.
% 3. ROI Assessment and Validation:
%    - After the ROI is extracted, the function assesses the image based on the volume
%      of the non-palm area (black area).
%
% Function Usage:
% roi = PreTangent(hand, database, figureFlag)
% - hand      : Input palmprint image.
% - database  : Specifies the database ('polyu', 'delhi', or 'casia').
% - figureFlg : Controls whether to display the processing results.
% - roi       : Extracted Region of Interest (ROI), or an empty matrix if extraction fails.
%
% Reference:
% - Zhang, D., Kong, W., You, J., & Wong, M. (2003). "Online Palmprint Identification". 
%   IEEE Transactions on Pattern Analysis and Machine Intelligence. 

% ***** FIGURE ***** ==> These lines contain the main commands for the figureFlag illustration option
% ----- FIGURE ----- ==> These lines depict parts used for debugging or% more information on the images
% Uncomment the lines containing the ----- FIGURE ----- to depict the complete stages of the ROI
% extraction, if necessary.

% ========================================================================
% ------------------------------ Database parameters
imSize = size(hand);

switch database
    case 'polyu' % --------------- Polytech database
        % Special window length for Gaussian filtering
        winLength = 11;   % Calibrated
        % The threshold for the gray to binary image conversion
        threshold = 0.09; % Calibrated
        % The coefficient that specifies the number of boundary points
        coefNo = 0.85;
    
    case 'delhi' % --------------- Delhi database
        winLength = 9;    % Calibrated
        threshold = 0.31; % Calibrated
        coefNo    = 1.1;  % Calibrated
    
    case 'casia' % --------------- CASIA database
        winLength = 7;    % Calibrated
        threshold = 0.2;  % Calibrated
        coefNo = 0.85;
end
%figure; imhist(hand); % ----- FIGURE -----

% ========================================================================

% ------------------------------ Applying the Gaussian filter
winGaus = fspecial('gaussian', winLength, 3); 
imfilt  = imfilter(hand, winGaus);
%figure; imshowpair(hand, imfilt, 'montage'); % ----- FIGURE -----

% ------------------------------ Thresholding   
imbw = imbinarize(imfilt, threshold);
%figure; imshowpair(hand, imbw, 'montage'); % ----- FIGURE -----

% ------------------------------ Eliminating the non-hand parts
palmbw = bwareafilt(imbw, 1, 'largest', 8);
%figure; imshowpair(hand, palmbw, 'montage'); % ----- FIGURE -----

% ------------------------------ Fill the holes in the palm area
palmbw = imfill(palmbw, "holes");
%figure; imshowpair(hand, palmbw, 'montage'); % ----- FIGURE -----

% ------------------------------ Edge detection 
imedge = edge(palmbw, 'canny');
if figureFlag
    figure; subplot(2,2,1); imshow(imedge); hold on; % ***** FIGURE *****
    title('Coordinate System and ROI Boundaries'); % ***** FIGURE *****
end

% ========================================================================
% ------------------------------ Obtaining the gaps between the fingers

% -------------------- Calculating the hand's area
% Finding the first non-zero left, right, top, and bottom points of the hand's binary image
% to eliminate the search area within left and right columns and top and bottom rows
[leftX, leftY]   = find(imedge == 1, 1, 'first');
%plot(leftY, leftX, 'rx', 'LineWidth', 2); % ----- FIGURE -----
%[rightX, rightY] = find(imedge == 1, 1, 'last');
%plot(rightY, rightX, 'rx', 'LineWidth', 2); % ----- FIGURE -----
[topY, topX] = find(imedge' == 1, 1, 'first'); % Utilizing the transposed image
%plot(topY, topX, 'rx', 'LineWidth', 2); % ----- FIGURE -----
[bottomY, bottomX] = find(imedge' == 1, 1, 'last'); % Utilizing the transposed image
%plot(bottomY, bottomX, 'rx', 'LineWidth', 2); % ----- FIGURE -----

% -------------------- Preliminary estimation of the hand orientation
% Calculationg the hand center point
centH = regionprops(palmbw, 'Centroid');
% The structuring element has only one element because there is only one connected component
centHand = round([centH.Centroid(2), centH.Centroid(1)]); % Hand center point
%plot(centHand(2), centHand(1), 'rx', 'LineWidth', 2); % ----- FIGURE -----

% Calculating the center point of the right side of the hand's center point (centHand)
centH = regionprops(palmbw(:,centHand(2):end), 'Centroid');
% The structuring element may have more than one element because there may be more than one
% connected component after cropping the image
% If there is more than one element, the desired center point is the closest to the hand center point
if length(centH) == 1
    centRist = round([centH.Centroid(2), centH.Centroid(1)+centHand(2)-1]); % Rist center point
else
    dist = zeros(size(centH));
    for i = 1:length(centH)
        % Using Manhattan distance (city block distance) to calculate the distance from the
        % hand center point
        dist(i) = abs(centHand(1) - centH(i).Centroid(2)) + abs(centHand(2) - centH(i).Centroid(1));
    end
    [~, id] = min(dist);
    centRist = round([centH(id).Centroid(2), centH(id).Centroid(1)+centHand(2)-1]); % Rist center point
end
%plot(centRist(2), centRist(1), 'rx', 'LineWidth', 2); % ----- FIGURE -----

% The preliminary estimation of the hand is the slope of a line that passes through these two
% points (Hand center point: centHand, rist center point: centRist)
handOri = (centHand(2)-centRist(2))/(centHand(1)-centRist(1)); % Hand's orientation
% Calculating the line that passes through these two points
% Line equation ==> Y = (Y2-Y1)/(X2-X1)*(X-X1) + Y1 = slope*(X-X1) + Y1
% Line equation ==> X = (X2-X1)/(Y2-Y1)*(Y-Y1) + X1 = 1/slope*(Y-Y1) + Y1
YhandLine = centHand(2):-1:leftY; % The columns of the line
XhandLine = 1/handOri*(YhandLine - centHand(2)) + centHand(1); % The rows of the line
%plot(YhandLine, XhandLine, 'g', 'LineWidth', 1); % ----- FIGURE -----

% -------------------- Finding the fingers except the thumb to specify the first and last
%                      points for the boundary tracing
% Computing perpendicular lines to the calculated hand's line, which also passes through the
% points on the left side of the hand's line. Then, these lines are utilized to find the
% four fingers. If the line passes through the four fingers, it must have three gaps and
% four connected components. Note that the first line that has these features is not the
% best one to use for specifying the first and last points for the boundary tracing.
% Therefore, an exhaustive search is used to find all the lines that pass through the four
% fingers. Afterward, the middle line in the series of found lines is used to specify the
% first and last points for the boundary tracing.
%
% Selecting the first and last points for boundaries
% In order to find these points, the exhaustive search is started from the closest line to the
% hand's center point. The only part of the line that includes positive column values is
% considered and mapped to the closest pixels. Then, the common part of the hand binary
% image and this mapped line is extracted. If the extracted array contains four fingers,
% there are four sets of 1s in this column. Therefore, the connected components are
% specified and labeled. If there are four connected components, then the corresponding
% index of the line is saved, including the number of the found lines so far. After the
% complete search, the middle element in the founded lines is considered to extract the
% boundaries' first and last points.

% If two lines are perpendicular, their slopes are opposite and symmetrical (slope1*slope2 = -1).
XfinLine = topX:bottomX; % The rows of the line that passes through four fingers
%YfinLine = -1/handOri*(XfinLine - centHand(1)) + centHand(2); % ----- FIGURE -----
%plot(YfinLine, XfinLine, 'g', 'LineWidth', 1); % ----- FIGURE -----

count = 0; % Number of the found lines
finLoc = []; % corresponding index of the line (Finger location)
for i = 1:length(YhandLine)
    YfinLine = round(-1/handOri*(XfinLine - XhandLine(i)) + YhandLine(i)); % The columns of the line
    id = YfinLine > 0 & YfinLine <= imSize(2); % Marking the column value withing the image size
    %plot(YfinLine(id), XfinLine(id), 'g', 'LineWidth', 1); % ----- FIGURE -----
    finID = sub2ind(imSize, XfinLine(id), YfinLine(id)); % Converting the line's X and Y to indexes
    % The common part of the hand binary image and this mapped line is extracted, and the
    % connected components and the number of them is specified
    [fingers, finNo] = bwlabeln(palmbw(finID)); % Label connected components in binary image/array
    if finNo == 4
        count = count + 1;
        finLoc(count) = i; % Saving the index of the point on the hand's line (XhandLine,YhandLine)
    end
end

% If the number of found lines is zero, return an error/warning/message and terminate;
% if not, specify the middle element 
if count == 0
    roi = [];
    %error('Please place the hand properly.')
    %warning('Please place the hand properly.')
    disp('Please place the hand properly.')
    return;
else
    id = finLoc(ceil(count/2));
end

% Calculating the line, its corresponding points, and the connected components on it
YfinLine = round(-1/handOri*(XfinLine - XhandLine(id)) + YhandLine(id));
id = YfinLine > 0;
%plot(YfinLine(id), XfinLine(id), 'g', 'LineWidth', 1); % ----- FIGURE -----
finID = sub2ind(imSize, XfinLine(id), YfinLine(id));
fingers = bwlabeln(palmbw(finID));

% The top boundary's starting point is the first point of the first connected component.
% The top boundary's last point is the last point of the second connected component.
% The bottom boundary's starting point is the first point of the fourth connected component.
% The bottom boundary's last point is the last point of the third connected component.

% The first point for the top boundary
[row, col] = ind2sub(imSize, finID(1));
id = find(fingers == 1, 1, 'last');
Fin1X = XfinLine(id+row-XfinLine(1));
Fin1Y = YfinLine(id+row-XfinLine(1));
% The last point for the top boundary
id = find(fingers == 2, 1, 'first');
Fin2X = XfinLine(id+row-XfinLine(1));
Fin2Y = YfinLine(id+row-XfinLine(1));
% The last point for the bottom boundary
id = find(fingers == 3, 1, 'last');
Fin3X = XfinLine(id+row-XfinLine(1));
Fin3Y = YfinLine(id+row-XfinLine(1));
% The first point for the bottom boundary
id = find(fingers == 4, 1, 'first');
Fin4X = XfinLine(id+row-XfinLine(1));
Fin4Y = YfinLine(id+row-XfinLine(1));
%plot(Fin1Y, Fin1X, 'yx', 'LineWidth', 2); % ----- FIGURE -----
%plot(Fin2Y, Fin2X, 'yx', 'LineWidth', 2); % ----- FIGURE -----
%plot(Fin3Y, Fin3X, 'yx', 'LineWidth', 2); % ----- FIGURE -----
%plot(Fin4Y, Fin4X, 'yx', 'LineWidth', 2); % ----- FIGURE -----

% -------------------- % Boundary tracing and gap selection
bound = cell(2,1); % Creating a cell array for the two boundaries coordinates
gap   = cell(2,1); % Creating a cell array for the two gaps coordinates

% ---------- Top boundary
% The calculated first point is not necessarily a point on the image edge. Therefore, the
% closest point on the edge from the calculated first point is considered. For this
% purpose, a small window around the first calculated point is searched for an image
% edge. If not found, the size of the window grows. The only limit is if the first
% calculated point is at the top of the image.
i = 0;
while imedge(Fin1X,Fin1Y) ~= 1
    i = i + 1;
    win1D = i; % Window size (win1D*win1D)
    if Fin1X > win1D
        if Fin1Y > win1D
            [row, col] = find(imedge(Fin1X-win1D:Fin1X+win1D, Fin1Y-win1D:Fin1Y+win1D) == 1, 1, 'first');
            if ~isempty(row)
                Fin1X = Fin1X - win1D + row - 1;
                Fin1Y = Fin1Y - win1D + col - 1;
            end
        else
            % The calculated first point is at the left side of the image
            [row, col] = find(imedge(Fin1X-win1D:Fin1X+win1D,           1:Fin1Y+win1D) == 1, 1, 'first');
            if ~isempty(row)
                Fin1X = Fin1X - win1D + row - 1;
                Fin1Y = col;
            end
        end
    else
        if Fin1Y > win1D
            % The calculated first point is at the top of the image
            [row, col] = find(imedge(          1:Fin1X+win1D, Fin1Y-win1D:Fin1Y+win1D) == 1, 1, 'first');
            if ~isempty(row)
                Fin1X = row;
                Fin1Y = Fin1Y - win1D + col - 1;
            end
        else
            % The calculated first point is at the top-left corner of the image
            [row, col] = find(imedge(          1:Fin1X+win1D,           1:Fin1Y+win1D) == 1, 1, 'first');
            if ~isempty(row)
                Fin1X = row;
                Fin1Y = col;
            end
        end
    end
end
if figureFlag, plot(Fin1Y, Fin1X, 'gx', 'LineWidth', 2); end % ***** FIGURE *****

% Calculating the number of boundary elements that are needed
% This number is approximated from the distance of the starting point and the hand center point
boundNo = round(coefNo*(abs(centHand(1)-Fin1X) + abs(centHand(2)-Fin1Y))); % Manhattan distance

% Obtaining top boundary
bound{1} = bwtraceboundary(imedge, [Fin1X, Fin1Y], 'E', 8, boundNo, 'clockwise');
%plot(bound{1}(:,2), bound{1}(:,1), 'y', 'LineWidth', 2); % ----- FIGURE -----

% Specifying top boundary's last point
% The calculated last point is not necessarily a point on the boundary. Therefore, the closest
% point on the boundary from the calculated last point is considered. 
[~, id] = min(abs(bound{1}(:,1)-Fin2X) + abs(bound{1}(:,2)-Fin2Y)); % Manhattan distance
if figureFlag, plot(bound{1}(id,2), bound{1}(id,1), 'gx', 'LineWidth', 2); end % ***** FIGURE *****

% Eliminating the extra boundary points
if id < length(bound{1})
    bound{1}(id+1:end,:) = [];
end
if figureFlag, plot(bound{1}(:,2), bound{1}(:,1), 'g', 'LineWidth', 2); end % ***** FIGURE *****

% ---------- Bottom boundary
% The calculated first point is not necessarily a point on the image edge. Therefore, the
% closest point on the edge from the calculated first point is considered. For this
% purpose, a small window around the first calculated point is searched for an image
% edge. If not found, the size of the window grows. The only limit is if the first
% calculated point is at the bottom of the image.
i = 0;
while imedge(Fin4X,Fin4Y) ~= 1
    i = i + 1;
    win1D = i; % Window size (win1D*win1D)
    if imSize(1) - Fin4X >= win1D
        if Fin4Y > win1D
            [row, col] = find(imedge(Fin4X-win1D:Fin4X+win1D, Fin4Y-win1D:Fin4Y+win1D) == 1, 1, 'first');
        else
            % The calculated first point is at the left side of the image
            [row, col] = find(imedge(Fin4X-win1D:Fin4X+win1D,           1:Fin4Y+win1D) == 1, 1, 'first');
        end
    else
        if Fin4Y > win1D
            % The calculated first point is at the bottom of the image
            [row, col] = find(imedge(Fin4X-win1D:imSize(1), Fin4Y-win1D:Fin4Y+win1D) == 1, 1, 'first');
        else
            % The calculated first point is at the bottom-left corner of the image
            [row, col] = find(imedge(Fin4X-win1D:imSize(1),           1:Fin4Y+win1D) == 1, 1, 'first');
        end
    end
    if ~isempty(row)
        Fin4X = Fin4X - win1D + row - 1;
        if Fin4Y > win1D
            Fin4Y = Fin4Y - win1D + col - 1;
        else
            Fin4Y = col;
        end
    end
end
if figureFlag, plot(Fin4Y, Fin4X, 'gx', 'LineWidth', 2); end % ***** FIGURE *****

% Calculating the number of boundary elements that are needed
% This number is approximated from the distance of the starting point and the hand center point
boundNo = round(coefNo*(abs(centHand(1)-Fin4X) + abs(centHand(2)-Fin4Y))); % Manhattan distance

% Obtaining bottom boundary
bound{2} = bwtraceboundary(imedge, [Fin4X, Fin4Y], 'E', 8, boundNo, 'counterclockwise');
%plot(bound{2}(:,2), bound{2}(:,1), 'y', 'LineWidth', 2); % ----- FIGURE -----

% Specifying bottom boundary's last point
% The calculated last point is not necessarily a point on the boundary. Therefore, the closest
% point on the boundary from the calculated last point is considered. 
[~, id] = min(abs(bound{2}(:,1)-Fin3X) + abs(bound{2}(:,2)-Fin3Y)); % Manhattan distance
if figureFlag, plot(bound{2}(id,2), bound{2}(id,1), 'gx', 'LineWidth', 2); end % ***** FIGURE *****

% Eliminating the extra boundary points
if id < length(bound{2})
    bound{2}(id+1:end,:) = [];
end
if figureFlag, plot(bound{2}(:,2), bound{2}(:,1), 'g', 'LineWidth', 2); end % ***** FIGURE *****

% -------------------- Gap selection
% Selecting an efficient part of the boundaries (gaps) for tangent calculation
FinY = [max(Fin1Y,Fin2Y), max(Fin3Y,Fin4Y)];
for i = 1:2
    Max = max(bound{i}(:,2));
    gapY = round(Max - (Max-FinY(i))/2);
    idMax1 = find(bound{i}(:,2) == Max, 1, 'first');
    idMax2 = find(bound{i}(:,2) == Max, 1, 'last');
    idGap1 = find(bound{i}(:,2) == gapY, 1, 'first');
    idGap2 = find(bound{i}(:,2) == gapY, 1, 'last');
    if ~isempty(idGap1) && ~isempty(idGap2) && idGap1 < idMax1 && idMax2 < idGap2
        gap{i} = bound{i}(idGap1:idGap2,:);
    else
        gap{i} = bound{i};
    end
    if figureFlag, plot(gap{i}(:,2), gap{i}(:,1), 'r', 'LineWidth', 2); end % ***** FIGURE *****
end

% ========================================================================
% ------------------------------ Calculating the tangent of the two gaps

% If less than two gaps are found, return an error/warning/message and terminate
if isempty(gap{1}) || isempty(gap{2})
    roi = [];
    %error('Please place the hand properly.')
    %warning('Please place the hand properly.')
    disp('Please place the hand properly.')
    return;
end

% .-----> Y coordinate / column
% |
% | 
% V
% X coordinate / row
gap1len = size(gap{1},1); % Length of first gap
gap2len = size(gap{2},1); % Length of second gap

% -------------------- Forming the most probable tangent line as the starting line
% Selecting the first point with the maximum y (column) in the first curve/gap
[Y1, id1] = max(gap{1}(:,2));
X1 = gap{1}(id1,1);
% Selecting the first point with the maximum y (column) in the second curve/gap
[Y2, id2] = max(gap{2}(:,2));
X2 = gap{2}(id2,1);

% If the line's slope is negative (Y1<Y2 ==> \), the last point with maximum y (column) in the
% second curve/gap must be considered. Also, if the line's slope is positive (Y1>Y2 ==> /),
% the last point with maximum y (column) in the first curve/gap must be considered.
if     Y1 < Y2
    [Y2, id2] = max(gap{2}(end:-1:1,2));
    id2 = gap2len-id2+1;
    X2 = gap{2}(id2,1);
elseif Y1 > Y2
    [Y1, id1] = max(gap{1}(end:-1:1,2));
    id1 = gap1len-id1+1;
    X1 = gap{1}(id1,1);
end
%plot(Y1, X1, 'rx', 'LineWidth', 2); % ----- FIGURE -----
%plot(Y2, X2, 'rx', 'LineWidth', 2); % ----- FIGURE -----
%plot([Y2 Y1], [X2 X1], 'g', 'LineWidth', 1); % ----- FIGURE -----

% -------------------- Calculating the exact tangent line of the two curves
% Check if the specified line is the actual tangent of the line. If it is the tangent line,
% do nothing. If not, choose another point on the gap1 and consider all the points on
% the second gap till the line is the actual tangent of the second gap. Again, test if
% the new line is the tangent of the first line. If yes, pass this part; if not, choose
% the next point on the first gap and repeat the previous process for the second gap.
% Continue this process till the exact tangent of gap1 and gap2 is obtained.
flag1 = 1; flag2 = 1;
while flag1
    line1 = (Y1 - Y2)/(X1 - X2)*(gap{1}(:,1) - X1) + Y1; % Sub-line 1 is the tangent of the gap1
    if gap{1}(:,2) <= line1
        flag1 = 0;
    else
        flag1 = 1;
        flag2 = 1;
        if id1 == gap1len
            id1 = 1;
        else
            id1 = id1 + 1;
        end
        X1 = gap{1}(id1,1); 
        Y1 = gap{1}(id1,2);
    end
    while flag2
        line2 = (Y1 - Y2)/(X1 - X2)*(gap{2}(:,1) - X2) + Y2; % Sub-line 2 is the tangent of the gap2
        if gap{2}(:,2) <= line2
            flag2 = 0;
        else
            flag1 = 1;
            flag2 = 1;
            if id2 == gap2len
                id2 = 1;
            else
                id2 = id2 + 1;
            end
            X2 = gap{2}(id2,1); 
            Y2 = gap{2}(id2,2);
        end
    end
end
if figureFlag, plot([Y2 Y1], [X2 X1], 'rx', 'LineWidth', 2); end % ***** FIGURE *****
tanSlope = (Y1 - Y2)/(X1 - X2); % Tangent line slope
theta = atan(tanSlope); % Tangent line angle in radians

if figureFlag
    % Calculating the line points that cross these two points
    X = [1, imSize(1)]; % ***** FIGURE *****
    Y = tanSlope*(X - X1) + Y1; % Calculating the Y value for the line % ***** FIGURE *****
    plot(Y, X, 'y', 'LineWidth', 1); % ***** FIGURE *****
end

% ========================================================================
% ------------------------------ Specifying the Region Of Interest (ROI)

% -------------------- Database parameters
switch database
    case 'polyu' % --------------- Polytech database
        % Size of the ROI minus one
        roiSize = [127 127]; % Calibrated
        % The X (row) coordinate distance of the top left point in the cropped square from the
        % tangent line midpoint after rotation.
        midRow(1) = ceil(roiSize(1)/2);
        % The X (row) coordinate distance of the bottom left point in the cropped square from the
        % tangent line midpoint after rotation.
        midRow(2) = roiSize(1) - midRow(1);
        % The Y (column) coordinate distance of the top left point in the cropped square from the
        % tangent line midpoint after rotation.
        midCol = 50; % Calibrated (PreTangent)
    
    case 'delhi' % --------------- Delhi database
        roiSize = [209 209];
        midRow(1) = ceil(roiSize(1)/2) - 10 ;
        midRow(2) = roiSize(1) - midRow(1);
        midCol = 30;
    
    case 'casia' % --------------- CASIA database
        roiSize = [255 255];
        midRow(1) = ceil(roiSize(1)/2);
        midRow(2) = roiSize(1) - midRow(1);
        midCol = 50;
end

% -------------------- Calculating the line midpoint coordinates in the original image
% Calculating the line midpoint
midLine = mean([X1 Y1 ; X2 Y2], 1);

if figureFlag
    plot(midLine(2), midLine(1), 'yx', 'LineWidth', 2); % ***** FIGURE *****
    
    % If two lines are perpendicular, their slopes are opposite and symmetrical (slope1*slope2 = -1).
    % Original ine equations       ==> Y = SLOPE*(X - X0) + Y0, X = 1/SLOPE*(Y - Y0) + X0
    % Perpendicular line equations ==> Y90 = -1/SLOPE*(X90 - X0) + Y0, X90 = -SLOPE*(Y90 - Y0) + X0
    Y90 = [1, imSize(2)]; % ***** FIGURE *****
    X90 = -tanSlope*(Y90 - midLine(2)) + midLine(1); % Perpendicular line % ***** FIGURE *****
    plot(Y90, X90, 'y', 'LineWidth', 1); % ***** FIGURE *****
    
    % -------------------- Depicting the coordinate axis and the corners' coordinates of the cropped ROI
    % .-----> Y coordinate / column
    % |
    % | 
    % V
    % X coordinate / row
    
    % Calculating X and Y of a point with a fixed distance from another point on a line
    % X, Y = desired point's coordinates, X0, Y0 = another point's coordinates, fixed distance = dist
    % Y = a*X + b = a*(X - X0) + Y0, theta = atan(a) ==> theta is in radians
    % First approach
    % X = +-cos(theta)*dist + X0, Y = +-sin(theta)*dist + Y0, Y = tan(theta)*X + b = a*X + b
    % Second approach
    % dist^2 = (X - X0)^2 + (Y - Y0)^2, Y = a*X + b ==> 
    % dist^2 = (X - X0)^2 + (aX + b - aX0 - b)^2     ==> 
    % dist^2 = (X - X0)^2 + a^2*(X - X0)^2           ==>
    % dist^2 = (1 + a^2)*(X - X0)^2                  ==>
    % dist = sqrt(1 + a^2)*(X - X0)                  ==>
    % (X - X0) = dist/sqrt(1 + a^2)                  ==>
    % X = +-dist/sqrt(1 + a^2) + X0, Y = +-a/sqrt(1 + a^2)*dist + Y0, Y = a*X + b
    
    corners = zeros(4,2); % ***** FIGURE *****
    for i = -1:2:1 % ***** FIGURE *****
        % (Xdist, Ydist): Coordinates of the two points on the tangent line with a fixed
        % distance from the midpoint.
        % "i" is used to specify the sign of the cosine function. To use the "i" for midRow indexes,
        % -1 must be mapped to 1 and 1 to 3. The equation to do so is: ind = i/2 + 1.5
        Xdist = i*cos(theta)*midRow(i/2+1.5) + midLine(1); % cos(theta) = cos(-theta) % ***** FIGURE *****
        Ydist = tanSlope*(Xdist - midLine(1)) + midLine(2); % ***** FIGURE *****
        %plot(Ydist, Xdist, 'yx', 'LineWidth', 2); % ----- FIGURE -----
        
        % Plot the perpendicular line to the tangent line that passes through this point
        %Y90 = [1, imSize(2)]; % ----- FIGURE -----
        %X90 = -tanSlope*(Y90 - Ydist) + Xdist; % Perpendicular line % ----- FIGURE -----
        %plot(Y90, X90, 'y', 'LineWidth', 1); % ----- FIGURE -----
        
        % Calculating two corners
        % sin(pi/2 + theta) = cos(theta)
        % Y = +-sin(pi/2+theta)*dist + Y0 = +-cos(theta)*dist + Y0
        j = i + 2; % ***** FIGURE *****
        % Corners on the left side
        corners(j,2) = abs(cos(theta))*midCol + Ydist; % ***** FIGURE *****
        corners(j,1) = -tanSlope*(corners(j,2) - Ydist) + Xdist; % ***** FIGURE *****
        %plot(corners(j,2), corners(j,1), 'yx', 'LineWidth', 2); % ----- FIGURE -----
        % Corners on the right side
        corners(j+1,2) = abs(cos(theta))*(midCol+roiSize(2)) + Ydist; % ***** FIGURE *****
        corners(j+1,1) = -tanSlope*(corners(j+1,2) - Ydist) + Xdist; % ***** FIGURE *****
        %plot(corners(j+1,2), corners(j+1,1), 'yx', 'LineWidth', 2); % ----- FIGURE -----
    end % ***** FIGURE *****
    plot([corners(1,2) corners(2,2)], [corners(1,1) corners(2,1)], 'g', 'LineWidth', 1); % ***** FIGURE *****
    plot([corners(1,2) corners(3,2)], [corners(1,1) corners(3,1)], 'g', 'LineWidth', 1); % ***** FIGURE *****
    plot([corners(4,2) corners(2,2)], [corners(4,1) corners(2,1)], 'g', 'LineWidth', 1); % ***** FIGURE *****
    plot([corners(4,2) corners(3,2)], [corners(4,1) corners(3,1)], 'g', 'LineWidth', 1); % ***** FIGURE *****
    hold off; % ***** FIGURE *****
end

% -------------------- Calculating the coordinates of the midpoint after rotating the image
% To crop the ROI image, firstly, the image is rotated using the "imrotate" function, and
% then a square or rectangle window is cropped from the original image. For this purpose,
% the coordinates of the tangent line's midpoint are needed in the rotated image.

% Calculating the tangent line rotation angle in degrees for the "imrotate" function 
alpha = -theta*180/pi; % alpha = -atand(tanSlope);

% Calculating the location of the line midpoint in the rotated image
midLine = round(midLine);
% Construct a new image by the size of the hand image. Then, mark the line midpoint in it.
% Afterward, rotate it and find the mark in the new image. Note that the mark starts from
% one pixel, and if the rotated image only contains zero pixels, the mark area gets expanded.
midMark = zeros(imSize); % Constructing the new image
markCord = []; % Line midpoint's new coordinates in the rotated image
i = 0;
while isempty(markCord)
    % Mark the line's midpoint area in the new image
    midMark(midLine(1)-i:midLine(1)+i, midLine(2)-i:midLine(2)+i) = 255;
    i = i + 1;
    midRotate = imrotate(midMark, alpha, 'nearest', 'crop'); % Rotating the new image
    [row, col] = find(midRotate == 255); % Finding a set of coordinates
    markCord = [row, col];
end
markCord = round(mean(markCord, 1)); % Calculating the center of the found marks

if figureFlag
    edgeRotate = imrotate(imedge, alpha, 'bicubic', 'crop'); % ***** FIGURE *****
    subplot(2,2,2); imshow(edgeRotate); hold on; % ***** FIGURE *****
    title('Rotated Coordinate System and ROI Boundaries'); % ***** FIGURE *****
    plot(markCord(2), markCord(1), 'yx', 'LineWidth', 2); % ***** FIGURE *****
    plot([markCord(2), markCord(2)], [1, imSize(1)], 'y', 'LineWidth', 1); % ***** FIGURE *****
    plot([1, imSize(2)], [markCord(1), markCord(1)], 'y', 'LineWidth', 1); % ***** FIGURE *****
end

% -------------------- Specifying the top left corner of the ROI and depicting the desired area
topLeft = [markCord(1)-midRow(1), markCord(2)+midCol]; % Top left corner of the cropped image
if figureFlag
    rectangle('Position', [topLeft(2) topLeft(1) roiSize(2) roiSize(1)], ... % ***** FIGURE *****
              'EdgeColor', 'g', 'LineWidth', 1); hold off; % ***** FIGURE *****
end

% -------------------- Extracting the ROI
handRotate = imrotate(hand, alpha, 'bicubic', 'crop'); % Rotationg the image
roi = imcrop(handRotate, [topLeft(2) topLeft(1) roiSize(2) roiSize(1)]); % Cropping the desired area
if figureFlag
    subplot(2,2,3); imshow(hand); title('Hand Image'); % ***** FIGURE *****
    subplot(2,2,4); imshow(roi);  title('Segmented ROI'); % ***** FIGURE *****
end

% -------------------- Assessing the cropped image
% After cropping the RIO image, the region of interest is assessed to see if it contains the
% desired part of the palm and is a suitable candidate for the learning or identification
% processes. Here, the quality of the ROI is assessed with the occupied area of the image.
% For this purpose, the ROI image is binarized, and the sum of the non-zero pixels is
% calculated. Then, the ratio of the ROI part against the whole image is computed. If it
% is less than a certain value, the function returns an empty image, which will also not
% be saved. Also, a message will be displayed to notify the user that the hand image was
% inappropriate and needs to be taken again.
switch database
    case 'polyu' % --------------- Polytech database
        maxRatio = 0.93; % Calibrated
    case 'delhi' % --------------- Delhi database
        maxRatio = 0.8;  % Calibrated
    case 'casia' % --------------- CASIA database
        maxRatio = 0.83; % Calibrated
end
roiRegion = sum(imbinarize(roi, threshold), 'all')/(size(roi,1)*size(roi,2));
if roiRegion < maxRatio
    roi = [];
    disp('Please place the hand properly.')
    return;
end

% -------------------- Resizing the cropped image
% In the case that the cropped image size is smaller than the desired ROI size, the final image
% must be expanded. For this purpose, if the columns of the cropped image are less than the
% specified value "roiSize(2)", columns of zeros will be added to the right side of the
% image. If the rows of the cropped image are less than the determined value "roiSize(1)",
% corresponding rows will be added to the top and bottom of the image to maintain the ROI
% in the center of the image's rows. The corresponding rows are calculated as follows. The
% difference between the desired and cropped image rows will be computed. Then, half of it
% will be added to the top of the image and the rest to the bottom. Therefore, the final
% ROI image stays at the center of the rows and the left side of the adjusted image. Note
% that the "roiSize" matrix's values are one point less than the desired ROI size.
outIMsize = size(roi);
if outIMsize(2) <= roiSize(2)
    diff = roiSize(2) - outIMsize(2) + 1;
    roi = [roi, zeros(outIMsize(1), diff)];
end
if outIMsize(1) <= roiSize(1)
    diff = roiSize(1) - outIMsize(1) + 1;
    roi = [zeros(ceil(diff/2), roiSize(2)+1) ; roi ; zeros(diff - ceil(diff/2), roiSize(2)+1)];
end

end