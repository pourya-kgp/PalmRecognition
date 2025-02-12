% ----------------------------------------------------------------------------------------------------
% Author        : Pourya Khodagholipour (P.KH)
% Project Name  : Implementation of a Palmprint Recognition System based on Image Processing (M.S. Thesis)
% File Name     : PreMorphological.m
% Description   : Processes a palmprint image and extracts the ROI using a morphological-based 
%                 segmentation method
% Creation Date : 2010/02
% Revision Date : 2025/02/10
% ----------------------------------------------------------------------------------------------------

function  roi = PreMorphological(hand, database, figureFlag)

% Palmprint ROI Extraction Using a Morphological-Based Method
%
% This function extracts the Region of Interest (ROI) from palmprint images using a novel 
% morphological-based preprocessing approach.
%
% Function Usage:
% roi = PreMorphological(hand, database, figureFlag)
% - hand       : Input palmprint image.
% - database   : Specifies the database ('polyu', 'delhi', or 'casia').
% - figureFlag : Enables visualization of the finger detection process & processing results.
% - roi        : Extracted Region of Interest (ROI), or an empty matrix if extraction fails.
%
% Reference:
% Morphological-based segmentation method:
%   - Khodagholipour, P., Hamed S. (2010). "A Novel Method for Locating Region-Of-Interest
%     in Palmprint Recognition Systems"
%     IEEE 2010 International Conference on Progress in Informatics and Computing (PIC),
%     Shanghai, China, (Accepted, Not Registered).
%     It can be found at the address:
%     https://github.com/pourya-kgp/PalmRecognition/blob/main/Docs/2010_IEEE_PIC_ROI_Palmprint.pdf

% ========================================================================
% ------------------------------ Visualization option

% First flag (figureFlag1) controls whether to display a figure showing the stages of
% finger detection. In this figure, the process of finding the fingers, as well as the
% detection of the index and little fingers, is presented in multiple subplots.  
%
% Second flag (figureFlag2) controls whether to display a figure showing the processing
% results. In this figure, the hand image, coordinate system with the ROI boundaries,
% rotated fingers with the ROI boundaries, and the segmented ROI are presented.  
%
% ===== FIGURE ===== ==> These lines contain the main commands for the finger detection process
%                        illustration. The associated flag to them is figureFlag1.
% ***** FIGURE ***** ==> These lines contain the main commands for the processing results
%                        illustration. The associated flag to them is figureFlag2.
% ----- FIGURE ----- ==> These lines depict parts used for debugging or more information on the image
% Uncomment the lines containing the ----- FIGURE ----- to depict the complete stages of the ROI
% extraction, if necessary.

if figureFlag
    figureFlag1 = true; % Options: true, false (Finger detection process flag)
    figureFlag2 = true; % Options: true, false (Processing results flag)
else
    figureFlag1 = false;
    figureFlag2 = false;
end

% ========================================================================
% ------------------------------ Database parameters

imSize = size(hand);

switch database
    case 'polyu' % --------------- Polytech database
        % Number of fingers in the image
        fingNo = 4;
        % The thresholds for the gray to binary image conversion
        threshold1 = 0.01;           % Calibrated
        threshold2 = 0.09;           % Calibrated
        % Special window length for Gaussian filtering
        winLength = 11;
        % structuring elements
        SE1 = strel('disk', 3, 8);  % Calibrated
        SE2 = strel('disk', 36);     % Calibrated
        SE3 = strel('disk', 39);     % Calibrated
    
    case 'delhi' % --------------- Delhi database
        fingNo = 5;
        threshold1 = 0.31;           % Calibrated
        threshold2 = 0.31;           % Calibrated
        winLength = 9; 
        SE1 = strel('disk', 11, 8);
        SE2 = strel('disk', 50);     % Calibrated
        SE3 = strel('disk', 55);     % Calibrated
    
    case 'casia' % --------------- CASIA database
        fingNo = 5;
        threshold1 = 0.2;            % Calibrated
        threshold2 = 0.2;            % Calibrated
        winLength = 7;
        SE1 = strel('disk', 9, 8);
        SE2 = strel('disk', 74);     % Calibrated
        SE3 = strel('disk', 79);     % Calibrated
end
%figure; imhist(hand); % ----- FIGURE -----

% ========================================================================
% ------------------------------ Background Noise Reduction

% -------------------- Thresholding   
imbw = imbinarize(hand, threshold1);
%figure; imshowpair(hand, imbw, 'montage'); % ----- FIGURE -----

% -------------------- Disconnecting the non-hand parts
imopn = imopen(imbw, SE1);
%figure; imshowpair(imbw, imopn, 'montage'); % ----- FIGURE -----

% -------------------- Eliminating the non-hand parts
mask = bwareafilt(imopn, 1, 'largest', 8);
%figure; imshowpair(imbw, mask, 'montage'); % ----- FIGURE -----

% -------------------- Constructing the final hand area
palm = uint8(mask).*hand;     
%figure; imshowpair(hand, palm, 'montage'); % ----- FIGURE -----

%figure; subplot(2,2,1); imshow(imbw); subplot(2,2,2); imshow(imopn); % ----- FIGURE -----
%        subplot(2,2,3); imshow(mask); subplot(2,2,4); imshow(palm);  % ----- FIGURE -----

% ========================================================================
% ------------------------------ Hand binarized image

% -------------------- Applying the Gaussian filter
winGaus = fspecial('gaussian', winLength, 3); 
imfilt  = imfilter(palm, winGaus);
%figure; imshowpair(palm, imfilt, 'montage'); % ----- FIGURE -----

% -------------------- Thresholding   
palmbw = imbinarize(imfilt, threshold2);
%figure; imshowpair(palm, palmbw, 'montage'); % ----- FIGURE -----

% -------------------- Filling the holes in the palm area
palmbw = imfill(palmbw, "holes");
%figure; imshowpair(palm, palmbw, 'montage'); % ----- FIGURE -----
  
% ========================================================================      
% ------------------------------ Specifying fingers: index, middle, ring, little, thumb
    
% -------------------- Separating the palm from the hand image 
% Specifying the palm area
imErod = imerode(palmbw, SE2);
imErod = bwareafilt(imErod, 1, 'largest', 8); 

% Calculationg the palm center point
centH = regionprops(imErod, 'Centroid');
% The structuring element has only one element because there is only one connected component
centPalm = round([centH.Centroid(2), centH.Centroid(1)]); % Palm center point

% -------------------- Removing the palm area from the hand image
imDilat = imdilate(imErod, SE3);
nonPalm = palmbw & ~imDilat;

% -------------------- Selecting fingers from the non-palm image
% Labeling connected components in the binary image
[conComp, ccNo] = bwlabel(nonPalm);
% Computing the number of pixels in connected components (connected components volume)
ccVol = zeros(ccNo, 1);
for i = 1:ccNo
    ccVol(i) = sum(conComp==i, 'all');
end
% Selecting 4 or 5 largest connected components in the non-hand image. In other words, eliminating
% any connected components that are not fingers, which are the largest connected components.
fingers = nonPalm;
if ccNo > fingNo % If the number of connected components is bigger than the predefined finger number
    [~, id] = sort(ccVol, 'descend');
    for i = fingNo+1:length(id)
        fingers(conComp==id(i)) = 0;
    end
elseif ccNo == 4
    fingNo = 4;
elseif ccNo < 4
    % Comment to view all the results, even the corrupted results
    roi = [];
    return;
end

if figureFlag1
    figure; % ===== FIGURE =====
    subplot(2,3,1); imshow(palmbw);  title('1. Binarized Hand Image'); % ===== FIGURE =====
    subplot(2,3,2); imshow(imErod);  title('2. Palm Central Area (Eroded Image)'); % ===== FIGURE =====
    subplot(2,3,3); imshow(imDilat); title('3. Palm Area (Dilated Image)'); % ===== FIGURE =====
    subplot(2,3,4); imshow(nonPalm); title('4. Nonpalm Area'); % ===== FIGURE =====
    subplot(2,3,5); imshow(fingers); title('5. Fingers'); hold on; % ===== FIGURE =====
end

% ========================================================================      
% ------------------------------ Specifying four fingers: index, middle, ring, and little

% Calculating fingers' center points
centF = regionprops(fingers, 'Centroid'); % a structure with "fingNo" elements

% Constructing a finger's related matrix (centFinger), in which every row in it has the following elements:
% [center's row, center's column, corresponding connected component label, corresponding connected component volume]
centFinger = zeros(fingNo, 4);
for i = 1:fingNo
    centFinger(i,1:2) = round([centF(i).Centroid(2), centF(i).Centroid(1)]); % Fingers' centers
    j = 0;
    while centFinger(i,3) == 0
        window = conComp(centFinger(i,1)-j:centFinger(i,1)+j, centFinger(i,2)-j:centFinger(i,2)+j);
        centFinger(i,3) = max(window, [], 'all'); % Finger's corresponding connected component label
        j = j + 1;
    end
    centFinger(i,4) = ccVol(centFinger(i,3)); % Finger's corresponding connected component volume
end
if figureFlag1, plot(centFinger(:,2), centFinger(:,1), 'rx', 'LineWidth', 2); end % ===== FIGURE =====

if fingNo > 4
    % -------------------- Calculating the mean value of the fingers' centers
    meanCent = mean(centFinger(:,1:2), 1);
    if figureFlag1, plot(meanCent(2), meanCent(1), 'bx', 'LineWidth', 2); end % ===== FIGURE =====

    % -------------------- Calculationg the mean center of the index, middle, and ring fingers
    % Computing the Euclidean distance of the fingers' centers to the mean value of them
    dist = sqrt(sum(bsxfun(@minus, centFinger(:,1:2), [meanCent(1) meanCent(2)]).^2, 2));
    % Choosing three fingers whose centers are the nearest to the "meanCent" value. It means
    % selecting three fingers out of all fingers, which are the index, middle, and ring fingers.
    [~, id] = sort(dist, 'ascend');
    % Calculating a new mean value for the fingers' center point, which is the mean value of the
    % index, middle, and ring fingers' centers.
    meanCent = mean(centFinger(id(1:3),1:2), 1);
    if figureFlag1, plot(meanCent(2), meanCent(1), 'gx', 'LineWidth', 2); hold off; end % ===== FIGURE =====

    % -------------------- Removing the thumb
    % Computing the Euclidean distance of the fingers' centers to the new mean value of them
    dist = sqrt(sum(bsxfun(@minus, centFinger(:,1:2), [meanCent(1) meanCent(2)]).^2, 2));
    % Choosing four fingers whose centers are the nearest to the new "meanCent" value. It means
    % selecting four fingers out of all fingers, which are the index, middle, ring, and
    % little fingers. In this case, 
    [~, id] = sort(dist, 'ascend');
    fingers(conComp==centFinger(id(end),3)) = 0; % Eliminating the thumb finger
    centFinger(id(end),:) = []; % Eliminating the center of the thumb finger in the cetfinger matrix
end

if figureFlag1
    subplot(2,3,6); imshow(fingers); hold on; % ===== FIGURE =====
    title('6. Index & Little Fingers'); % ===== FIGURE =====
    plot(centPalm(2), centPalm(1), 'rx', 'LineWidth', 2); % ===== FIGURE =====
    plot(centFinger(:,2), centFinger(:,1), 'rx', 'LineWidth', 2); % ===== FIGURE =====
end

% To view the results of the four-finger segmentation (Index, middle, ring, and little fingers)
%roi = fingers;
%return;

% ========================================================================      
% ------------------------------ Rearranging the fingers' centers matrix

% -------------------- Specifying the index and little fingers
% In order to specify the index and little fingers, the cosine theorem in triangles is used.
% The angle between the lines from the centers of the index and little fingers to the palm
% center is the biggest angle of all others.
% Calculating every two possible combinations of fingers' centers
combCent = nchoosek(1:length(centFinger), 2);
combleng = length(combCent);
% Computing angles between the palm center and every two possible combinations of fingers' centers
Gama = 180*ones(combleng, 1); % Angles in degrees
for i = 1:combleng
    triSide = pdist([centPalm ; centFinger(combCent(i,1),1:2) ; centFinger(combCent(i,2),1:2)]); % Triangle sides length
    Gama(i) = acosd((triSide(1)^2 + triSide(2)^2 - triSide(3)^2)/(2*triSide(1)*triSide(2))); % Cosine theorem in triangle
end
% Specifying the biggest angle between two fingers' centers and the palm center
[~, id] = max(Gama); 

% -------------------- Rearranging the fingers' centers matrix
% Rearrange the fingers' center matrix (centFinger) so that the index and little fingers will
% be at the first or last rows of the matrix
middRingID = 1:4;
middRingID(combCent(id,:)) = [];
centFinger = centFinger([combCent(id,1), middRingID, combCent(id,2)], :);
% Rearrange the fingers' center matrix (centFinger) so that the new arrangement will be
% Index-Middle-Ring-Little or Little-Ring-Middle-Index. For this purpose, the Euclidean
% distance of the first center from the second and third centers is used. The nearest
% one is the closest finger to the first center.
dist = sqrt(sum(bsxfun(@minus, centFinger(1,1:2), centFinger(2:3,1:2)).^2, 2));
[~, id] = sort(dist, 'ascend');
centFinger = centFinger([1, id'+1, 4], :);

if figureFlag1
    plot([centFinger(1,2) centPalm(2)], [centFinger(1,1) centPalm(1)], 'g', 'LineWidth', 1); % ===== FIGURE =====
    plot([centFinger(4,2) centPalm(2)], [centFinger(4,1) centPalm(1)], 'g', 'LineWidth', 1); % ===== FIGURE =====
    hold off; % ===== FIGURE =====
end

% ========================================================================      
% ------------------------------ Specifying two key points for rotation

% -------------------- Specifying key points
% The nearest pixels of the middle and ring fingers to the gravity center of the palm are
% used to find the roots of the index and little fingers located in the gaps between the
% index-middle fingers and ring-little fingers, respectively. The nearest pixel of the
% index finger from the obtained middle finger point is considered one of the key points,
% and the same method is used to specify the other key point.
keyPoint = zeros(2);
for i = 1:2
    [row1, col1] = find(conComp==centFinger(i^2, 3)); % Index  | Little fingers
    [row2, col2] = find(conComp==centFinger(i+1, 3)); % Middle | Ring   fingers
    dist = sqrt(sum(bsxfun(@minus, [row2 col2], centPalm).^2, 2));
    [~, id] = min(dist);
    middRingID = [row2(id), col2(id)]; % The nearest Middle/Ring finger element to the palm center
    dist = sqrt(sum(bsxfun(@minus, [row1 col1], middRingID).^2, 2));
    [~, id] = min(dist);
    keyPoint(i,:) = [row1(id), col1(id)]; % The nearest Index/Little finger element to middRingID
end

% For compatibility with the "PreTangent" script, which can be found in the PalmIdentify
% repository at https://github.com/pourya-kgp/PalmIdentify
X1 = keyPoint(1,1);
Y1 = keyPoint(1,2);
X2 = keyPoint(2,1);
Y2 = keyPoint(2,2);

% Note: Most of the rest of the code is duplicated from the "PreTanget"
% script which can be find at https://github.com/pourya-kgp/PalmIdentify

tanSlope = (Y1 - Y2)/(X1 - X2); % Line slope
theta = atan(tanSlope); % Line angle in radians

if figureFlag2
    figure; subplot(2,2,1); imshow(fingers); hold on; % ***** FIGURE *****
    title('Coordinate System & ROI Boundaries'); % ***** FIGURE *****
    plot([Y2 Y1], [X2 X1], 'rx', 'LineWidth', 2); % ***** FIGURE *****
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
        % calculated line midpoint after rotation.
        midRow(1) = ceil(roiSize(1)/2);
        % The X (row) coordinate distance of the bottom left point in the cropped square from the
        % calculated line midpoint after rotation.
        midRow(2) = roiSize(1) - midRow(1);
        % The Y (column) coordinate distance of the top left point in the cropped square from the
        % calculated line midpoint after rotation.
        midCol = 55; % Calibrated (PreMorphological)
    
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

if figureFlag2
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
        % (Xdist, Ydist): Coordinates of the two points on the calculated line with a fixed distance
        % from the midpoint.
        % "i" is used to specify the sign of the cosine function. To use the "i" for midRow indexes,
        % -1 must be mapped to 1 and 1 to 3. The equation to do so is: ind = i/2 + 1.5
        Xdist = i*cos(theta)*midRow(i/2+1.5) + midLine(1); % cos(theta) = cos(-theta) % ***** FIGURE *****
        Ydist = tanSlope*(Xdist - midLine(1)) + midLine(2); % ***** FIGURE *****
        %plot(Ydist, Xdist, 'yx', 'LineWidth', 2); % ----- FIGURE -----
        
        % Plot the perpendicular line to the calculated line that passes through this point
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
% the coordinates of the calculated line's midpoint are needed in the rotated image.

% Calculating the calculated line rotation angle in degrees for the "imrotate" function 
alpha = -theta*180/pi; % alpha = -atand(tanSlope);

% Calculating the location of the line midpoint in the rotated image
midLine = round(midLine);
% Construct a new image by the size of the hand image. Then, mark the line midpoint in it.
% Afterward, rotate it and find the mark in the new image. Note that the mark starts from
% one pixel, and if the rotated image only contains zero pixels, the mark area gets expanded.
midMark = zeros(size(hand)); % Constructing the new image
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

if figureFlag2
    edgeRotate = imrotate(fingers, alpha, 'bicubic', 'crop'); % ***** FIGURE *****
    subplot(2,2,2); imshow(edgeRotate); hold on; % ***** FIGURE *****
    title('Rotated Coordinate System and ROI Boundaries'); % ***** FIGURE *****
    plot(markCord(2), markCord(1), 'yx', 'LineWidth', 2); % ***** FIGURE *****
    plot([markCord(2), markCord(2)], [1, imSize(1)], 'y', 'LineWidth', 1); % ***** FIGURE *****
    plot([1, imSize(2)], [markCord(1), markCord(1)], 'y', 'LineWidth', 1); % ***** FIGURE *****
end

% -------------------- Specifying the top left corner of the ROI and depicting the desired area
topLeft = [markCord(1)-midRow(1), markCord(2)+midCol]; % Top left corner of the cropped image
if figureFlag2
    rectangle('Position', [topLeft(2) topLeft(1) roiSize(2) roiSize(1)], ... % ***** FIGURE *****
              'EdgeColor', 'g', 'LineWidth', 1); hold off; % ***** FIGURE *****
end

% -------------------- Extracting the ROI
handRotate = imrotate(hand, alpha, 'bicubic', 'crop'); % Rotationg the image
roi = imcrop(handRotate, [topLeft(2) topLeft(1) roiSize(2) roiSize(1)]); % Cropping the desired area
if figureFlag2
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
roiRegion = sum(imbinarize(roi, threshold2), 'all')/(size(roi,1)*size(roi,2));
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