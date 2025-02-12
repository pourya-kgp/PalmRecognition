% ----------------------------------------------------------------------------------------------------
% Author        : Pourya Khodagholipour (P.KH)
% Project Name  : Implementation of a Palmprint Recognition System based on Image Processing (M.S. Thesis)
% File Name     : FeatureExtraction.m
% Description   : Extracts a set of discriminative features from a palmprint ROI image
% Creation Date : 2010/01
% Revision Date : 2025/02/10
% ----------------------------------------------------------------------------------------------------

function [reGaborBin, imGaborBin, mask] = FeatureExtraction(roi, vectSize, database, figureFlag)

% This function performs feature extraction on a region of interest (ROI) using a circular Gabor filter.
% It outputs three binary matrices representing the processed features: the real part, the imaginary 
% part, and a mask. The parameters of the Gabor filter are customized for specific databases.
% 
% Inputs:
% - roi        : Region of interest (grayscale image)
% - vectSize   : Desired vector size for feature matrices
% - database   : Name of the database ("polyu", "delhi", or "casia")
%
% Outputs:
% - reGaborBin : Binarized and resized Gabor filter output (real part)
% - imGaborBin : Binarized and resized Gabor filter output (imaginary part)
% - mask       : Binarized and resized mask of the ROI
%
% References:
% Morphological-based segmentation method:
%   - Khodagholipour, P., Hamed S. (2010). "A Novel Method for Locating Region-Of-Interest
%     in Palmprint Recognition Systems"
%     IEEE 2010 International Conference on Progress in Informatics and Computing (PIC),
%     Shanghai, China, (Accepted, Not Registered).
%     It can be found at the address:
%     https://github.com/pourya-kgp/PalmRecognition/blob/main/Docs/2010_IEEE_PIC_ROI_Palmprint.pdf
%
% Tangent-based segmentation, feature extraction, and matching methods:
%   - Zhang, D., Kong, W. K, You, J., & Wong, M. (2003). "Online Palmprint Identification". 
%     IEEE Transactions on Pattern Analysis and Machine Intelligence.

% -------------------- Database parameters
switch database
    % Gabor filter parameters:
    % theta : Orientation of the function (in radians)
    % u     : Frequency of the sinusoidal wave
    % sigma : Standard deviation of the Gaussian envelope
    % filSz : Size of the filter (It is recommended to select an odd size for the filter)
    %
    % Thresholding parameters:
    % threshold: The threshold for the gray to binary image conversion
    % sensValue = Sensitivity value for adaptive thresholding
    
    case 'polyu' % --------------- Polytech
        theta = pi/4;     % According to the paper
        u     = 0.0916;   % According to the paper
        sigma = 5.6179;   % According to the paper
        filSz = 7;        % Calibrated
        threshold = 0.09; % Calibrated
        sensValue = 0.6;  % Calibrated
    
    case 'delhi' % --------------- Delhi
        theta = pi/4;     % According to the paper
        u     = 0.0916;   % According to the paper
        sigma = 5.6179;   % According to the paper
        filSz = 13;       % Calibrated
        threshold = 0.31; % Calibrated
        sensValue = 0.6;  % Calibrated
    
    case 'casia' % --------------- CASIA
        theta = pi/4;     % According to the paper
        u     = 0.0916;   % According to the paper
        sigma = 5.6179;   % According to the paper
        filSz = 9;        % Calibrated
        threshold = 0.2;  % Calibrated
        sensValue = 0.6;  % Calibrated
end
% Note: In the PolyU and CASIA databases, it is a better approach to have a filter bank with theta 
% equal to pi/4 & -pi/4. For the right hand with lines from the top left to the bottom right, the
% theta must be pi/4. Also, for the left hand with lines from the bottom left to the top right,
% the theta must be -pi/4. In the Delhi database, all the palmprints were rotated so that they
% had approximately the same lines orientation from the top left to the bottom right.

% -------------------- Feature Extraction
% Apply the Circular Gabor Filter to the ROI
[reGabor, imGabor] = CircularGaborFilter(roi, filSz, theta, u, sigma);

% -------------------- Generating Binary Feature Matrices
% Normalize and binarize the real part of the Gabor filter output
reGabor = mat2gray(reGabor);
reGaborBin = imbinarize(reGabor, "adaptive", ForegroundPolarity="dark", Sensitivity=sensValue);

% Normalize and binarize the imaginary part of the Gabor filter output
imGabor = mat2gray(imGabor);
imGaborBin = imbinarize(imGabor, "adaptive", ForegroundPolarity="dark", Sensitivity=sensValue);

% Create a mask by binarizing the ROI with a fixed threshold
imBin = imbinarize(roi, threshold);

% -------------------- Resizing the Feature Matrices
% Ensure the vector size is even
if mod(vectSize, 2) ~= 0
    vectSize = vectSize + 1;
end

% Calculate the size for resizing the feature matrices
gabSize = round(sqrt(vectSize / 2));

% Resize the feature matrices and mask
reGaborBin = imresize(reGaborBin, [gabSize, gabSize], 'bilinear');
imGaborBin = imresize(imGaborBin, [gabSize, gabSize], 'bilinear');
mask       = imresize(imBin,      [gabSize, gabSize], 'bilinear');

% -------------------- Visualization
if figureFlag
    figure;
    subplot(2,3,1); imshow(roi);        title('ROI');
    subplot(2,3,2); imshow(reGabor);    title('Gabor Filter (Real Part)');
    subplot(2,3,3); imshow(imGabor);    title('Gabor Filter (Imaginary Part)');
    subplot(2,3,4); imshow(mask);       title('Mask');
    subplot(2,3,5); imshow(reGaborBin); title('Binarized & Resized (Real Part)');
    subplot(2,3,6); imshow(imGaborBin); title('Binarized & Resized (Imaginary Part)');
end

end