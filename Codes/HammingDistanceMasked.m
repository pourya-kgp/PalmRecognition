% ----------------------------------------------------------------------------------------------------
% Author        : Pourya Khodagholipour (P.KH)
% Project Name  : Implementation of a Palmprint Recognition System based on Image Processing (M.S. Thesis)
% File Name     : HammingDistanceMasked.m
% Description   : Computes the minimum Hamming distance between two sets of binary feature matrices,
%                 accounting for mask regions
% Creation Date : 2010/01
% Revision Date : 2025/02/10
% ----------------------------------------------------------------------------------------------------

function minHamDist = HammingDistanceMasked(reMatBin1, imMatBin1, mask1, ...
                                            reMatBin2, imMatBin2, mask2)

% HammingDistanceMasked computes the minimum Hamming distance between two sets of binary 
% feature matrices. The calculation considers potential translations of one matrix relative
% to the other within a range of -2 to 2 pixels.
% 
% Inputs:
% - reMatBin1, imMatBin1: Real and imaginary binary matrices of the first feature vector set.
% - mask1: Mask matrix for the first feature vector set, indicating valid regions.
% - reMatBin2, imMatBin2: Real and imaginary binary matrices of the second feature vector set.
% - mask2: Mask matrix for the second feature vector set, indicating valid regions.
% 
% Outputs:
% - minHamDist: The minimum Hamming distance computed between the two feature vector sets.
%
% Reference:
% - Zhang, D., Kong, W., You, J., & Wong, M. (2003). "Online Palmprint Identification". 
%   IEEE Transactions on Pattern Analysis and Machine Intelligence. 

% Get the size of the input matrices
matSz = size(mask1); % Assuming all input matrices have the same dimensions

% Define translation range dynamically for flexibility
translationRange = -2:2;
numTranslations = length(translationRange);

% Define translation ranges for rows (dim1) and columns (dim2)
dim1 = {3:matSz(1), 2:matSz(1), 1:matSz(1), 1:matSz(1)-1, 1:matSz(1)-2}; % Vertical shifts [-2, -1, 0, +1, +2]
dim2 = {3:matSz(2), 2:matSz(2), 1:matSz(2), 1:matSz(2)-1, 1:matSz(2)-2}; % Horizontal shifts [-2, -1, 0, +1, +2]

% Initialize Hamming distance matrix to store results for each translation
hammingDist = ones(numTranslations, numTranslations); % Initialize with ones as default high values for normalized distances

% Check if both masks fully cover the matrix area
if all(mask1(:)) && all(mask2(:))
    % Case 1: Full coverage (no invalid regions in the masks)
    for i = 1:numTranslations
        for j = 1:numTranslations
            % Translate matrices according to the offsets (i,j)
            reMatBin1Trans = reMatBin1(dim1{i}, dim2{j});
            imMatBin1Trans = imMatBin1(dim1{i}, dim2{j});
            reMatBin2Trans = reMatBin2(dim1{end-i+1}, dim2{end-j+1});
            imMatBin2Trans = imMatBin2(dim1{end-i+1}, dim2{end-j+1});
            
            % Compute Hamming distance without masks
            numerator = xor(reMatBin1Trans, reMatBin2Trans) + ...
                        xor(imMatBin1Trans, imMatBin2Trans);
            denominator = 2*numel(reMatBin1Trans); % Total number of bits
            hammingDist(i,j) = sum(numerator(:))/denominator; % Normalized Hamming distance
        end
    end
else
    % Case 2: Partial coverage (consider masks)
    for i = 1:numTranslations
        for j = 1:numTranslations
            % Translate matrices and masks according to the offsets (i,j)
            reMatBin1Trans = reMatBin1(dim1{i}, dim2{j});
            imMatBin1Trans = imMatBin1(dim1{i}, dim2{j});
            mask1Trans     =     mask1(dim1{i}, dim2{j});
            reMatBin2Trans = reMatBin2(dim1{end-i+1}, dim2{end-j+1});
            imMatBin2Trans = imMatBin2(dim1{end-i+1}, dim2{end-j+1});
            mask2Trans     =     mask2(dim1{end-i+1}, dim2{end-j+1});
            
            % Compute Hamming distance with masks applied
            validMask = mask1Trans & mask2Trans; % Overlap of both masks
            numerator = (validMask & xor(reMatBin1Trans, reMatBin2Trans)) + ...
                        (validMask & xor(imMatBin1Trans, imMatBin2Trans));
            denominator = 2*sum(validMask(:)); % Total valid bits based on mask overlap
            if denominator > 0
                hammingDist(i,j) = sum(numerator(:))/denominator; % Normalized Hamming distance
            else
                hammingDist(i,j) = 1; % Assign maximum distance if no overlap
            end
        end
    end
end

% Extract the minimum Hamming distance from all translations
minHamDist = min(hammingDist(:));