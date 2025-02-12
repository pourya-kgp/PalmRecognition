% ----------------------------------------------------------------------------------------------------
% Author        : Pourya Khodagholipour (P.KH)
% Project Name  : Implementation of a Palmprint Recognition System based on Image Processing (M.S. Thesis)
% File Name     : FeatureVectorGenerator.m
% Description   : Reads palmprint ROI images from the specified database and extracts feature vectors
%                 using the selected segmentation method (Morphological-based or Tangent-based).
% Creation Date : 2010/01
% Revision Date : 2025/02/10
% ----------------------------------------------------------------------------------------------------

% This script extracts feature vectors from palmprint ROI images obtained
% from three different databases (PolyU, Delhi, and CASIA) using the selected
% segmentation method (Morphological-based or Tangent-based).
% It provides several options for customization and analysis:
%
% 1. ROI Extraction Method: The user can select between the Morphological-based or
%    Tangent-based method for segmentation.
% 2. Database Selection: The user can choose between three available palmprint databases:
%    PolyU, Delhi, or CASIA.
% 3. Progress Monitoring: The script provides feedback on the processed images and
%    identifies any feature extractions that were unsuccessful.
% 4. Saving Options: Feature matrices are always saved in `.mat` format, and binary
%    images can be saved optionally.
% 5. Visualization Option: The user can choose whether to display binary feature images
%    for analysis and debugging.
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

% Clear workspace and close figures
close all; clear; clc;

% ========================================================================
% ------------------------------ Choose the Method (Morphological-based / Tangent-based)

% Select the desired method for ROI extraction:
% - 'Morph'   -> Morphological-based method
% - 'Tangent' -> Tangent-based method
method = 'Morph'; % Options: 'Morph', 'Tangent'

% ========================================================================
% ------------------------------ Choose the database (PolyU / IIT Delhi / CASIA)

% Select the desired database for feature extraction:
% - 'polyu' -> Hong Kong Polytechnic University (PolyU) Palmprint Database (The Second Version)
% - 'delhi' -> IIT Delhi Touchless Palmprint Database version 1.0
% - 'casia' -> CASIA Palmprint Image Database
database = 'polyu'; % Options: 'polyu', 'delhi', 'casia'

% ========================================================================
% ------------------------------ Choose the desired outputs

% Specify the desired outputs:
% - true  -> Saves the binary feature images (real and imaginary parts) for visualization, 
%            along with the feature matrices in .mat format
% - false -> Saves only the feature matrices in .mat format, without binary images
imageFlag = false; % Options: true, false

% ========================================================================
% ------------------------------ Choose the visualization option

% Specify the visualization mode:
% - true  -> Displays the binary feature images (real and imaginary parts) for visualization
% - false -> Disables visualization
%  
% Note: Enable this option only for analyzing a limited number of ROIs,  
% as displaying too many images may cause system crashes.
figureFlag = false; % Options: true, false

% ========================================================================
% ------------------------------ Configure paths and parameters for the chosen method & database

% Common base path for datasets ("Data" directory path)
ComPath = 'D:\GitHub\PalmRecognition\Data\';

% Define input/output directories based on the selected method
readDir = sprintf('Pre%s', method);
saveDir = sprintf('Pre%s_Feature_Vectors', method);

% Set database-specific paths and parameters
switch database
    case 'polyu' % --------------- Polytech
        readPath = fullfile(ComPath, 'PolyU', readDir); % Input directory path
        savePath = fullfile(ComPath, 'PolyU', saveDir); % Output directory path
        persons = 386; % Number of individuals
        samples = 17; % Maximum samples per individual
        spec = 'FS'; % Specifications: F = First Edition, S = Second Edition
        vectSize = 2048; % Feature vector size (According to the paper)
        
    case 'delhi' % --------------- Delhi
        readPath = fullfile(ComPath, 'Delhi', readDir); % Input directory path
        savePath = fullfile(ComPath, 'Delhi', saveDir); % Output directory path
        persons = 234; % Number of individuals
        samples = 11; % Maximum samples per individual
        spec = {'Left  Hand', 'Right Hand'}; % Hand specifications
        vectSize = 5408; % Feature vector size
        
    case 'casia' % --------------- CASIA
        readPath = fullfile(ComPath, 'CASIA', readDir); % Input directory path
        savePath = fullfile(ComPath, 'CASIA', saveDir); % Output directory path
        persons = 312; % Number of individuals
        samples = 18; % Maximum samples per individual
        spec = ['ml'; 'mr'; 'fl'; 'fr']; % Specifications: m = Male, f = Female, l = Left, r = Right
        vectSize = 8192; % Feature vector size
        
    otherwise
        error('Unsupported database. Choose from "polyu", "delhi", or "casia".');
end

% Ensure the selected method is valid
if ~(strcmp(method, 'Morph') || strcmp(method, 'Tangent'))
    error('Unsupported method. Choose from "Morph" or "Tangent".');
end

% ========================================================================
% ------------------------------ Feature Extraction and Saving

% Loop over all individuals, hand specifications, and samples
for i = 1:persons          % Individual No.
    for j = 1:length(spec) % Hand specifications
        for k = 1:samples  % Sample No.

            % Generate the file name based on the database structure            
            switch database
                case 'polyu' % --------------- Polytech
                    Fname = sprintf('PolyU_%03d_%c_%02d.bmp', i, spec(j), k);
                case 'delhi' % --------------- Delhi
                    Fname = sprintf('%s\\%03d_%d.bmp', spec{j}, i, k);
                otherwise    % --------------- CASIA
                    Fname = sprintf('%04d_%c_%c_%02d.jpg', i, spec(j,1), spec(j,2), k);
            end

            % Construct full file path and verify existence
            filePath = fullfile(readPath, Fname);
            if ~isfile(filePath)
                continue; % Skip to the next iteration if the file does not exist
            end

            % Read the ROI image
            roi = imread(filePath);
            
            % Extract feature matrices and corresponding mask
            [reGaborBin, imGaborBin, mask] = FeatureExtraction(roi, vectSize, database, figureFlag);
            
            % Save binary feature images if visualization is enabled
            if imageFlag
                filePath = fullfile(savePath, [Fname(1:end-4)  '_Re' Fname(end-3:end)]);
                imwrite(reGaborBin, filePath, 'bmp')
                filePath = fullfile(savePath, [Fname(1:end-4)  '_Im' Fname(end-3:end)]);
                imwrite(imGaborBin, filePath, 'bmp')
            end

            % Save extracted feature vectors in `.mat` format
            filePath = fullfile(savePath, Fname(1:end-4));
            save(filePath, 'reGaborBin', 'imGaborBin', 'mask');
            
            % Check if the .mat extension has been saved successfully (Log success or failure)
            if ~isfile([filePath, '.mat'])
                fprintf('Method: %s-based, Database: %s, %s ==> Feature matrices not saved\n', ...
                         method, database, Fname(1:end-4));
            else
                fprintf('Method: %s-based, Database: %s, %s ==> Feature matrices saved\n', ...
                         method, database, Fname(1:end-4));
            end
            
        end
    end
end

% Cleanup
clearvars