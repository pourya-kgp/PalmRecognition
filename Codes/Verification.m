% ----------------------------------------------------------------------------------------------------
% Author        : Pourya Khodagholipour (P.KH)
% Project Name  : Implementation of a Palmprint Recognition System based on Image Processing (M.S. Thesis)
% File Name     : Verification.m
% Description   : Constructs the result matrix for verification assessment, and visualizes results
% Creation Date : 2010/01
% Revision Date : 2025/02/10
% ----------------------------------------------------------------------------------------------------

% Verification Assessment for Palmprint Recognition
%
% This script evaluates the verification performance of a palmprint recognition system
% across three different databases (PolyU, Delhi, and CASIA) using the selected
% segmentation method (Morphological-based or Tangent-based).
% It constructs a result matrix, evaluates biometric parameters, and visualizes:
% - Genuine vs. Imposter Distributions
% - Receiver Operating Characteristic (ROC) Curve
% - False Acceptance Rate (FAR) & False Rejection Rate (FRR)
% - Area Under the Curve (AUC)
% The results are saved periodically to allow for resumption in case of long computations.
% It provides several options for customization and analysis:
%
% 1. ROI Extraction Method: The user can select between the Morphological-based or
%    Tangent-based method for segmentation.
% 2. Database Selection: The user can choose between three available palmprint databases:
%    PolyU, Delhi, or CASIA.
% 3. Computation Mode:
%    - Compute fresh results from scratch.
%    - Resume calculations from previously saved results, allowing continuation when
%      long processing is required.
% 4. Visualization Mode: The user can choose to display previously computed results
%    without recomputation.
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
% ------------------------------ Choose the desired process

% Specify the desired process mode:
% - true  -> Continues from the last saved results
% - false -> Constructs a new result
continueFlag = true; % Options: true, false

% ========================================================================
% ------------------------------ Choose the visualization option

% Specify the visualization mode:
% - true  -> Visualizes results without recomputation
% - false -> Runs verification process
figureFlag = true; % Options: true, false

% ========================================================================
% ------------------------------ Configure paths and parameters for the chosen method & database

% Common base path for datasets ("Data" directory path)
ComPath = 'D:\GitHub\PalmRecognition\Data\';

% Define input/output directories based on the selected method
readDir = sprintf('Pre%s_Feature_Vectors', method);
saveDir = sprintf('Pre%s_Results', method);

% Set database-specific paths and parameters
switch database
    case 'polyu' % --------------- Polytech
        readPath = fullfile(ComPath, 'PolyU', readDir); % Input directory path
        savePath = fullfile(ComPath, 'PolyU', saveDir); % Output directory path
        persons = 386; % Number of individuals
        samples = 17; % Maximum samples per individual
        spec = 'FS'; % Specifications: F = First Edition, S = Second Edition
        maxPalm = 7752; % Maximum number of palmprints in the database
        
    case 'delhi' % --------------- Delhi
        readPath = fullfile(ComPath, 'Delhi', readDir); % Input directory path
        savePath = fullfile(ComPath, 'Delhi', saveDir); % Output directory path
        persons = 234; % Number of individuals
        samples = 11; % Maximum samples per individual
        spec = {'Left  Hand', 'Right Hand'}; % Hand specifications
        maxPalm = 2769; % Maximum number of palmprints in the database

    case 'casia' % --------------- CASIA
        readPath = fullfile(ComPath, 'CASIA', readDir); % Input directory path
        savePath = fullfile(ComPath, 'CASIA', saveDir); % Output directory path
        persons = 312; % Number of individuals
        samples = 18; % Maximum samples per individual
        spec = ['ml'; 'mr'; 'fl'; 'fr']; % Specifications: m = Male, f = Female, l = Left, r = Right
        maxPalm = 5486; % Maximum number of palmprints in the database

    otherwise
        error('Unsupported database. Choose from "polyu", "delhi", or "casia".');
end

% Ensure the selected method is valid
if ~(strcmp(method, 'Morph') || strcmp(method, 'Tangent'))
    error('Unsupported method. Choose from "Morph" or "Tangent".');
end

% The File name for saving the results
Fsave = 'VerificationResults.mat';

% ========================================================================
% ------------------------------ Matching and Saving

% Predefine the result matrix
maxSize = nchoosek(maxPalm, 2);
results = inf(maxSize, 8);

% Predefine the individual who algorithm starts at
if figureFlag
    startPerson = persons;
    startSpec   = length(spec);
    startSample = samples;
else
    startPerson = 1;
    startSpec   = 1;
    startSample = 1;
end

% Cunter for the current row of the result matrix
matRow = 0;

% Load previous results if required (Setting the figure flag or continue flag)
% Assign value to the correspondig variables (results matrix, matRow, startPerson)
% Construct full file path
filePath = fullfile(savePath, Fsave);
if (figureFlag || continueFlag) && isfile(filePath)
    load(filePath, 'results', 'matRow');
    if ~figureFlag && continueFlag
        startPerson = results(matRow,1) + 1;
    end
end

% Loop over all individuals, hand specifications, and samples
for i = startPerson:persons         % Individual No.
    for j = startSpec:length(spec)  % Hand specifications
        for k = startSample:samples % Sample No.
            
            % Generate the file name based on the database structure            
            switch database
                case 'polyu' % --------------- Polytech
                    Fname1 = sprintf('PolyU_%03d_%c_%02d.mat', i, spec(j), k);
                case 'delhi' % --------------- Delhi
                    Fname1 = sprintf('%s\\%03d_%d.mat', spec{j}, i, k);
                otherwise    % --------------- CASIA
                    Fname1 = sprintf('%04d_%c_%c_%02d.mat', i, spec(j,1), spec(j,2), k);
            end

            % Construct full file path and verify existence
            filePath = fullfile(readPath, Fname1);
            if ~isfile(filePath)
                continue; % Skip to the next iteration if the file does not exist
            end

            % Load the feature matrices and their mask
            load(filePath, 'reGaborBin', 'imGaborBin', 'mask');
            reMatBin1 = reGaborBin;
            imMatBin1 = imGaborBin;
            mask1     = mask;
            
            % Compare with the rest of individuals, hand specifications, and samples
            for ii = i:persons          % Individual No.
                for jj = 1:length(spec) % Hand specifications
                    for kk = 1:samples  % Sample No.

                        % Skip the repeated comparison for the ith person
                        if ii == i
                            if jj < j || (jj == j && kk <= k)
                                continue;
                            end
                        end

                        % Generate the file name based on the database structure
                        switch database
                            case 'polyu' % --------------- Polytech
                                Fname2 = sprintf('PolyU_%03d_%c_%02d.mat', ii, spec(jj), kk);
                            case 'delhi' % --------------- Delhi
                                Fname2 = sprintf('%s\\%03d_%d.mat', spec{jj}, ii, kk);
                            otherwise    % --------------- CASIA
                                Fname2 = sprintf('%04d_%c_%c_%02d.mat', ii, spec(jj,1), spec(jj,2), kk);
                        end
            
                        % Construct full file path and verify existence
                        filePath = fullfile(readPath, Fname2);
                        if ~isfile(filePath)
                            continue; % Skip to the next iteration if the file does not exist
                        end
                        
                        % Load the feature matrices
                        load(filePath, 'reGaborBin', 'imGaborBin', 'mask');
                        reMatBin2 = reGaborBin;
                        imMatBin2 = imGaborBin;
                        mask2     = mask;
                        
                        % Calculate the minimum normalized hamming distance between the transfered images
                        minHamDist = HammingDistanceMasked (reMatBin1, imMatBin1, mask1, ...
                                                            reMatBin2, imMatBin2, mask2);
                        
                        % Calculate the correspondig row number to save the results
                        matRow = matRow + 1; % Increament the matRow cunter
                        
                        % i : person1,  j: Hand specification from person1,  k: Sample No. from person1 ==> Palmprint1
                        % ii: person2, jj: Hand specification from person2, kk: Sample No. from person2 ==> Palmprint2
                        % minhamDist: Match value between Palmprint1 and Palmprint2
                        % isMatch: genuine or imposter person
                        % Determine genuine or imposter match
                        switch database
                            case 'polyu' % --------------- Polytech
                                isMatch = i==ii;
                            otherwise    % --------------- Delhi & CASIA
                                isMatch = i==ii & j==jj;
                        end
                        results(matRow,:) = [i, j, k, ii, jj, kk, minHamDist, isMatch];
                        fprintf(['Method: %s-based, Database: %s, %s vs. %s ==> ' ...
                                 'Hamming Distance = %1.4f, Match:%d, No.:%d\n'], ...
                                  method, database, Fname1(1:end-4), Fname2(1:end-4), ...
                                  minHamDist, isMatch, matRow);
                    end
                end
            end
        end
    end
    % Save the results if the purpose is not to only see the figures
    if ~figureFlag
        % Construct full file path
        filePath = fullfile(savePath, Fsave);
        save(filePath, 'results', 'matRow')
    end
end

% Trim the inf part of the results matrix
results(matRow+1:end,:) = [];

% Depict the figures accossiated with the results (Geniune vs. Imposter histogram, ROC, FAR, FRR)
[area, cutoff] = BioFigures(results(:,7:8));

% Cleanup unnecessary variables
clearvars -except results matRow area cutoff;