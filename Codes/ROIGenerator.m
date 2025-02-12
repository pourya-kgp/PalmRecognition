% ----------------------------------------------------------------------------------------------------
% Author        : Pourya Khodagholipour (P.KH)
% Project Name  : Implementation of a Palmprint Recognition System based on Image Processing (M.S. Thesis)
% File Name     : ROIGenerator.m
% Description   : Reads images from the specified database and extracts the Region of Interest (ROI)
%                 images based on the selected method (Morphological-based or Tangent-based).
% Creation Date : 2010/01
% Revision Date : 2025/02/10
% ----------------------------------------------------------------------------------------------------

% This script reads palmprint images from three different databases (PolyU, Delhi, and CASIA),
% extracts the Regions of Interest (ROIs), and saves them in a specified directory.
% It provides several options for customization and analysis:
%
% 1. ROI Extraction Method: The user can select between the Morphological-based or
%    Tangent-based method for segmentation.
% 2. Database Selection: The user can choose between three available palmprint databases:
%    PolyU, Delhi, or CASIA.
% 3. Progress Monitoring: The script provides feedback on the processed images and
%    identifies any ROI extractions that were unsuccessful.
% 4. Unsaved ROI Detection: The user can enable a flag to determine how many and which
%    ROIs were not saved successfully.
% 5. Visualization Option: The user can choose whether to display figures illustrating
%    different steps of the ROI extraction process for analysis and debugging.
% 6. Success Report: At the end of the code, a report is included on the performance of
%    the selected segmentation method, detailing the number of images that could not be
%    processed correctly for each database.
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

% Select the desired database for ROI extraction:
% - 'polyu' -> Hong Kong Polytechnic University (PolyU) Palmprint Database (The Second Version)
% - 'delhi' -> IIT Delhi Touchless Palmprint Database version 1.0
% - 'casia' -> CASIA Palmprint Image Database
database = 'polyu'; % Options: 'polyu', 'delhi', 'casia'

% ========================================================================
% ------------------------------ Choose the desired process

% Specify the desired process mode:
% - true  -> Skips preprocessing and identifies unsaved files without reprocessing
% - false -> Processes new data and extracts ROIs
unsavedFlag = false; % Options: true, false

% ========================================================================
% ------------------------------ Choose the visualization option

% Specify whether to display figures during execution:
% - true  -> The figure will be shown during execution
% - false -> The process runs without displaying it
%
% Note: Enable this option only for analyzing a limited number of ROIs,  
% as displaying too many images may cause system crashes.
figureFlag = false; % Options: true, false

% Morphological-based method:
% This flag controls whether to display two figures showing the stages of the finger
% detection process & ROI alignment and cropping. In the first figure, the process
% of finding the fingers, as well as the detection of the index and little fingers,
% is presented in multiple subplots. In the second figure, the hand image, coordinate
% system with the ROI boundaries, rotated fingers with the ROI boundaries, and the
% segmented ROI are presented. The option to show one or both of these figures can
% be found on the top of the PreMorphological.m function.
%
% Tangent-based method:
% This flag controls whether to display a figure showing the processing results.  
% In this figure, the hand image, tangent line with the ROI boundaries, rotated
% hand image with the ROI boundaries, and the segmented ROI are presented.  

% ========================================================================
% ------------------------------ Configure paths and parameters for the chosen method & database

% Common base path for datasets ("Data" directory path)
ComPath = 'D:\GitHub\PalmRecognition\Data\';

% Define input/output directories based on the selected method
readDir = 'Original';
saveDir = sprintf('Pre%s', method);

% Set database-specific paths and parameters
switch database
    case 'polyu' % --------------- Polytech
        readPath = fullfile(ComPath, 'PolyU', readDir); % Input directory path
        savePath = fullfile(ComPath, 'PolyU', saveDir); % Output directory path
        persons = 386; % Number of individuals
        samples = 17; % Maximum samples per individual
        spec = 'FS'; % Specifications: F = First Edition, S = Second Edition
    
    case 'delhi' % --------------- Delhi
        readPath = fullfile(ComPath, 'Delhi', readDir); % Input directory path
        savePath = fullfile(ComPath, 'Delhi', saveDir); % Output directory path
        persons = 234; % Number of individuals
        samples = 11; % Maximum samples per individual
        spec = {'Left  Hand', 'Right Hand'}; % Hand specifications
    
    case 'casia' % --------------- CASIA
        readPath = fullfile(ComPath, 'CASIA', readDir); % Input directory path
        savePath = fullfile(ComPath, 'CASIA', saveDir); % Output directory path
        persons = 312; % Number of individuals
        samples = 18; % Maximum samples per individual
        spec = ['ml'; 'mr'; 'fl'; 'fr']; % Specifications: m = Male, f = Female, l = Left, r = Right
    
    otherwise
        error('Unsupported database. Choose from "polyu", "delhi", or "casia".');
end

% ========================================================================
% ------------------------------ Function handler for ROI extraction 

% Assigns the appropriate function for the chosen segmentation method.
switch method
    case 'Morph'   % --------------- Morphological-based method
        roiExtractor = @(hand, database, figureFlag) ...
                         PreMorphological(hand, database, figureFlag);
    case 'Tangent' % --------------- Tangent-based method
        roiExtractor = @(hand, database, figureFlag) ...
                         PreTangent(hand, database, figureFlag);
    otherwise
        error('Unsupported method. Choose from "Morph" or "Tangent".');
end

% ========================================================================
% ------------------------------ Read the image file, preprocess it, and save it

roi = []; % Initialize ROI variable
readFiles = 0; % Counter for successfully processed images
unsavedFiles = 0; % Counter for unsaved ROIs due to segmentation failures

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
                    Fname = sprintf('%04d\\%04d_%c_%c_%02d.jpg', i, i, spec(j,1), spec(j,2), k);
            end

            % Construct full file path and verify existence
            filePath = fullfile(readPath, Fname);
            if ~isfile(filePath)
                continue; % Skip to the next iteration if the file does not exist
            end

            % Read the image if preprocessing is not skipped and preprocess it
            if ~unsavedFlag
                hand = imread(filePath);
                if strcmp(database, 'delhi') % Rotate images for Delhi database
                    hand = imrotate(hand, 90);
                end
                roi = roiExtractor(hand, database, figureFlag); % Preprocess the image
            end
            
            % Construct save path
            % Adjust file name for CASIA for save operation
            if strcmp(database, 'casia')
                Fname = Fname(6:end);
            end
            % Construct full file path for save operation
            filePath = fullfile(savePath, Fname);
            
            % Save the extracted ROI if it is valid
            if isempty(roi) ~= 1
                switch database
                    case 'polyu' % --------------- Polytech
                        imwrite(roi, filePath, 'bmp')
                    case 'delhi' % --------------- Delhi
                        imwrite(roi, filePath, 'bmp')
                    otherwise    % --------------- CASIA
                        imwrite(roi, filePath, 'jpg')
                end
            end

            % Check if the ROI has been saved successfully (Log success or failure)
            if ~isfile(filePath)
                fprintf('Method: %s-based, Database: %s, ROI: %s ==> ROI not saved\n', ...
                         method, database, Fname);
                unsavedFiles = unsavedFiles + 1; % Increment unsaved counter
            else
                if ~unsavedFlag
                    fprintf('Method: %s-based, Database: %s, ROI: %s ==> ROI saved\n', ...
                             method, database, Fname);
                end
            end
            readFiles = readFiles + 1; % Increment read counter

        end
    end
end

% Display summary of processing
fprintf('Total processed images: %d\nTotal unsaved ROIs: %d\n', readFiles, unsavedFiles);

% Cleanup
clearvars

% ========================================================================
% ------------------------------ The Report (Morphological-based method)

% -------------------- PolyU database results
% 9 of 7,752 palmprint images can not be segmented properly. 
% %0.012 of the database images are not recognizable at the preprocessing stage.
%
% The 9 images that can not be segmented are:
% PolyU_004_F_07/09/10 | PolyU_022_F_05 | PolyU_122_S_03/09 | PolyU_162_F_04 | PolyU_200_F_10 |
% PolyU_225_F_07 |

% -------------------- delhi database results
% 23 of 2,769 palmprint images can not be segmented properly.
% %0.83 of the database images are not recognizable at the preprocessing stage.
%
% The 23 images that can not be segmented are:
% Left Hand\  ... 008_3/4/5/6/7 | 061_1 | 107_2 | 138_4/6 | 163_6 | 190_1/2 |
% Right Hand\ ... 014_5 | 041_1/2/4/5/6 | 054_6 | 069_6 | 140_4 | 201_3/4 |

% -------------------- CASIA database results
% 167 of 5,486 palmprint images can not be segmented properly.
% %3.0 of the database images are not recognizable at the preprocessing stage.
%
% The 167 images that can not be segmented are:
% 0001 ~ 0099 (24 images)
% 0003_m_r_05 | 0030_m_l_05 | 0030_m_r_01/~/08 | 0036_m_l_03 | 0052_f_r_04/06/07 | 0059_m_l_06 |
% 0076_f_l_03 | 0076_f_r_01/02 | 0077_m_r_04/05/06 | 0078_m_l_07/08 | 0093_m_l_03 |
% 0100 ~ 0199 (46 images)
% 0101_m_l_09/10 | 0116_m_l_07/08/09 | 0117_f_r_01/~/08 | 0135_m_r_01/02/03/05/~/08 |
% 0142_m_l_06 | 0151_m_r_04 | 0166_f_l_01/~/04 | 0166_f_r_04/~/09 | 0185_m_r_01/~/08 |
% 0188_f_r_01 | 0189_m_l_07 | 0191_m_r_07/08 | 0193_m_r_08 | 0199_m_l_01 |
% 0200 ~ 0299 (29 images)
% 0214_m_l_01/~/05/08/09/10 | 0215_m_l_08 | 0234_m_l_02 | 0235_m_l_01/06/07/10 | 0237_m_l_06 |
% 0237_m_r_05 | 0247_m_r_09 | 0259_m_r_05 | 0261_m_l_01/02/03 | 0270_m_r_11 | 0285_m_l_02/03/08 |
% 0285_m_r_10 | 0290_m_l_04 | 0292_m_r_01/05 |
% 0300 ~ 0312 (68 images)
% 0304_m_l_05/06/07 | 0304_m_r_06 | 0306_f_l_01/~/04/09/~/12 | 0306_f_r_01/~/06/09/~/14 |
% 0307_m_l_01/02/04/~/18 | 0307_m_r_01/~/09/12 | 0308_f_l_02 | 0308_f_r_04/07/~/11 |
% 0309_m_r_12/~/15 | 0310_m_l_01/~/04 | 0311_m_l_01 | 0311_m_r_01 |

% ========================================================================
% ------------------------------ The Report (Tanget-based Method)

% -------------------- PolyU database results
% 26 of 7,752 palmprint images can not be segmented properly.
% %0.034 of the database images are not recognizable at the preprocessing stage.
%
% The 24 images that can not be segmented are:
% PolyU_065_S_01/04/07/08/09 | PolyU_073_S_08/09/10 | PolyU_109_F_02/03/04/05/08/09/10 |
% PolyU_139_F_03/04/05       | PolyU_245_F_01       | PolyU_257_S_06                   |
% PolyU_287_F_05/09/10       | PolyU_293_F_09/10    | PolyU_379_F_09                   |

% -------------------- delhi database results
% 48 of 2,769 palmprint images can not be segmented properly.
% %1.7 of the database images are not recognizable at the preprocessing stage.
%
% The 48 images that can not be segmented are:
% Left Hand\  ... 026_5 | 035_6 | 040_2/4 | 051_6 | 058_7 | 077_4 | 078_6 | 136_1/~/4 |
% Left Hand\  ... 138_1/3 | 141_3/9 | 155_3 | 190_2 | 212_5/6 | 216_2 | 226_3 |
% Right Hand\ ... 033_4 | 040_3 | 045_4 | 058_3 | 111_3/4 | 121_5 | 130_1 | 131_4 | 132_2 |
% Right Hand\ ... 133_3 | 136_4 | 138_1/~/5 | 140_4 | 141_4/~/8 | 212_5/6 | 226_6 |

% -------------------- CASIA database results
% 381 of 5,486 palmprint images can not be segmented properly.
% %6.9 of the database images are not recognizable at the preprocessing stage.
%
% The 381 images that can not be segmented are:
% 0001 ~ 0099 (67 images)
% 0001_m_l_02 | 0003_m_r_01/05 | 0012_m_r_02 | 0014_m_r_07/08 | 0017_m_r_07 | 0021_m_r_04/05 |
% 0027_f_r_01/02/03/06/07/08 | 0030_m_l_05 | 0030_m_r_01/~/08 | 0034_f_r_05 | 0036_m_l_06 |
% 0037_f_r_01/03 | 0039_m_r_08 | 0040_m_l_05 | 0044_m_r_04/05 | 0046_f_r_06 | 0047_m_r_04 |
% 0048_m_r_06 | 0049_m_r_06/08 | 0050_m_r_08 | 0052_f_l_03 | 0053_m_r_01/02/03/06/07 |
% 0055_m_l_01/02 | 0055_m_r_05/06 | 0056_m_r_06 | 0064_m_r_06 | 0066_m_l_04 | 0069_f_r_08 |
% 0073_m_l_05 | 0076_f_l_01/03/05/~/08 | 0076_f_r_01/02 | 0079_m_r_03/05 | 0085_f_l_01/06/07 |
% 0089_m_l_01 |
% 0100 ~ 0199 (105 images)
% 0101_m_l_01/~/10 | 0101_m_r_07/~/10 | 0105_m_r_01 | 0110_m_r_04 | 0113_m_r_02/07/08 |
% 0117_f_l_03/04 | 0119_m_l_04 | 0123_f_l_01/~/04/08 | 0127_m_r_01/02 | 0129_m_r_04 | 0130_m_r_01 |
% 0133_m_l_03 | 0137_m_r_02 | 0143_m_l_04 | 0143_m_r_01/02 | 0144_m_r_01/05/08 | 0145_f_l_03 |
% 0145_f_r_05 | 0147_m_l_02 | 0147_m_r_01 | 0149_m_r_08 | 0150_f_r_02/03 | 0153_f_l_07/08 |
% 0154_f_l_06 | 0155_f_r_01/04/05 | 0159_f_l_03 | 0162_m_r_06 | 0164_m_l_04 | 0166_f_l_04 |
% 0166_f_r_04/07/09 | 0168_m_l_01/02/03 | 0168_m_r_01/~/07 | 0170_m_l_06 | 0172_m_l_01/02 |
% 0182_f_r_01 | 0185_m_l_01 | 0188_f_r_01 | 0189_m_r_01/~/05 | 0190_f_l_01 | 0191_m_r_01/~/08 |
% 0192_f_l_07/08 | 0195_f_l_01/02 | 0196_f_r_02/03 | 0197_f_l_01/04 | 0199_m_r_01/03/~/06/08 |
% 0200 ~ 0299 (105 images)
% 0201_f_r_03/04/05 | 0202_f_r_05 | 0207_m_l_03 | 0211_m_l_07 | 0214_m_l_03/04/06/~/09 |
% 0218_m_l_05 | 0220_m_l_01 | 0223_f_l_01/02 | 0225_m_l_01/02/03 | 0226_m_l_01/~/06 |
% 0232_m_r_04 | 0233_m_l_08 | 0234_m_r_04/08 | 0236_m_l_01/02 | 0237_m_r_09 | 0238_m_l_02 |
% 0239_m_l_05 | 0240_m_l_01/~/04/06 | 0248_m_l_01 | 0249_m_l_01/02/03 | 0250_m_r_04/05/09 |
% 0251_m_l_01/04 | 0254_m_r_01 | 0255_m_l_07/10 | 0256_m_l_06 | 0258_m_l_01/02/08/09 |
% 0259_m_l_01/02 | 0267_m_l_02/~/05 | 0269_m_l_01/~/09 | 0269_m_r_01 | 0270_m_r_05 |
% 0278_m_l_01/02 | 0280_m_l_07 | 0283_m_l_09/10 | 0283_m_r_06/08/09/10 | 0285_m_l_02/03 |
% 0287_m_l_08/09/10 | 0287_m_r_03/04 | 0288_m_l_02 | 0290_m_l_01/~/05/07 | 0290_m_r_03 |
% 0291_m_r_07 | 0294_m_r_01 | 0295_m_l_03 | 0295_m_r_08 | 0297_m_r_01 | 0299_m_l_01/02/06 |
% 0300 ~ 0312 (104 images)
% 0302_m_l_09/10 | 0304_m_l_04/05/06/10 | 0304_m_r_06 | 0306_f_l_01/~/13 | 0306_f_r_01/~/14 |
% 0307_m_l_01/~/18 | 0307_m_r_01/~/13 | 0308_f_l_01 | 0308_f_r_01/02/04/05/07/~/11 |
% 0309_m_l_01/03/08/12 | 0309_m_r_01/04/07/~/13/15 | 0310_m_l_01/~/06/08/~/11 | 0310_m_r_11 |
% 0311_m_l_01/11 | 0311_m_r_01/11 |

% ========================================================================
% ------------------------------ Databases

% -------------------- PolyU database
% In the tangent-based method, most of the unsuccessful preprocessing in this database is 
% the cause of changing the algorithm to be capable of segmenting three different databases
% with a single method. In the case of the PolyU database, changing the method of finding
% the first point of the boundary tracing function will solve this problem. The substitute
% algorithm is as simple as finding the two nearest edge points to the top left and bottom
% left corners of the image to achieve the first points for gap1 and gap2, respectively.

% -------------------- delhi database
% This database has pictures of the whole hand with all the fingers and sometimes bracelets or
% watches. The most challenging problem in this database is the image binarizing because the
% hand images were taken in a semi-controlled environment. This means the illumination varies
% extremely. Therefore, converting the image from gray to binary with a single threshold values
% is not the best solution. Even the Otsu method does not result in acceptable hand images.

% -------------------- CASIA database
% This database is the most challenging one among the others. Images do not contain the whole
% hand. Therefore, all the fingers are not recognizable in them. Also, in some of them, a tiny
% portion of the little finger is included, which makes it hard to specify the hand orientation
% or the four fingers. In this database, the hand orientation varies extremely, and with the
% incomplete hand image, finding the four fingers except the thumb seems difficult with this
% method. Another problem in this database is image binarizing. Although the images were 
% captured in a controlled environment, the images' illumination varies tremendously.
% Therefore, converting the image from gray to binary with the global threshold values is
% not the best solution. Even the Otsu method does not result in acceptable hand images.