% ----------------------------------------------------------------------------------------------------
% Author        : Pourya Khodagholipour (P.KH)
% Project Name  : Implementation of a Palmprint Recognition System based on Image Processing (M.S. Thesis)
% File Name     : BioFigures.m
% Description   : Computes biometric performance metrics and evaluation graphs, including
%                 genuine/imposter distributions, ROC, and FAR/FRR curves
% Creation Date : 2010/01
% Revision Date : 2025/02/10
% ----------------------------------------------------------------------------------------------------

function [area, cutoff] = BioFigures(results)

% BioFigures computes biometric performance metrics and visualizations.
%
% This function processes genuine and imposter match data to generate:
% - Hamming distance distributions (Genuine vs. Imposter)
% - ROC (Receiver Operating Characteristic) curve
% - FAR (False Acceptance Rate) vs. FRR (False Rejection Rate)
% It also calculates the Area Under the Curve (AUC) and determines the optimal threshold (cutoff point).
%
% Input:
% - results : Nx2 matrix where:
%             Column 1: Hamming distance values
%             Column 2: Label (1 = Genuine match, 0 = Imposter match)
%
% Outputs:
% - area    : Area Under the Curve (AUC) for the ROC curve
% - cutoff  : Optimal threshold that balances FAR and FRR
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

% ========================================================================
% ------------------------------ Genuine & Imposter Matrices

% Extract genuine and imposter match results based on the last column indicator
genuine  = results(results(:,2)==1,1); % All genuine matches
imposter = results(results(:,2)==0,1); % All imposter matches

% Compute mean Hamming distances for genuine and imposter matches
genDistMean = mean(genuine);
impDistMean = mean(imposter);

% ========================================================================
% ------------------------------ Genuine & Imposter Distributions

% Define the X-axis range for histogram bins
xAxis = 0:0.002:1;

% -------------------- Plot genuine distribution
[counts, edges] = histcounts(genuine, xAxis); % Histogram bin counts
centers = edges(1:end-1) + diff(edges)/2; % Compute bin centers
counts = 100*counts/numel(genuine); % Normalize counts as percentage
figure, plot(centers, counts, 'b-', 'LineWidth', 1), hold on
xlim([0 0.6])

% -------------------- Plot imposter distribution
[counts, edges] = histcounts(imposter, xAxis);
centers = edges(1:end-1) + diff(edges)/2;
counts = 100*counts/numel(imposter);
plot(centers, counts, 'r--', 'LineWidth', 1), hold off

% Configure plot labels and legend
xlabel('Hamming Distance')
ylabel('Percentage (%)')
title('Normalized Hamming Distance Distributions (Genuine vs. Imposter)')
legend('Genuine', 'Imposter', 'Location', 'northwest')

% ========================================================================
% ------------------------------ Computing Area Under the Curve and Cutoff Point

% GA = Genuine Acceptance
% GR = Genuine Rejection
% FA = False Acceptance 
% FR = False Rejection
%
% genDistMean < impDistMean
%   GA : Genuine  & Hamming distance <  Thershold
%   GR : Imposter & Hamming distance >= Thershold
%   FA : Imposter & Hamming distance <  Thershold
%   FR : Genuine  & Hamming distance >= Thershold
% genDistMean >= impDistMean
%   GA : Genuine  & Hamming distance >  Thershold
%   GR : Imposter & Hamming distance <= Thershold
%   FA : Imposter & Hamming distance >  Thershold
%   FR : Genuine  & Hamming distance <= Thershold

% -------------------- Compute Receiver Operating Characteristic (ROC) values
uniqDist = unique(results(:,1)); % Unique Hamming distances
uniqDistNo = numel(uniqDist);
ROC = zeros(uniqDistNo,3); % Preallocate matrix for ROC values

for i = 1:uniqDistNo
    if genDistMean < impDistMean
        GA = sum(genuine  <  uniqDist(i));
        GR = sum(imposter >= uniqDist(i));
        FA = sum(imposter <  uniqDist(i));
        FR = sum(genuine  >= uniqDist(i));
    else
        GA = sum(genuine  >  uniqDist(i));
        GR = sum(imposter <= uniqDist(i));
        FA = sum(imposter >  uniqDist(i));
        FR = sum(genuine  <= uniqDist(i));
    end
    % Compute False Acceptance Rate (FAR), Genuine Acceptance Rate (GAR), and False Rejection Rate (FRR)
    totalGenuine  = GA + FR; % Total genuine attempts
    totalImposter = GR + FA; % Total imposter attempts
    FAR = FA / totalImposter; % False Acceptance Rate
    GAR = GA / totalGenuine; % Genuine Acceptance Rate
    FRR = FR / totalGenuine; % False Rejection Rate
    ROC(i,:) = [FAR, GAR, FRR];
end

% Flip ROC data if necessary
if genDistMean >= impDistMean
    ROC = flipud(ROC);
    uniqDist = flipud(uniqDist);
end

% Specify ROC points
xROC = [0 ; ROC(:,1) ; 1]; % False Acceptance Rate (FAR)
yROC = [0 ; ROC(:,2) ; 1]; % Genuine Acceptance Rate (GAR)
zROC = [0 ; ROC(:,3) ; 1]; % False Rejection Rate (FRR)

% -------------------- Compute AUC and optimal cutoff threshold
% Compute Area Under the Curve (AUC) using numerical integration
area = trapz(xROC, yROC);

% Determine the optimal cutoff point as the closest point to (0,1)
Dist = sqrt(xROC.^2 + (1-yROC).^2); % Euclidean distance from (0,1)
[~, id] = min(Dist); % Find the index of the optimal threshold
cutoff = uniqDist(id);

% ========================================================================
% ------------------------------ Plotting Graphs

% -------------------- ROC (Receiver Operating Characteristic) Curve
xROCuniq = unique(xROC);
yROCuniq = zeros(size(xROCuniq));
for i = 1:length(xROCuniq)
    yROCuniq(i) = max(yROC(xROC==xROCuniq(i)));
end
% Specify the y-coordinate of the ROC data below a certain value for elimination
idLow = find(yROCuniq<0.7, 1, 'last');
if isempty(idLow), idLow = 0; end

figure
semilogx(100*xROCuniq(idLow+1:end), 100*yROCuniq(idLow+1:end), 'b.-', 'LineWidth', 1.5);
hold on
plot(100*xROC(id), 100*yROC(id), 'rx', 'LineWidth', 2);
hold off
xlabel('False Acceptance Rate (%)')
ylabel('Genuine Acceptance Rate (%)')
title('ROC Curve')
legend('ROC Curve', 'Cutoff Point', 'Location', 'southeast')
ylim([0 100])

% -------------------- FAR (False Acceptance Rate) & FRR (False Rejection Rate) Plot
figure
plot(uniqDist, xROC(2:end-1), 'g', 'LineWidth', 1.5) % FAR
hold on
plot(uniqDist, zROC(2:end-1), 'b', 'LineWidth', 1.5) % FRR
hold off
xlabel('Threshold')
title('FAR & FRR')
legend('FAR', 'FRR', 'Location', 'north')
axis square