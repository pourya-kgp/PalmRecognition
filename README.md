# PalmRecognition
Implementation of a Palmprint Recognition System based on Image Processing (M.S. Thesis)
****************************************************************************************************
## Overview
This project implemented a palmprint recognition system using image processing techniques. The system
followed these main steps:

1. Accurate implementation of Zhang's "Online Palmprint Identification" (2003):
   - Used the tangent-based method for Region-Of-Interest (ROI) segmentation (preprocessing).
   - Extracted features using a 2D Circular Gabor filter.
   - Performed matching using the Hamming distance metric.
   - Was initially tested on the PolyU database.

2. Enhancements to improve preprocessing across multiple databases:
   - Applied modifications to support PolyU, Delhi, and CASIA palmprint databases.
   - Optimized algorithm execution and result accuracy.

   Key Enhancements:
   - Hand Orientation Estimation: Developed a method to determine hand orientation before ROI extraction.
   - Boundary Tracing & Gap Selection: Improved starting point selection for contour tracing and enhanced
     finger gap selection criteria for improved tangent line computation.
   - ROI Assessment & Validation: Evaluates the extracted ROI based on non-palm area (black region) volume.

3. Development of a novel morphological-based ROI segmentation method:
   - Introduced morphological image processing techniques for accurate ROI extraction.
   - improved ROI extraction robustness across different databases by adapting to variations in hand 
     positioning and rotation. 

4. Evaluation of preprocessing methods across three databases:
   - Computed biometric performance metrics and generated visualizations, including:
     - Hamming distance distributions (Genuine vs. Imposter)
     - ROC (Receiver Operating Characteristic) curve
     - FAR (False Acceptance Rate) vs. FRR (False Rejection Rate)
****************************************************************************************************
## Included Research Papers
- Morphological-based segmentation paper:
  - Title      : A Novel Method for Locating the Region-Of-Interest in Palmprint Recognition Systems
  - Authors    : Pourya Khodagholipour, Hamed Sadjedi
  - Conference : IEEE PIC 2010 (Accepted but not registered)
  - File       : `2010_IEEE_PIC_ROI_Palmprint.pdf`

- Seminar on Palmprint Recognition (Persian/Farsi):
  - Title      : A Review on Image Processing Techniques for Palmprint Verification/Identification
  - Authors    : Pourya Khodagholipour, Hamed Sadjedi
  - File       : `2009_Seminar_Palmprint_Recognition_Research.pdf`

These documents can be found in the `Docs/` directory.
****************************************************************************************************
## M-Files Overview

### ROI Segmentation
- `ROIGenerator.m`          : Reads images from the database and extracts ROI images.
- `PreMorphological.m`      : Uses morphological-based segmentation for ROI extraction.
- `PreTangent.m`            : Uses tangent-based segmentation for ROI extraction.

### Feature Extraction
- `FeatureVectorGenerator.m`: Extracts feature vectors from ROI images.
- `FeatureExtraction.m`     : Extracts discriminative features from the ROI image.
- `CircularGaborFilter.m`   : Applies a Circular Gabor Filter (zero DC component).

### Matching & Results
- `Verification.m`          : Performs verification assessment and visualizes results.
- `Identification.m`        : Performs identification assessment and visualizes results.
- `BioFigures.m`            : Computes biometric performance metrics, including ROC, FAR, FRR.
- `HammingDistanceMasked.m` : Computes minimum Hamming distance between binary feature matrices.
****************************************************************************************************
## Data Directory Structure
This directory contains database-related files, including input images, processed results, and 
evaluation metrics. The structure includes three subdirectories for PolyU, Delhi, and CASIA databases.
Each of these subdirectories is structured as follows:

- `Original/`                   – Contains raw images from each database.
- `PreMorph/`                   – Stores morphological-based segmented ROI images.
- `PreMorph_Feature_Vectors/`   – Stores extracted feature vectors from morphological-based segmentation.
- `PreMorph_Figures/`           – Contains visualized results of morphological-based segmentation.
- `PreMorph_Results/`           – Contains verification & identification results from morphological-based segmentation.
- `PreTangent/`                 – Stores tangent-based segmented ROI images.
- `PreTangent_Feature_Vectors/` – Stores extracted feature vectors from tangent-based segmentation.
- `PreTangent_Figures/`         – Contains visualized results of tangent-based segmentation.
- `PreTangent_Results/`         – Contains verification & identification results from tangent-based segmentation.

- `Palmprint Database Description.pdf` – Provides detailed database information.
****************************************************************************************************
## References

### 1. Morphological-based segmentation method
- Title      : A Novel Method for Locating the Region-Of-Interest in Palmprint Recognition Systems
- Authors    : Pourya Khodagholipour, Hamed Sadjedi
- Conference : IEEE 2010 International Conference on Progress in Informatics and Computing (PIC)
- Link       : https://github.com/pourya-kgp/PalmRecognition/blob/main/Docs/2010_IEEE_PIC_ROI_Palmprint.pdf

### 2. Tangent-based segmentation, feature extraction & matching
- Title      : Online Palmprint Identification
- Authors    : D. Zhang, W. K. Kong, J. You, M. Wong
- Journal    : IEEE Transactions on Pattern Analysis and Machine Intelligence, Vol. 25, No. 9, pp. 1041-1050, 2003
- Link       : https://ieeexplore.ieee.org/document/1227981
****************************************************************************************************
## Databases
This project evaluates three palmprint databases:

### 1. PolyU (Hong Kong Polytechnic University) Palmprint Database
- Contains : 7752 grayscale palmprint images from 386 palms
- Format   : BMP images
- Source   : http://www.comp.polyu.edu.hk/~biometrics/

### 2. Delhi (IIT Delhi Touchless Palmprint Database v1.0)
- Contains : Left & right hand images from 230+ subjects
- Format   : BMP images
- Source   : http://web.iitd.ac.in/~ajaykr/Database_Palm.htm

### 3. CASIA (CASIA Palmprint Image Database)
- Contains : 5502 grayscale images from 312 subjects
- Format   : JPEG images
- Source   : http://www.cbsr.ia.ac.cn/PalmDatabase.htm
****************************************************************************************************
## Fun Fact
During parameter adjustment for Zhang's "Online Palmprint Identification" method, two specific figures
from the paper (Fig. 3 & Fig. 7) were found to have distinguishing palm marks (moles/scars). By analyzing
the PolyU database, these images were identified as PolyU_002_F_01 and PolyU_107_X_XX. Based on these
observations, initial parameters such as:

- Gaussian filter size
- Global threshold value
- ROI size
- Gabor filter size

were calibrated using these figures. Later, the parameters were further optimized for better performance.

These images and their relevant details are included in the repository with the same filenames as in the
PolyU database. They can be found in the `Docs/2003_Online_Palmprint_Identification/` directory.
****************************************************************************************************
### Final Thoughts
This project provides a comprehensive palmprint recognition pipeline, supporting two preprocessing methods
(morphological & tangent-based) and three major databases. The repository includes all necessary MATLAB scripts,
documentation, and references to ensure reproducibility and further research.

For any suggestions, feel free to contribute! 😊