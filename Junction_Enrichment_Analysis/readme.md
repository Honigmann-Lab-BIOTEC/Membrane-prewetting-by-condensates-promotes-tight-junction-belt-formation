Programming language: MATLAB R2019a, Toolbox Image Processing required

OS Name: Windows 10 Enterprise LTSC

OS Type: 64-bit

Analysis code for evaluation protein enrichment in ZO1 condensates from 2-color live cell movies of mN-ZO1 and mS-Client protein after Ca-switch (see Fig. S2)

Make sure to add the folder "Matlab helper functions" to your MATLAB path.

The code quantifies protein enrichment kinetics for mN-ZO1 and mS-Client from the input data in the folder "Data" and stores the results in folder "Analysis Results".
Inputs (Run Code and Select Folder "Data"): 
- Maximum intensity projections over time of the ZO1 channel.
- Maximum intensity projections over time of the Client channel.
- Segmentation of the ZO1 channel (done with FIJI, Rolling Ball Background removal Size 3, Binarize)
- Cell Segmentation of the ZO1 Channel (done with CellPose, Cyto2 model)

Processing:
The code performs the following steps.
- Loading the data
- Filtering the ZO1 segmentation to remove pixel noise, i.e. segmented regions should be larger than 5 pixels.
- For each timepoint all segmented cells are evaluated sequentially.
  - The corresponding segmented ZO1 condensates for each cell are found by dilating the cell region and getting the overlap with ZO1 segmentation.
  - The signal intensities in the overlapping ZO1 condensates are measured for both channels
  - The signal intensity of the cell cytoplasm is measured from the cell region (excluding the ZO1 segmentation)
  - The junction enrichment factor is calculated for each cell mean(intensity_ZO1_condensates) / mean(intensity_cytoplasm)
- For each timepoint the mean enrichment factors for all cells are calculated and normalized.
- The data and plots are stored in folder "Analysis Results".


