Programming language: MATLAB R2019a, Toolbox Image Processing required
OS Name: Windows 10 Enterprise LTSC
OS Type: 64-bit

Analysis code for STED images of ZO1 and PATJ at tight junctions

Make sure to add the folder "Matlab helper functions" to your MATLAB path.

The code produces a distance histogram of the potein PATJ to the most apical strand of the protein ZO1 (See figure 3 in the paper).

1. Use "Run_Distance_Analysis.m" to analyze STED images (.tif files) in the folder "Processed Images". Images are rotated such that the apical side is at the top. The code first segments and skeletonizes the ZO1 channel. Second the PATJ channel is segmented and the centers of PATJ clusters are identified. If necessary tweak segmentation input parameters to optain ideal segmentation. Finally, the nearest neighbor distance of each PATJ cluster to the most apical ZO1 strand is determined. The code saves the distance results (-Dist_Analysis.mat) for each image into the folder "Analysis Results".

2. Use "Pool_distance_measurements.m" to open all analyzed mat files in the folder "Analysis Results" and plot a histogram of the distances of PATJ to the apical interface. The code saves the pooled distance measurements and the Histogram-Figure into the same folder.
    
