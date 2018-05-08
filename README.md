Dynamic Connectivity States (DCS) identified from k-means clustering (k=3) of dynamic functional connectivity (DFC) windows computed from resting-state data across 52 adolescents (age range = 15-19y)
=====================================================================================

Files description
1) DFC-centroids-k3.mat  :  MATLAB data file with cluster centroids information
	C1 = centroid 1 (intermediate state)
	C2 = centroid 2 (Low Arousal State DCS)
	C3 = centroid 3 (High Arousal State DCS) 
Centroid is arranged in a 122x122 matrix, where C(x,y) represent connectivity strength between node x and y.
The list of ROIs is provided in ROIslist.txt in the same order, i.e. C(10,2) is the connectivity between the 10th ROI (lh.17Networks_LH_DorsAttnB_FEF) and 2nd ROI (lh.17Networks_LH_VisCent_ExStr)

2) ROIslist.txt  :  the list of ROIs included in DFC computation
There are 122 ROIs in total, comprises of 114 ROIs in the cortical areas (data in surface space) and 8 ROIs in the subcortical regions. 
114 cortical ROIs were adapted from the 17-network parcellation described in Yeo et al., 2011. 
Subcortical ROIs (including left and right hemispherical thalamus, striatum, hippocampus, and amygdala) were extracted from Freesurfer segmentation of the FSL MNI152 template. 


=====================================================================================
Visualization of correlation matrix
Run MATLAB script PlotCorrelationMatrix.m to plot the correlation matrix.
Example:
>> load DFC-centroids-k3.mat
>> PlotCorrelationMatrix(C1, [-0.6 0.6], 'C1')
This will plot the connectivity pattern of the first centroid and save it with the name "C1_minsc-0.60_maxsc0.60.jpg"
Read the script's help function for more details.
The final arrangement of ROIs is listed in ROIs_reordered.txt.


=====================================================================================
Scan parameters
Participants underwent two runs of a 6-min resting state (eyes open, black screen) with the following scan parameters: TR = 2000ms, TE = 30 ms, FA = 90 degree, FoV = 192x192 mm, isotropic voxel size = 3x3x3 mm. 
High resolution structural images were also collected using MPRAGE sequence (TR = 2300 ms, TI = 900 ms, FA = 8 degree, voxel dimension = 1x1x1 mm, FoV = 256 x 240 mm)


=====================================================================================
fMRI preprocessing
1) Discard the first four frames of each run
2) Correct for slice acquisition-dependent time shifts in each volume with SPM
3) Correcting for head motion using rigid body translation and rotation parameters with FSL
4) Remove linear trends of each run and apply low pass filter to retain frequency below 0.08 Hz
5) Remove spurious variance from head motion, whole brain signal, ventricle signal, white matter signal, and their derivatives using linear regression
6a) Project functional data onto the Freesurfer surface space (2 mm mesh), smooth with a 6 mm FWHM kernel, and then downsample to a 4 mm surface mesh
6b) Normalize functional data to MNI152 template using nonlinear registration in Freesurfer, downsample to 2 mm voxels, and then smooth with a 6 mm FWHM kernel
(Structural processing was done using Freesurfer 4.5.0)


=====================================================================================
Dynamic functional connectivity (DFC)
DFC was computed using a sliding window approach following Allen, et al. Specifically, BOLD time series from the 122 ROIs (see list of ROIs in ROIs_list.txt) were first de-spiked and de-meaned.  A tapered window was constructed by convolving a rectangular window (20 TRs) and a Gaussian function (sigma = 3 TRs). Covariance among all possible ROIs pairs within the tapered window were estimated using the regularized precision matrix. The graphical LASSO method with L1 norm penalty (regularization parameter  = 0.1) was applied to promote sparsity. This process was repeated by shifting the tapered window by 1 TR. For each functional run, we obtained 156 covariance matrices, each with 7381 (122 x 121 / 2) unique correlation values.
Covariance matrices from 52 subjects were collapsed together and k-means clustering was performed to classify each DFC matrix using L1 distance as the cost function. 

