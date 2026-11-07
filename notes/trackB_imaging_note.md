# Track B: biofilm image analysis (practice)

- Dataset: BiofilmQ sample data (Sample_data_v1.0.0), file Pos0-1_ch1_frame000001_Nz21.tif. Source and citation: see notes/trackB_data_source.txt.
- Software and steps: Fiji. Calibration from the dataset's metadata file (pixel size 0.063 um, slice spacing 0.4 um). Otsu threshold (16-bit histogram, dark background) calculated on slice 1, lower limit 8358, applied to all 22 slices. Thresholded area measured per slice; biovolume = sum of area x 0.4 um.
- What worked: One stack segmented and measured. Approximate biovolume 45.2 um^3. Saved results/trackB_measurements.csv, figures/trackB_projection.png and figures/trackB_area_vs_depth.png.
- What failed: The first threshold (2778) selected the whole image until I switched on the 16-bit histogram. Stack histogram and 16-bit histogram could not be used together, so the threshold came from slice 1 only.
- Limitations: One image at one time point, practice data (not my own experiment). Threshold was set from slice 1 and checked by eye, so it may miss dimmer cells in other slices. Most of the thresholded area (about 92 of 113 um^2) is in slice 1, and I have not tested why. Biovolume from a single threshold is approximate.
- Next to learn: Checking 3D segmentation quality, and comparing with BiofilmQ's own segmentation.