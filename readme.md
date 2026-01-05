# About

This repo analyze the loci tracking data by combining the output from oufti & u-track. It calculates the normalized position of foci, and does diffusion analysis for tracks.

All images are taken using microscope B, PC lens (NA=1.45, 100x), Andor camera (pixel=160nm)

- imaging setting: YFP_yh, 5mW 514nm, 20-200ms, 101 frames, Gain 3, 300x, angle=8400

- Current dataset: SK734, 727, 729, 731 (2025 Fall, taken by Liam)



### Before Running the Code

1. change current working directory to the folder that contains `tracking00x` folders (uTrack output) and `mesha.mat` (oufti output)

2. change `varPath` to your own local folder to store analysis results (all results will be saved under the `lociPath`)

3. make sure the oufti output file `mesh.mat` is cleaned up and saved as `mesha.mat`



#### Major scripts

- `lociAnalysis_spotNorm.m`: combine oufti & uTrack, run spotNorm, diffusion analysis, save tf & Loci files
- `plot_spotNorm.m`: plot xNorm & lNorm from the data



#### Raw data

- time-lapse fluorescent images of gene loci

- phase contrast images of cells for each FOV



#### u-track analysis parameters

- std = 1 pix, alpha = 0.01, alpha = 0.01
- frame 0 gap, 40+ frame
- 20ms exposure, 200ms interval, search radius = 2 pix



#### oufti analysis parameter

- `Phase_SPT_alvin.set`, post processed by `meshCleanup_andor.m`

