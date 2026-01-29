## About

This repo analyze the loci tracking data by combining the output from oufti & u-track. It calculates the normalized position (spotNorm) of foci, and does diffusion analysis for tracks.

All images are taken using microscope B, PC lens (NA=1.45, 100x), Andor camera (pixel=160nm)

- imaging setting: YFP_yh, 5mW 514nm, 20-200ms, 101 frames, Gain 3, 300x, angle=8400
- Current dataset: SK734, 727, 729, 731 (2025 Fall, taken by Liam)

### Major scripts

- `lociPrepare_spotNorm.m`: prepares the raw loci tracking data to be ready for uTrack/oufti analysis. 
- `lociAnalysis_spotNorm.m`: combines oufti & uTrack, run spotNorm, diffusion analysis, save `tf oufti` & `Loci oufti` files.
  - this scripts call function: `spotNorm_loci`, `reNumCells_yh`, `setCameraGain`, `getCellInfo`.
- `lociCombine_spotNorm.m`: combines multiple single-day analysis files saved by `lociAnalysis_spotNorm.m`.
- `meshCleanup_andor.m`: cleans up the cell meshes with abnormal cell size and shape.
- `plot_spotNorm.m`: plots xNorm & lNorm from the data.
- `plot_MSD_fit.m`: plots MSD & fitting of the tracking data.
- `plot_all.m`: plot basic properties of the data: track length, MSD, bleaching, D & alpha & locErr histogram.

### Minor scripts

- `plot_subMSD_cell.m`: plots MSD of subpopulations binned by cell geometry (length/width).
- `plot_subMSD_spotNorm.m`: plots MSD of subpopulations binned by spotNorm (xNorm/LNorm).
- `plot_subMSD_amp.m`: plots MSD of subpopulations binned by signal intensity (amplitude).
- `plot_MSDxy_cell.m`: plots MSD along cellular long & short axes, two plots with multiple dataset overlaid
- `plot_MSDxy_compare.m`: plots MSD along cellular long & short axes, one plot for each dataset
- `plot_corr_tracks.m`: plots correlation of any two track quantities using scatter plot.
- `plot_MSD_truncate_fit.m`: plot TA-MSD using truncated tracks.


## Before Running the Code

1. Change current working directory to the folder that contains `tracking00x` folders (uTrack output) and `mesha.mat` (oufti output).
2. Always change `varPath` in the code to your own local folder to store analysis results (all results will be saved under the `lociPath`).
3. Make sure the oufti output file `mesh.mat` is cleaned up and saved as `mesha.mat`.

## Running workflow

1. When capturing images using the NIS Element, follow naming rules
   1. folder name: `date-strain extra` (i.e. `260101-SK1 Gly`, `extra` is optional)
   2. tracking file: `tracking0xx.nd2` (xx: 001-099)
   3. phase contrast images: `epi0xx.nd2` 
2. Export ND2 files to Tiff files in NIS Element (follow instruction)
   1. all phase contrast images should be exported into one folder (same as input)
   2. all tracking images should be exported into individual subfolder
3. Transfer the data to the hard drive
4. Run `lociPrepare_spotNorm.m` to prepare files for oufti & u-track analysis
5. Run oufti analysis using the parameter file `phase_SPT_alvin.set`
   1. save the output file `mesh.mat` inside the data folder 
6. Run u-track analysis using the parameters below
   1. save the output file the same as the input folder (tracking00x)
7. Run `lociAnalysis_spotNorm.m` to combine oufti & uTrack results, and for further spotNorm & diffusion analysis (results will be saved under the `lociPath`)
8. Run plotting scripts to visualize the analysis results
9. (optional) To combine multiple single-day results, run `lociCombine_spotNorm`. It will combined single-day `tf oufti` & `Loci oufti` files and move the individual files to the `single day` subfolder. 

## Analysis Parameters

### u-track analysis

- std = 1 pix, alpha = 0.01, alpha = 0.01
- frame 0 gap, 40+ frame
- 20ms exposure, 200ms interval, search radius = 2 pix

### oufti analysis

- `Phase_SPT_alvin.set`, post processed by `meshCleanup_andor.m`

## Folder structure

This is the local folder created by `lociPrepare_spotNorm.m`. the folder is used to save analysis results and is named according to the hard drive folder containing raw ND2 data.

```text
251108-SK731/
├── phase/
│   ├── epi001.tif
│   ├── ...
│   └── epi00x.tif
├── tracking001/
├── ...
├── tracking00x/
```

This is the folder structure after running oufti & u-track analysis. oufti results are saved under the folder as `mesh.mat`. u-track results are saved inside each `tracking00x` folder. `mesha.mat` are saved by `meshCleanup_andor.m` based on `mesh.mat`.

```text
251108-SK731/
├── phase/
│   ├── epi001.tif
│   ├── ...
│   └── epi00x.tif
├── tracking001/
│   ├── backups/
│   ├── TrackingPackage/
|       ├── GaussianMixtureModels/
|           └── Channel_1_detection_result.mat
│       └── tracks/
|           └── Channel_1_tracking_result.mat
|   └── movieData.mat
├── ...
├── tracking00x/
├── mesh.mat
└── mesha.mat
```