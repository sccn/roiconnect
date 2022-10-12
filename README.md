[![GitHub stars](https://img.shields.io/github/stars/arnodelorme/roiconnect?color=green&logo=GitHub)](https://github.com/arnodelorme/roiconnect/issues) [![GitHub issues](https://img.shields.io/github/issues/arnodelorme/roiconnect?color=%23fa251e)](https://github.com/arnodelorme/roiconnect/issues) [![GitHub pulls](https://img.shields.io/github/issues-pr-raw/arnodelorme/roiconnect)](https://github.com/arnodelorme/roiconnect/issues-pr-raw) [![GitHub forks](https://img.shields.io/github/forks/arnodelorme/roiconnect?style=social)](https://github.com/arnodelorme/roiconnect/forks) [![GitHub contributors](https://img.shields.io/github/contributors/arnodelorme/roiconnect?style=social)](https://github.com/arnodelorme/roiconnect/contributors)

# What is ROIconnect?

ROIconnect is a freely available open-source plugin to [EEGLAB](https://github.com/sccn/eeglab) for EEG data analysis. It allows you to perform functional connectivity analysis between and within regions of interests (ROIs) on source level.  The results can be visualized in 2-D and 3-D. ROIs are defined based on popular fMRI atlases, and source localization is performed through LCMV beamforming and eLORETA. Connectivity analysis is performed between all pairs of brain regions using Granger Causality, Directed Transfer Entropy, and many other methods. 

You can choose to access the core functions from the EEGLAB GUI. Experienced users can access additional utilities from the command line. If you do decide to run a function from the command line, please refer to the respective documentation provided in the code. 


Code developed by Stefan Haufe and team, with EEGLAB interface, coregistration, 3-D vizualisation and Fieldtrip integration performed by Arnaud Delorme.

# Installation

# Key features
The features of the toolbox are implemented in the following three main functions: `pop_roi_activity`, `pop_roi_connect` and `pop_roi_connectplot`.

## Source reconstruction
`pop_roi_activity` asks for the following inputs: an EEG struct containing EEG sensor activitiy, a pointer to the headmodel and a source model, the atlas name and the number of PCs for the dimensionality reduction. In addition, this function also supports [FOOOF analysis](https://fooof-tools.github.io/fooof/). Here is a command line example including FOOOF:

```matlab
EEG = pop_roi_activity(EEG, 'leadfield',EEG.dipfit.sourcemodel,'model','LCMV','modelparams',{0.05},'atlas','LORETA-Talairach-BAs','nPCA',3, 'fooof', 'on', 'fooof_frange', [1 30]);
```

The function performs source reconstruction by calculating a source projection filter and applying it to the sensor data. Power is calculated with the Welch method on the voxel time series and then summed across voxels within regions. To enable region-wise FC computation, the function applies PCA to the time series of every region. It then selects the *n* strongest PCs for every region. The resulting time series is stored in `EEG.roi.source_roi_data`, power is stored in `EEG.roi.source_roi_power`.

# Available features

- ROI connectivity analysis using TRGC and CS
- Use of surface or volumetric atlases
- Compatibile with Fieldtrip, Brainstorm and NFT head models


