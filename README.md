[![GitHub stars](https://img.shields.io/github/stars/arnodelorme/roiconnect?color=green&logo=GitHub)](https://github.com/arnodelorme/roiconnect/issues) [![GitHub issues](https://img.shields.io/github/issues/arnodelorme/roiconnect?color=%23fa251e)](https://github.com/arnodelorme/roiconnect/issues) [![GitHub pulls](https://img.shields.io/github/issues-pr-raw/arnodelorme/roiconnect)](https://github.com/arnodelorme/roiconnect/issues-pr-raw) [![GitHub forks](https://img.shields.io/github/forks/arnodelorme/roiconnect?style=social)](https://github.com/arnodelorme/roiconnect/forks) [![GitHub contributors](https://img.shields.io/github/contributors/arnodelorme/roiconnect?style=social)](https://github.com/arnodelorme/roiconnect/contributors)

# What is ROIconnect?

ROIconnect is a freely available open-source plugin to [EEGLAB](https://github.com/sccn/eeglab) for EEG data analysis. It allows you to perform functional connectivity analysis between and within regions of interests (ROIs) on source level.  The results can be visualized in 2-D and 3-D. ROIs are defined based on popular fMRI atlases, and source localization is performed through LCMV beamforming and eLORETA. Connectivity analysis is performed between all pairs of brain regions using Granger Causality, Time-reversed Granger Causality, Multivariate Interaction Measure, Maximized Imaginary Coherency, and other methods. Additionally, this plugin is compatible with Fieldtrip, Brainstorm and NFT head models.

📚 Check out the following paper to learn about recommended methods and pipelines for connectivity experiments:
> Pellegrini, F., Delorme, A., Nikulin, V. & Haufe, S., 2022. Identifying best practices for detecting inter-regional functional connectivity from EEG. <preprint_location>.

You can choose to access the core functions from the EEGLAB GUI. Experienced users can access additional utilities from the command line. If you do decide to run a function from the command line, please refer to the respective documentation provided in the code. 


Code developed by Stefan Haufe and team, with EEGLAB interface, coregistration, 3-D vizualisation and Fieldtrip integration performed by Arnaud Delorme.

# Installation
First of all, please make sure to have [EEGLAB](https://github.com/sccn/eeglab#installingcloning) installed. Your project layout should look as follows

```
eeglab/	 				
├── functions/ 							
├── plugins/ 								
├── sample_data/ 					
├── sample_locs/		
├── tutorial_scripts/						
└── /.../ 					
```
The ROIconnect plugin should be installed in the `plugins` folder. The easiest way to do so is to simply clone our repository. Navigate to the `plugins` folder by typing 
```
cd <eeglab_install_location>/eeglab/plugins
``` 
Next, clone the repository
```
git clone https://github.com/arnodelorme/roiconnect.git
```
That's it! If you want to run the plugin, please start [EEGLAB](https://github.com/sccn/eeglab#to-use-eeglab) first. You may need to add EEGLAB to the [MATLAB path](https://de.mathworks.com/help/matlab/ref/addpath.html).  

# Key features
The features of the toolbox are implemented in the following three main functions: `pop_roi_activity`, `pop_roi_connect` and `pop_roi_connectplot`.

## Source reconstruction
`pop_roi_activity` asks for the following inputs: an EEG struct containing EEG sensor activitiy, a pointer to the headmodel and a source model, the atlas name and the number of PCs for the dimensionality reduction. In addition, this function also supports [FOOOF analysis](https://fooof-tools.github.io/fooof/). Here is a command line example including FOOOF:

```matlab
EEG = pop_roi_activity(EEG, 'leadfield',EEG.dipfit.sourcemodel,'model','LCMV','modelparams',{0.05},'atlas','LORETA-Talairach-BAs','nPCA',3, 'fooof', 'on', 'fooof_frange', [1 30]);
```

The function performs source reconstruction by calculating a source projection filter and applying it to the sensor data. Power is calculated with the Welch method on the voxel time series and then summed across voxels within regions. To enable region-wise FC computation, the function applies PCA to the time series of every region. It then selects the *n* strongest PCs for every region. The resulting time series is stored in `EEG.roi.source_roi_data`, and power is stored in `EEG.roi.source_roi_power`.

## Connectivity analysis
`pop_roi_connect` accepts the following inputs: the EEG struct computed by `pop_roi_activity` and the names of the FC metrics. In addition, you can tell the function to estimate FC on time snippets of 60 s length (default) which can be averaged (default) or used as input for later statistical analyses. The following command line example asks the function to perform FC analysis on snippets using default values (explicitely passed as input parameters in this example). 

```matlab
EEG = pop_roi_connect(EEG, 'methods', { 'MIM', 'TRGC'}, 'snippet', 'on', 'snip_length', 60, 'fcsave_format', 'mean_snips');
```

The function computes all FC metrics in a frequency-resolved way, i.e., the output contains FC scores for every region-region frequency combination. We include the option to estimate FC on time snippets to avoid biases due to different data lengths. The output of this function is stored in `EEG.roi.<fc_metric_name>`.

## Visualization
You can visualize power and FC in different modes by calling `pop_roi_connectplot`. Below, we show results of one subject from the real data example in [[1]](#1). You can find the MATLAB code and corresponding analyes [here](https://github.com/fpellegrini/MotorImag). The left plots show power or FC in left motor imagery, the right plots indicate results in right motor imagery. Due to the nature of the task, we show results in the alpha frequency band but you are free to choose any frequency or frequency band you want. 

:pushpin: If any of the images are too small for you, simply click on them, they will open in full size in another tab.

### Power as a region-wise bar plot
If you wish to visualize power as a barplot only, please make sure to explicitely turn `plotcortex` off because it is turned on by default. 
```matlab
EEG = pop_roi_connectplot(EEG, 'measure', 'roipsd', 'plotcortex', 'off', 'plotbarplot', on, 'freqrange', [8 13]) % alpha band;
```
<p float="middle">
  <img src="https://github.com/Hiyeri/roiconnect/blob/master/resources/power_barplot_left.jpg?raw=true" width="400"/>     
  &nbsp; &nbsp;
  <img src="https://github.com/Hiyeri/roiconnect/blob/master/resources/power_barplot_right.jpg?raw=true" width="400"/>
</p>

### Power as a source-level cortical surface topography
```matlab
EEG = pop_roi_connectplot(EEG, 'measure', 'roipsd', 'plotcortex', 'on', 'freqrange', [8 13]);
```
<p float="middle">
  <img src="https://github.com/Hiyeri/roiconnect/blob/master/resources/power_cortex_left.jpg?raw=true" width="400"/>     
  &nbsp; &nbsp;
  <img src="https://github.com/Hiyeri/roiconnect/blob/master/resources/power_cortex_right.jpg?raw=true" width="400"/>
</p>

### FC as region-to-region matrix 
Again, if you do not wish to see the cortex plot, you should explicitely turn `plotcortex` off. 
```matlab
EEG = pop_roi_connectplot(EEG, 'measure', 'mim', 'plotcortex', 'off', 'plotmatrix', 'on', 'freqrange', [8 13]);
```
<p float="middle">
  <img src="https://github.com/Hiyeri/roiconnect/blob/master/resources/FC_MIM_matrix_left.jpg?raw=true" width="400"/>     
  &nbsp; &nbsp;
  <img src="https://github.com/Hiyeri/roiconnect/blob/master/resources/FC_MIM_matrix_right.jpg?raw=true" width="400"/>
</p>

You can additionally filter by hemispheres and regions belonging to specific brain lobes. As an example, let us see how FC of the left/right hemisphere look like for left/right motor imagery.

```matlab
pop_roi_connectplot(EEG, 'measure', 'mim', 'plotcortex', 'off', 'plotmatrix', 'on', 'freqrange', [8 13], 'hemisphere', 'left') % left hemisphere, left motor imagery;
```
<p float="middle">
  <img src="https://github.com/Hiyeri/roiconnect/blob/master/resources/FC_MIM_matrix_lefthem_left.jpg?raw=true" width="400"/>     
  &nbsp; &nbsp;
  <img src="https://github.com/Hiyeri/roiconnect/blob/master/resources/FC_MIM_matrix_righthem_right.jpg?raw=true" width="400"/>
</p>

### Net FC as a cortical surface topography
Here, the mean FC from all regions to all regions is visualized.
```matlab
pop_roi_connectplot(EEG, 'measure', 'mim', 'plotcortex', 'on', 'freqrange', [8 13]);  
```
<p float="middle">
  <img src="https://github.com/Hiyeri/roiconnect/blob/master/resources/FC_MIM_cortex_left.jpg?raw=true" width="400"/>     
  &nbsp; &nbsp;
  <img src="https://github.com/Hiyeri/roiconnect/blob/master/resources/FC_MIM_cortex_right.jpg?raw=true" width="400"/>
</p>

### Seed FC as a cortical surface topography
Here, the FC of a seed region to all other regions is visualized.
```matlab
pop_roi_connectplot(EEG, 'measure', 'mim', 'plotcortex', 'on', 'freqrange', [8 13], 'plotcortexseedregion', 49); 
```
<p float="middle">
  <img src="https://github.com/Hiyeri/roiconnect/blob/master/resources/FC_MIM_cortex_seed49_left.jpg?raw=true" width="400"/>     
  &nbsp; &nbsp;
  <img src="https://github.com/Hiyeri/roiconnect/blob/master/resources/FC_MIM_cortex_seed50_right.jpg?raw=true" width="400"/>
</p>


# References
<a id="1">[1]</a> 
Pellegrini, F., Delorme, A., Nikulin, V. & Haufe, S. (2022). 
Identifying best practices for detecting inter-regional functional connectivity from EEG. 
<preprint_location>.

<a id="2">[2]</a> 
https://github.com/fpellegrini/FCsim

<a id="3">[3]</a> 
https://github.com/fpellegrini/MotorImag

