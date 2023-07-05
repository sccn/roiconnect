function fp_plot_fc_shuffletest
%Function that generates p-values from the MIM null distribution and plots
%them 
% Copyright (c) 2022 Franziska Pellegrini and Stefan Haufe

DIRIN = '~/Dropbox/Franziska/Data_MEG_Project/RDE_shuffletest/right_MI/';
DIRFIG = '~/Desktop/';

nsub = 26;

%% generate p-values by comparing true MIM to null distribution 
for isb = 1:nsub    
    in = load([DIRIN 'RDE_shuf_' num2str(isb) '.mat']);
    
    %average over one region dimension to obtain netMIM 
    MIM_s = squeeze(mean(in.MIM_s,2)); 
    MIM_pn(:,isb) = sum(MIM_s(:,1)< MIM_s(:,2:end),2)./(size(in.MIM_s,3)-1);   
end

%% Use Stouffer's method to aggregate p-values 
nroi = size(MIM_pn,1);
for iroi = 1:nroi
    MIM_pn_s(iroi) = fp_stouffer(squeeze(MIM_pn(iroi,:)));
end

%% Use FDR-correction for multiple comparison's correction 
[p_fdr, ~] = fdr(MIM_pn_s, 0.05);
MIM_pn_s(MIM_pn_s>p_fdr)=1;

%% Plot 
load cm17;
load('bs_results.mat'); % load cortex 
MIM_pn_s(MIM_pn_s==0) = 0.001; % 1/nshuf
data = -log10(MIM_pn_s);
allplots_cortex_BS(cortex_highres,data, [35 52], cm17a ,'-log(p)', 0.3,...
    [DIRFIG 'netMIM_right'])