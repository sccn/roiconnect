function fp_fc_shuffletest(isb)
%calculates true MIM and generates MIM null distribution 
%isb = subject index 
% Copyright (c) 2022 Franziska Pellegrini and Stefan Haufe

DIRIN = './Data_MI/';

%subjects with high performance classification
subs = [3 4 5 8 9 11 12 14 15 16 17 18 19 21 22 23 25 27 28 29 30 31 33 34 35 37];
nshuf = 1001; %first shuffle = true MIM, then 1000 shuffles

%% load preprocessed EEG
sub = ['vp' num2str(subs(isb))];
load([DIRIN sub '.mat'])

%% shuffle
%one shuffle ~ 1 min on local machine
%first shuffle is true MIM 
npcs = repmat(3,1,nroi);
MIM_s = fp_shuffle_MIM(data,npcs,fs,filt1, nshuf); 

%% save 
save([DIRIN 'RDE_shuf_' num2str(isb) '.mat'],'-v7.3')
