% pop_toglleroiandica() - toggle between ROI and ICA processing
%
% EEG = pop_toglleroiandica(EEG, roiFlag)
%
% Inputs:
%   EEG     - CURRENT EEGLAB dataset
%   roiFlag - [0|1]
%
% Authors: Arnaud Delorme

% Copyright (C) Arnaud Delorme, arnodelorme@gmail.com
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
% THE POSSIBILITY OF SUCH DAMAGE.

function [STUDY, ALLEEG, com] = pop_toggleroiandica(STUDY, ALLEEG, roiFlag)

com = '';
if nargin < 1
    help pop_toggleroiandica;
    return;
end

if ~isfield(ALLEEG(1), 'roi') || ~isfield(ALLEEG(1).roi, 'source_roi_data')
    error('Cannot find ROI data - compute ROI data first');
end

if isfield(ALLEEG(1).roi, 'eeglab_using_roi')
    try
        curroiFlag = unique(cellfun(@(x)getfield(x, 'eeglab_using_roi'), { ALLEEG(:).roi }));
    catch
        curroiFlag = [0 1];
    end
    if length(unique(curroiFlag)) > 1
         warndlg2( strvcat( ...
            'Some datasets have different ROI settings.', ...
            'Restoring EEGLAB to use ICA. Call this', ...
            'menu again to use ROIs instead.'), 'Toggle ROI and ICA');
        [ALLEEG,com] = pop_toggleroiandica(ALLEEG, 0);
        return;
    end
else
    curroiFlag = false;
end

% check number of dims
% --------------------
dimsOK = true;
for iALLEEG = 1:length(ALLEEG)
    if isfield(ALLEEG, 'roi') && isfield(ALLEEG(1).roi, 'nPCA')
        if ALLEEG(1).roi.nPCA > 1
            error(strvcat('You must use only one PCA component when computing ROI', 'Check that all of your datasets are compliant'));
        end
    end
end

% check global workspace

if nargin < 2
    if curroiFlag % currently using ROI
        if ~isequal(ALLEEG(1).icaact, ALLEEG(1).roi.source_roi_data)
            res = questdlg2(strvcat('ROI data contained in EEG.icaact has been modified', ...
                'This could be a new ICA decomposition', 'Are you sure you want to erase it', ...
                'with the ICA solution stored by this function?'), 'Toggle ROI and ICA', ...
                'Cancel', 'Yes', 'Yes');
            if strcmpi(res, 'Cancel')
                return;
            end
        end
        res = questdlg2( strvcat( ...
            'Have EEGLAB use and plot ICA components instead of ROIs?'), 'Toggle ROI and ICA', ...
            'Cancel', 'Yes', 'Yes');
        if strcmpi(res, 'Yes')
            roiFlag = 0;
        else
            return; 
        end
    else
        res = questdlg2( strvcat( ...
            'Have EEGLAB use and plot ROIs instead of ICA components?', ...
            'EEGLAB can easily plot ICA components and can be made to plot', ...
            'ROI activity by temporarily using the ICA component matrices.', ...
            'ICA weights will be stored by this function and you may use', ...
            'this menu item again to restore them.'), 'Toggle ROI and ICA', ...
            'Cancel', 'Yes', 'Yes');
        if strcmpi(res, 'Yes')
            roiFlag = 1;
        else
            return; 
        end
    end
else
    if isequal(curroiFlag, roiFlag)
        if roiFlag
            disp('EEGLAB already using ROIs, not changing anything');
        else
            disp('EEGLAB already using ICA (not ROI), not changing anything');
        end
    end
end

disp('Updating EEGLAB menus...')
fig  = findobj('tag', 'EEGLAB');
changemenus( fig, roiFlag);

if roiFlag % use ROI
    disp('Exporting ROI activity into EEGLAB ICA matrices...');
    disp('We recommend to avoid saving datasets in that state');
    disp('and instead use this menu again after reloading the dataset(s).');
    for iDat = 1:length(ALLEEG)
        ALLEEG(iDat).roi.ica.icaweights = ALLEEG(iDat).icaweights;
        ALLEEG(iDat).roi.ica.icasphere = ALLEEG(iDat).icasphere;
        ALLEEG(iDat).roi.ica.icawinv   = ALLEEG(iDat).icawinv;
        ALLEEG(iDat).roi.eeglab_using_roi = true;
        if isfield(ALLEEG(iDat).roi, 'source_roi_data')
            ALLEEG(iDat).icaact = ALLEEG(iDat).roi.source_roi_data;
            ALLEEG(iDat).icasphere = eye(ALLEEG(iDat).nbchan,ALLEEG(iDat).nbchan);
            ALLEEG(iDat).icachaninds = 1:ALLEEG(iDat).nbchan;
            ALLEEG(iDat).icaweights = ones(size(ALLEEG(iDat).icaact, 1), ALLEEG(iDat).nbchan);
            ALLEEG(iDat).icawinv = pinv(ALLEEG(iDat).icaweights * ALLEEG(iDat).icasphere);
            ALLEEG(iDat).icalabels = { ALLEEG(iDat).roi.atlas.Scouts.Label };
        else
            ALLEEG(iDat).icaact = [];
            ALLEEG(iDat).icasphere = [];
            ALLEEG(iDat).icachaninds = [];
            ALLEEG(iDat).icaweights = [];
            ALLEEG(iDat).icawinv = [];
            ALLEEG(iDat).icalabels = { };
        end
    end
    if ~isempty(STUDY)
        STUDY.etc.icacluster = STUDY.cluster;
        
        % rebuild STUDY cluster (use the same subjects with different runs)
        labels = ALLEEG(1).icalabels;
        STUDY.cluster = [];
        STUDY.cluster(1).name = 'Parentcluster 1';
        STUDY.cluster(1).sets  = repmat([1:length(STUDY.datasetinfo)]', length(labels));
        STUDY.cluster(1).comps = 1:length(labels);
        STUDY.cluster(1).child  = {};
        STUDY.cluster(1).parent = {};

        for iLab  = 1:length(labels)
            for iSubj = 1:length(STUDY.subject)
                indSubj = strmatch(STUDY.subject{iSubj}, { STUDY.datasetinfo.subject}, 'exact')';
                STUDY.cluster(iLab+1).name  = labels{iLab};
                STUDY.cluster(iLab+1).sets  = [1:length(STUDY.datasetinfo)]';
                STUDY.cluster(iLab+1).comps = iLab; %repmat(iLab, 1, size(STUDY.cluster(1).sets,2));
                STUDY.cluster(iLab+1).child  = {};
                STUDY.cluster(iLab+1).parent = {};
                STUDY.cluster(1).child{end+1} = labels{iLab};
            end
        end
    end
    set(findobj(fig, 'tag', 'toggleroi'), 'Text', 'Have EEGLAB use ICA instead of ROIs');
else
    disp('Putting stored ICA activity (if any) into EEGLAB ICA matrices...');
    for iDat = 1:length(ALLEEG)
        ALLEEG(iDat).icaweights = ALLEEG(iDat).roi.ica.icaweights;
        ALLEEG(iDat).icasphere  = ALLEEG(iDat).roi.ica.icasphere;
        ALLEEG(iDat).icawinv    = ALLEEG(iDat).roi.ica.icawinv;
        ALLEEG(iDat).icaact     = [];
        ALLEEG(iDat).icalabels  = {};
        ALLEEG(iDat).roi.eeglab_using_roi = false;
    end
    if ~isempty(STUDY)
        STUDY.clusters = STUDY.etc.icaclusters;
    end
    set(findobj(fig, 'tag', 'toggleroi'), 'Text', 'Have EEGLAB use ROIs instead of ICA');
end

com = sprintf('ALLEEG = pop_toggleroiandica(ALLEEG, %d); EEG = ALLEEG(CURRENTSET);', roiFlag);

% -------------
% replace menus
% -------------
function changemenus( fig, roiFlag)
if isempty(fig), return; end;
allmenus = findobj(gcf, 'type', 'uimenu');

replaceText = { 'component', 'ROI';
                'Component', 'ROI';
                'cluster', 'accross participant' };
% disableIfText = { 'ica', 'component' };
if ~roiFlag
    replaceText = replaceText(:, [2 1]);
end
    
for iMenu = 1:length(allmenus)
    userDat = get(allmenus(iMenu), 'userdata');
    label = get(allmenus(iMenu), 'Text');
%     for iText = 1:length(disableIfText)
%         if ~isempty(strfind(label, disableIfText{iText}))
%             set(allmenus(iMenu), 'Enable', 'off');
%         end
%     end
    if isempty(strfind(userDat, 'roi:off')) % contains only introduced in Matlab 2016
        label = get(allmenus(iMenu), 'Text');
        for iRep = 1:size(replaceText,1)
            label = strrep(label, replaceText{iRep,1}, replaceText{iRep,2});
            if strcmpi(replaceText{iRep,2}, 'component')
                if length(label) > length(replaceText{iRep,2})
                    if isequal(label(1:9), 'component')
                        label(1) = 'C';
                    end
                end
            end
        end
        set(allmenus(iMenu), 'Text', label);
    end
end

