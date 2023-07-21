function MIM_s = shuffle_MIM(data,npcs,fres, nshuf)
    %data is chan x l_epo x trials 
    % Copyright (c) 2022 Franziska Pellegrini and Stefan Haufe
    
    [nchan,l_epo,nepo]= size(data);
    
    [inds, PCA_inds] = fp_npcs2inds(npcs);
    ninds = length(inds);
    
    warning('One iteration takes about 90 seconds.')
    fprintf('Generating null distribution using %d shuffles...\n', nshuf)
    fprintf('Progress of %d:', nshuf);
    for ishuf = 1:nshuf %one iteration takes ~90 sec on my local laptop
        if mod(ishuf, 10) == 0
            fprintf('%d', ishuf);
        elseif mod(ishuf, 2) == 0
            fprintf('.');
        end
        
        %shuffle trials
        if ishuf == 1
            shuf_inds = 1:nepo; %true MIM values
        else
            shuf_inds = randperm(nepo);   
        end
        
        clear MIM2 CS COH2
        CS = fp_tsdata_to_cpsd(data, fres, 'WELCH', 1:nchan, 1:nchan,1:nepo,shuf_inds);
        
        for ifreq = 1:size(CS,3)
            clear pow
            pow = real(diag(CS(:,:,ifreq)));
            COH2(:,:,ifreq) = CS(:,:,ifreq)./ sqrt(pow*pow');
        end
        
        % loop over sender/receiver combinations to compute time-reversed GC
        for iind = 1:ninds
            if ~isequal(inds{iind}{1}, inds{iind}{2})
                %ind configuration
                subset = [inds{iind}{1} inds{iind}{2}];
                subinds = {1:length(inds{iind}{1}), length(inds{iind}{1}) + (1:length(inds{iind}{2}))};
                
                %MIC and MIM
                [~ , MIM2(:, iind)] =  roi_mim2(COH2(subset, subset, :), subinds{1}, subinds{2});
            end
        end        
          
        % extract measures out of the conn struct
        clear conn
        conn.MIM = MIM2;
        conn.inds = inds;  
        nroi = nchan/npcs(1);
        MIM_s(:, :, :, ishuf) = get_connect_mat(conn.MIM, nroi, +1);
        % [MIM_s(:,:,ishuf), ~, ~, ~, ~, ~] = fp_unwrap_conn(conn,nroi,filt1,PCA_inds);
    end
    fprintf('\n');
end

function measure = get_connect_mat( measureOri, nROI, signVal)
% create a ROI x ROI connectivity matrix, if needed
% TRGCmat(f, ii, jj) is net TRGC from jj to ii
    measure = [];
    iinds = 0;
    for iroi = 1:nROI
        for jroi = (iroi+1):nROI
            iinds = iinds + 1;
            measure(:, iroi, jroi) = signVal * measureOri(:, iinds);
            measure(:, jroi, iroi) = measureOri(:, iinds);
        end
    end
end
