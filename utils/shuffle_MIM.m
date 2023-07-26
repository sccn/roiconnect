function conn = shuffle_MIM(data, npcs, output, fres, nshuf)
    % TO DO: 
    %   - add proper documentation
    %   - add g.freqresolution

    %data is chan x l_epo x trials 
    % Copyright (c) 2022 Franziska Pellegrini and Stefan Haufe
    
    [nchan, ndat, nepo] = size(data);

    [inds, PCA_inds] = fp_npcs2inds(npcs);
    ninds = length(inds);

    CSpara = [];
    CSpara.subave = 0;
    CSpara.mywindow = hanning(ndat) ./ sqrt(hanning(ndat)' * hanning(ndat));
%     CSpara.freqresolution = g.freqresolution;
    CSpara.freqresolution = 0; % for now
    
%     warning('One iteration takes about 90 seconds.')
    fprintf('Generating null distribution using %d shuffles...\n', nshuf)
    fprintf('Progress of %d:', nshuf);
    for ishuf = 1:nshuf 
        if mod(ishuf, 10) == 0
            fprintf('%d', ishuf);
        elseif mod(ishuf, 2) == 0
            fprintf('.');
        end
        
        % shuffle trials
        if ishuf == 1
            shuf_inds = 1:nepo; % true MIM values
        else
            shuf_inds = randperm(nepo);   
        end
        
        clear MIM2 CS

        data_shuf = data(:, :, shuf_inds);
        [CS, ~, wPLI, ~] = data2cs_event_shuf(data(:, :)', data_shuf(:, :)', ndat, floor(ndat/2), ndat, [], CSpara);
%         CS = fp_tsdata_to_cpsd(data, fres, 'WELCH', 1:nchan, 1:nchan,1:nepo,shuf_inds);
        
        if ~isempty(intersect(output, {'MIM', 'MIC', 'cCOH' 'iCOH', 'aCOH'}))
            clear cCOH iCOH aCOH
            for ifreq = 1:size(CS,3)
                clear pow 
                pow = real(diag(CS(:,:,ifreq)));
                cCOH(:,:,ifreq) = CS(:,:,ifreq) ./ sqrt(pow*pow');
                aCOH(:,:,ifreq) = abs(CS(:,:,ifreq)) ./ sqrt(pow*pow');
                iCOH(:,:,ifreq) = abs(imag(CS(:,:,ifreq) ./ sqrt(pow*pow')));
            end
        end
        
        % loop over sender/receiver combinations to compute time-reversed GC
        for iind = 1:ninds
            if ~isequal(inds{iind}{1}, inds{iind}{2})
                %ind configuration
                subset = [inds{iind}{1} inds{iind}{2}];
                subinds = {1:length(inds{iind}{1}), length(inds{iind}{1}) + (1:length(inds{iind}{2}))};
                
                %MIC and MIM
                [~ , MIM2(:, iind)] =  roi_mim2(cCOH(subset, subset, :), subinds{1}, subinds{2});
            end
        end         
        nroi = nchan/npcs(1);
        
        % reshape (in the case of MIM) or only keep the first principal component (other metrics)
        MIM_s(:, :, :, ishuf) = get_connect_mat(MIM2, nroi, +1);
        CS_s(:, :, :, ishuf) = rm_components(permute(CS, [3 1 2 4]), npcs(1));
        cCOH_s(:, :, :, ishuf) = rm_components(permute(cCOH, [3 1 2 4]), npcs(1));
        aCOH_s(:, :, :, ishuf) = rm_components(permute(aCOH, [3 1 2 4]), npcs(1));
        iCOH_s(:, :, :, ishuf) = rm_components(permute(iCOH, [3 1 2 4]), npcs(1));
        wPLI_s(:, :, :, ishuf) = rm_components(wPLI, npcs(1));
    end

    % save methods in a struct
    clear out
    conn.inds = inds;
    for iout = 1:length(output)
        eval(['conn.' output{iout} ' = ' output{iout} '_s;'])
    end
    fprintf('\n');
end
