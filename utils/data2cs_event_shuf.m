% Calculates surrogate cross-spectra from data for event-related
% measurements using a permutation-based approach by shuffling the trials of the original data. 
%
% The off-diagonals of the cross-spectrum are calculated by multiplying the 
% unshuffled time-series with the trial-shuffled conjugate time-series. The diagonals (power) 
% are calculated by multiplying the *unshuffled* time-series with its conjugate self.
%
% Usage:
%   [cs, coh, wpli, nave] = data2cs_event(data_orig, data_shuf, segleng, segshift, epleng, maxfreqbin, para)
%
% Inputs:
%   data_orig  - (ndat x nchan) original (unshuffled) data matrix, each colum is the time-series in one channel
%   data_shuf  - (ndat x nchan) trial-shuffled data matrix
%   segleng    - length of each segment in bins, e.g. segleng = 1000
%   segshift   - numer of bins by which neighboring segments are shifted, 
%                e.g. segshift = segleng/2 makes overlapping segments
%   epleng     - length of each epoch
%   maxfreqbin - max. frequency in bins
%   para       - optional structure:
%                   para.segave = 0   -> no averaging across segments
%                   para.segave neq 0 -> averaging across segments (default is 0)
%                   para.subave = 1   -> subtracts the average across epochs,
%                   para.subave ~= 1  -> no subtraction (default is 1)
%                   IMPORTANT: if you only have one epoch (e.g. for continuous data), set para.subave = 0
%                   
%                   para.proj -> must be a set of vector in channel space, if it exists then the raw output 
%                   contains the single trial Fourier-transform in that channel
%
% Outputs:
%   cs   - (nchan x chan x maxfreqbin x nseg) tensor, cs(:,:,f,i) contains the cross-spectrum at frequency f and segment i
%   coh  - complex coherency calculated from cs
%   nave - number of averages
%
% Authors: 
%   Guido Nolte, g.nolte@uke.de
%   Tien Dung Nguyen, tien-dung.nguyen@charite.de

function [cs, coh, wpli, nave] = data2cs_event_shuf(data_orig, data_shuf, segleng, segshift, epleng, maxfreqbin, para)

    subave = 1;
    
    if nargin < 6
        para = [];
    end
    
    segave = 0;
    mydetrend = 0;
    proj = [];
    if isfield(para, 'segave')
        segave = para.segave;
    end
    if isfield(para, 'detrend')
        mydetrend = para.detrend;
    end
    if isfield(para, 'proj')
        proj = para.proj;
    end
    if isfield(para, 'subave')
        subave = para.subave;
    end
    if isfield(para, 'freqresolution')
        desired_nfreq = para.freqresolution;
    end
    
    if desired_nfreq == 0
        maxfreqbin = min([maxfreqbin,floor(segleng/2)+1]);
    else
        maxfreqbin = desired_nfreq + 1;
    end
    fres = segleng / 2;
    
    [ndum, npat] = size(proj);
    [ndat, nchan] = size(data_orig);
    
    if npat > 0
        data_orig = data_orig * proj;
        nchan = npat;
    end
    
    nep = floor(ndat/epleng);
    
    nseg = floor((epleng-segleng)/segshift)+1; % total number of segments
    
    if segave == 0
        cs = zeros(nchan, nchan, maxfreqbin, nseg);
        av = zeros(nchan, maxfreqbin, nseg);
    else
        cs = zeros(nchan, nchan, maxfreqbin);
        av = zeros(nchan, maxfreqbin);
    end
    
    if npat > 0
        if segave == 0
            cs = zeros(nchan, nchan, maxfreqbin, nep, nseg);
            av = zeros(nchan, maxfreqbin, nep, nseg);
        else
            cs = zeros(nchan, nchan, maxfreqbin, nep);
            av = zeros(nchan, maxfreqbin, nep);
        end
    end
    
    mywindow = repmat(hanning(segleng), 1, nchan);
    if isfield(para, 'mywindow')
        mywindow = repmat(para.mywindow, 1, nchan);
    end
    
    % figure; plot(mywindow);
    nave = 0;
    wpli_numer = cs; wpli_denom = cs; % initialize in the same way as CS
    
    for j = 1:nep
        dataep_orig = data_orig((j-1) * epleng+1:j * epleng, :);
        dataep_shuf = data_shuf((j-1) * epleng+1:j * epleng, :);
        for i = 1:nseg % average over all segments;
            datalocfft_orig = calc_fft(dataep_orig, i);
            datalocfft_shuf = calc_fft(dataep_shuf, i);
    
            for f = 1:maxfreqbin % for all frequencies
                if npat == 0
                    if segave == 0
                        cs(:,:,f,i) = cs(:,:,f,i) + calc_cs(datalocfft_orig(f,:), datalocfft_shuf(f,:));
                        av(:,f,i) = av(:,f,i) + conj(datalocfft_orig(f,:)');
                        wpli_numer(:,:,f,i) = wpli_numer(:,:,f,i) + imag(calc_cs(datalocfft_orig(f,:), datalocfft_shuf(f,:)));
                        wpli_denom(:,:,f,i) = wpli_denom(:,:,f,i) + abs(imag(calc_cs(datalocfft_orig(f,:), datalocfft_shuf(f,:))));
                    else
                        %disp([i,f,size(datalocfft)])
                        cs(:,:,f) = cs(:,:,f) + calc_cs(datalocfft_orig(f,:), datalocfft_shuf(f,:));
                        av(:,f) = av(:,f) + conj(datalocfft_orig(f,:)');
                        wpli_numer(:,:,f) = wpli_numer(:,:,f) + imag(calc_cs(datalocfft_orig(f,:), datalocfft_shuf(f,:)));
                        wpli_denom(:,:,f) = wpli_denom(:,:,f) + abs(imag(calc_cs(datalocfft_orig(f,:), datalocfft_shuf(f,:))));
                    end
                else
                    if segave == 0
                        cs(:,:,f,j,i) = calc_cs(datalocfft_orig(f,:), datalocfft_shuf(f,:));
                        av(:,f,j,i) = conj(datalocfft_orig(f,:)');
                        wpli_numer(:,:,f,j,i) = imag(calc_cs(datalocfft_orig(f,:), datalocfft_shuf(f,:)));
                        wpli_denom(:,:,f,j,i) = abs(imag(calc_cs(datalocfft_orig(f,:), datalocfft_shuf(f,:))));
                    else
                        %disp([i,f,size(datalocfft)])
                        cs(:,:,f,j) = cs(:,:,f,j) + calc_cs(datalocfft_orig(f,:), datalocfft_shuf(f,:));
                        av(:,f,j) = av(:,f,j) + conj(datalocfft_orig(f,:)');
                        wpli_numer(:,:,f,j) = wpli_numer(:,:,f,j) + imag(calc_cs(datalocfft_orig(f,:), datalocfft_shuf(f,:)));
                        wpli_denom(:,:,f,j) = wpli_denom(:,:,f,j) + abs(imag(calc_cs(datalocfft_orig(f,:), datalocfft_shuf(f,:))));
                    end
                end
            end
        end
        nave = nave+1;
    end
    
    if segave == 0
        cs = cs / nave;
        av = av / nave;
        wpli_numer = wpli_numer / nave;
        wpli_denom = wpli_denom / nave;
    else
        nave = nave * nseg;
        cs = cs / nave;
        av = av / nave;
        wpli_numer = wpli_numer / nave;
        wpli_denom = wpli_denom / nave;
    end
    
    % compute wPLI
    wpli = abs(wpli_numer) ./ wpli_denom;
    wpli = permute(wpli, [3, 1, 2]);
    
    for f = 1:maxfreqbin
        if subave == 1
            if npat == 0
                if segave == 0
                    for i = 1:nseg
                        cs(:,:,f,i) = cs(:,:,f,i) - av(:,f,i) * av(:,f,i)';
                    end
                else
                    cs(:,:,f) = cs(:,:,f) - av(:,f) * av(:,f)';
                end
            else
                if segave == 0
                    for i = 1:nseg
                        for j = 1:nep
                            cs(:,:,f,j,i) = cs(:,:,f,j,i) - av(:,f,j,i) * av(:,f,j,i)';
                        end
                    end
                else
                    for j = 1:nep
                        cs(:,:,f,j) = cs(:,:,f,j) - av(:,f,j) * av(:,f,j)';
                    end
                end
            end
        end
    end
    
    ndim = length(size(cs));
    if ndim == 3
        [n1, n2, n3] = size(cs);
        coh = cs;
        for i = 1:n3
            c = squeeze(cs(:, :, i));
            coh(:,:,i) = c ./ sqrt(diag(c) * diag(c)');
        end
    elseif ndim == 4
        [n1, n2, n3, n4] = size(cs);
        coh = cs;
        for i = 1:n3
            for j = 1:n4
                c = squeeze(cs(:, :, i, j));
                coh(:,:,i,j) = c ./ sqrt(diag(c) * diag(c)');
            end
        end
    end

    function datalocfft = calc_fft(dataep, i) % nested function that has access to all variables
        dataloc = dataep((i-1) * segshift + 1:(i-1) * segshift + segleng, :);
    
        % apply Hanning window
        if mydetrend == 1
            dataloc_win = detrend(dataloc,0) .* mywindow;
        else
            dataloc_win = dataloc .* mywindow;
        end
    
        % zero padding if necessary
        if desired_nfreq ~= 0
            required_zeros = desired_nfreq - fres;
            if required_zeros < 0
                error('Desired frequency resolution cannot be lower than the actual resolution of the signal.')
            end
            pad = zeros(required_zeros, size(dataloc_win, 2));
            dataloc_win = cat(1, pad, dataloc_win, pad);
        end
        datalocfft = fft(dataloc_win);
    end

    function cs_ = calc_cs(datalocfft_orig, datalocfft_shuf)
        cs_ = conj(datalocfft_orig' * datalocfft_shuf);
        power = conj(datalocfft_orig' * datalocfft_orig); 
        
        % use power from the original (unshuffled) cross-spectrum 
        cs_ = cs_ - diag(diag(cs_)) + diag(diag(power));
    end
end


