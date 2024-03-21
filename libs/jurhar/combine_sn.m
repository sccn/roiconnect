function [X,Y] = combine_sn(rawX,noiseX,noiseY,delay,samplefreq,snr,theta,beta)
%   Combines signals and noise vectors with specified SNR (alpha), theta,
%   beta.

    [rawX,rawY] = mk_delay(delay,rawX);

    % zero-center signal and noise
    rawX = rawX-mean(rawX);
    rawY = rawY-mean(rawY);
    noiseX = noiseX-mean(noiseX);
    noiseY = noiseY-mean(noiseY);

    [~,~,~,pX] = obw(rawX,samplefreq);
    [~,~,~,pY] = obw(rawY,samplefreq);
    [~,~,~,pnoiseX] = obw(noiseX,samplefreq);
    [~,~,~,pnoiseY] = obw(noiseY,samplefreq);
    
    % Normalize signal and noise powers
    rawX = rawX / sqrt(pX);
    rawY = rawY /sqrt(pY);
    noiseX = noiseX / sqrt(pnoiseX);
    noiseY = noiseY /sqrt(pnoiseY);

    if nargin < 8 | isempty(beta)
        beta =1;
    end
    
    % Combine signal and noise with specified SNR
        X = snr *        rawX + (1-snr) * (noiseX + theta(1) * noiseY);
        Y = snr * beta * rawY + (1-snr) * (noiseY + theta(2) * noiseX);
end