% Implements Bispectral Time Delay Estimations from Nikias Paper
% (1988). Inputs are precomputed univariate bispectra.
%
% Inputs:
%   B_xxx, B_xyx, B_yyy - (n_freq x n_freq) univariate bispectra
%   method              - determines the TDE method, must be between 1:4.
%   frqs_idx            - (n_freq x 1) vector of frequency bins?
%   zeropad             - zero-padding, used for frequency-domain bispectral estimates
%
% Outputs:
%   T - (1      x 2 * n_freq - 1) time-domain time delay estimate (bispectral hologram)
%   I - (n_freq x 2 * n_freq - 1) method-dependent coefficient to estimate bispectral delays

function [T,I] = bispectral_TD_est(B_xxx, B_xyx, B_yyy, method, frqs_idx, zeropad)

    if nargin < 6
        zeropad = 1;
    end

    assert(ismember(method,1:4), 'Method must be [1;4].');
    [m,n,boots] = size(B_xxx);
    
    switch method 
        case 1
            phi = angle(B_xyx)-angle(B_xxx);
            I = exp(1i.*phi);
        case 2
            phi = angle(B_xyx)- 0.5 * (angle(B_xxx) + angle(B_yyy));
            I = exp(1i.*phi);
        case 3
            I = B_xyx ./ B_xxx;
        case 4
            phi = angle(B_xyx)- 0.5 * (angle(B_xxx) + angle(B_yyy));
            I = abs(B_xyx).* exp(1i.*phi) ./ sqrt(abs(B_xxx).*abs(B_yyy));
    end
    
    % translate back to time domain.
    if nargin > 4 & ~isempty(frqs_idx)
        frqs_idx = [frqs_idx zeros(1,n-length(frqs_idx))];
        I = frqs_idx.' * frqs_idx .* I;
    end
    
    % use for frequency-domain bispec estimates (not time domain)
    if zeropad 
        I = [I zeros(m,n-1,boots)];
    end
    
    T = sum(I,1);
    T = abs(fftshift(ifft(T),2));
    
end