% Function to generate simulated signals. Signals are generated either as
% as samples from an exponential distribution, a Gaussian distribution or
% from a Gaussian pink noise process. The option to model the signal as
% spikes is disabled here.
%
% Inputs:
%   t           - [integer] length of signal in datapoints 
%   source_type - [String] underlying model of the signal: 'exp', 'G', 'pink', 'pinksq'
%   param       - [integer] noise parameter
%   noise_mask  - [boolean] if 1, additional Gaussian noise is added to the signal
%
% Output:
%   signal - (1 x t) simulated signal of length t

function [signal] = mk_series(t,source_type,param,noise_mask)
     switch source_type
        case 'exp'
            signal= exprnd(param,1,t);
        case 'G' % Gaussian Noise
            signal = randn(1,t);
        case 'pink'
            signal = mkpinknoise(t,1,param).';
        case 'pinksq'
            signal = mkpinknoise(t,1,param).';
            signal = signal.^2;
%         case 'spiky'
%             load('spike_signal.mat');
%             signal = spikes_signal(randperm(t)).';
     end

     if nargin >3 & noise_mask % add Gaussian noise for numerical stability
        signal = signal./norm(signal,'fro');
        noise= randn(1,t);
        noise = noise./ norm(noise,'fro');
        
        signal = 0.95 * signal + 0.05 * noise;
     end
end
