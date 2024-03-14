function [X,Y] = mk_delay(delay,sig1,sig2)
%% Syntax
%     [X,Y] = mk_delay(delay,sig1,sig2)
%% Arguments
% _input_
%   delay         integer, introduced delay
%   sig1          (1,N1) time series 1
%   sig2          (1,N2) time series 2 (optional)
% _output_
%     X           (1,N1-delay), sig1 derived time series
%     Y           (1,N2-delay), sig2 derived time series 
%% Description
% Takes time series sig1 and sig2 and time shifts them by delay bins. delay
% can take positive and negative values. If only sig1 is provided, then X,Y
% are time-shifted instances of sig1.
%
% (C) Tin Jurhar, 2022

    if nargin < 3
        sig2 = sig1;
    end

    if delay < 0
        X = sig1(1:end-abs(delay));
        Y = sig2(max(abs(delay))+1:end);
    else
        X = sig1(max(abs(delay))+1:end);
        Y = sig2(1:end-abs(delay));
    end
end

