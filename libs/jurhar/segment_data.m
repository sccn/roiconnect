function Xs = segment_data(X,seglen,shift)

%% Syntax
%     Xs = segment_data(X,seglen,shift)
%% Arguments
% _input_
%   X               (1,N) Time series of length N
%   seglen          integer, desired segment length
%   shift         integer, number of bins shifted for next segment,(optional)
% _output_
%     Xs            (nseg,seglen) output matrix of segmeneted timeseries
%% Description
% Segments input time series Xs (row vector) into segments of length seglen
% with desired shift.
%
% (C) Tin Jurhar, 2022
    if nargin < 3 | isempty(shift)
        shift = seglen;
    end

    assert(shift<=seglen,"shift must be geq seglen.")
    assert(length(size(X))==2 & (size(X,1)==1 | size(X,2)==1),"X must be 1-dimensional.")
    if size(X,2)==1
        X = X.';
    end
    assert(seglen <= length(X),'seglen must be geq than the input vector.')

    [~,lenX] = size(X);
    Xs = [];
    i = 1;
    while i+seglen-1 <= lenX
        Xs = [Xs; X(i:i+seglen-1)];
        i = i+shift;
    end
end