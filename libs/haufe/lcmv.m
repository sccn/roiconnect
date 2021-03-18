function [P, W, S, out] = lcmv(data, L, para)

% 2017 Stefan Haufe

data = data';
[M_ N NDUM]=size(L);

if nargin < 3
  para = [];
end

if ~isfield(para, 'alpha')
    para.alpha = nan;
end

if ~isfield(para, 'car')
    para.car = 0;    
end

if ~isfield(para, 'timeseries')
    para.timeseries = 1;
end

if ~isfield(para, 'onedim')
    para.onedim = 1;
end


if para.car == 1
    H = eye(M_)-ones(M_) ./ M_;
    [UH SH VH] = svd(H);
    out.H =  SH(1:end-1, :)*VH';
else
    out.H = eye(M_);
end

if ~isfield(para, 'noise')
    out.noisecov = eye(M_);
else
  if isequal(para.noise, para.noise')
      out.noisecov = para.noise';
  else
      if isnan(para.alpha)
        out.noisecov = shrinkage(para.noise', struct('Gamma', 'auto'));
      else
        out.noisecov = shrinkage(para.noise', struct('Gamma', para.alpha));
      end
  end
end

% out.noisecov=out.noisecov+para.alpha*eye(size(out.noisecov))*norm(out.noisecov);
out.spatfilt = real((out.H*out.noisecov*out.H')^(-0.5))*out.H;
% out.spatfilt = out.H;

if isequal(data, data')
    out.datacov = data;
    para.timeseries = 0;
else
    if isnan(para.alpha)
      [out.datacov out.alpha] = shrinkage(data, struct('Gamma', 'auto'));
    else
      out.datacov = shrinkage(data, struct('Gamma', para.alpha));
      out.alpha = para.alpha;
    end
    data = out.spatfilt*data;
end
    
% out.datacov=out.datacov+para.alpha*eye(size(out.datacov))*norm(out.datacov);

out.noisecov = out.spatfilt*out.noisecov*out.spatfilt';
out.datacov = out.spatfilt*out.datacov*out.spatfilt';
    
M = size(out.spatfilt, 1);
out.L = reshape(out.spatfilt*reshape(L, M_, N*NDUM), M, N, NDUM);

CI = pinv(out.datacov);
CInoise = pinv(out.noisecov);
if para.timeseries && nargout > 2
  CIdata = out.datacov\data;
end

P = zeros(N, 1);
if para.onedim
  W = zeros(M, N);
  if para.timeseries && nargout > 2
    S = zeros(size(data, 2), N);
  end
else
  W = zeros(M, N, 3);
  if para.timeseries && nargout > 2
    S = zeros(size(data, 2), N, 3);
  end
end

for ivox = 1:N
    Lvox= reshape(out.L(:, ivox, :), M, NDUM);
    
    E1 = inv(Lvox'*CI*Lvox);
    E2 = inv(Lvox'*CInoise*Lvox);
    
    W_ = (E1*Lvox'*CI*out.spatfilt)';
    out.W(:, ivox, :) = W_;

    if para.onedim
      [V, D] = eig(E1, E2);
      [ma in] = max(diag(D));
      V = V(:, in);
%       out.V1D(:, ivox) = V;
%       out.V3D(:, ivox, :) = W_;

      W(:, ivox) = W_*V./sqrt(V'*E2*V);   
      P(ivox) = D(in, in);
    else
      W(:, ivox, :) = W_./repmat(sqrt(diag(E2))', M, 1);
      P(ivox) = trace(E1)/trace(E2);
    end

    if para.timeseries && nargout > 2
      if para.onedim
        S(:, ivox) = (W(:, ivox)'*data)';
      else
        S(:, ivox, :) = (E1*Lvox'*CIdata)'./repmat(sqrt(diag(E2))', size(CIdata, 2), 1);    
      end
    else
        S = [];
    end
end
