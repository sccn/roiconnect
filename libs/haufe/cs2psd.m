function ps = cs2psd(cs)

% (C) 2018 Stefan Haufe
  
[nfreq nchan nchan njack] = size(cs);
  
coh = cs;
for ijack = 1:njack
  for ifreq = 1:nfreq
    ps(ifreq, :, ijack) = diag(squeeze(cs(ifreq, :, :, ijack)));
  end
end

            
