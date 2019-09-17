function coh = cs2coh(cs)

% (C) 2018 Stefan Haufe
  
[nfreq nchan nchan njack] = size(cs);
  
coh = cs;
for ijack = 1:njack
  for ifreq = 1:nfreq
    coh(ifreq, :, :, ijack) = squeeze(cs(ifreq, :, :, ijack)) ...
      ./ sqrt(diag(squeeze(cs(ifreq, :, :, ijack)))*diag(squeeze(cs(ifreq, :, :, ijack)))');
  end
end

            
