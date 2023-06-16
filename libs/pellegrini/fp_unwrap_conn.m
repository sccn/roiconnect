function [MIM_, MIC_, GC_, DIFFGC_, iCOH_, aCOH_] = fp_unwrap_conn(conn,nroi,filt1,PCA_inds)

% Copyright (c) 2022 Franziska Pellegrini and Stefan Haufe

%%
%initialize output variables 
MIM_ = [];
MIC_ = [];
GC_=[];
DIFFGC_ = [];
iCOH_ =[];
aCOH_ = [];

% filt1
% filt1.band_inds

iinds = 0;
for iroi = 1:nroi
    for jroi = (iroi+1):nroi
        iinds = iinds + 1;
        
        if isfield(conn,'MIM')
            MIM_(iroi, jroi) = mean(conn.MIM(filt1.band_inds, iinds));
            MIM_(jroi,iroi) = MIM_(iroi,jroi);
        end
        if isfield(conn,'MIC')
            MIC_(iroi, jroi) = mean(conn.MIC(filt1.band_inds, iinds));
            MIC_(jroi,iroi) = MIC_(iroi,jroi);
        end
        if isfield(conn,'TRGC')
            DIFFGC_(iroi,jroi) = mean(squeeze(conn.TRGC(filt1.band_inds,iinds,1) - conn.TRGC(filt1.band_inds,iinds,2)));
            DIFFGC_(jroi,iroi) = -DIFFGC_(iroi,jroi);
        end
        if isfield(conn,'GC') 
           GC_(iroi,jroi) = mean(squeeze(conn.GC(filt1.band_inds,iinds,1) - conn.GC(filt1.band_inds,iinds,2)));
           GC_(jroi,iroi) = -GC_(iroi,jroi);
        end
    end
end

if isfield(conn,'COH')
    for iroi = 1:nroi
        for jroi = 1:nroi
            iCOH_(iroi, jroi) = mean(mean(mean(abs(imag(conn.COH(filt1.band_inds, PCA_inds{iroi}, PCA_inds{jroi}))),1), 2), 3);
            iCOH_(jroi,iroi) = iCOH_(iroi,jroi);
            aCOH_(iroi, jroi) = mean(mean(mean(abs(conn.COH(filt1.band_inds, PCA_inds{iroi}, PCA_inds{jroi})),1), 2), 3);
            aCOH_(jroi,iroi) = aCOH_(iroi,jroi);
        end
    end
end