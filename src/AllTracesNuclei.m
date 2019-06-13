function [AllTracesVector,AllTracesAP,AllTracesDV]=AllTracesNuclei(FrameInfo,Nuclei,varargin)

%Order all traces in a time array. Also create a vector with the
%corresponding AP positions

AllTracesVector=nan(length(FrameInfo),length(Nuclei));
AllTracesAP=nan(length(Nuclei),1);
AllTracesDV=nan(length(Nuclei),1);

for nuc=1:length(Nuclei)
    for frame=1:length(Nuclei(nuc).Frames)
        
        AllTracesVector(Nuclei(nuc).Frames(frame),nuc)=...
        Nuclei(nuc).FluoMax(Nuclei(nuc).Frames==Nuclei(nuc).Frames(frame));
    
        if length(varargin)==1
            if strcmp(varargin{1},'NoAP')
                %do nothing
            end
        else
            AllTracesAP(nuc)= Nuclei(nuc).MeanAP;
            AllTracesDV(nuc)= Nuclei(nuc).MeanDV;
        end
    end
   
end

