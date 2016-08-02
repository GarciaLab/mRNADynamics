function [AllTracesVector,AllTracesAP]=AllTracesNuclei(FrameInfo,Nuclei)

%Order all traces in a time array. Also create a vector with the
%corresponding AP positions

AllTracesVector=zeros(length(FrameInfo),length(Nuclei));
AllTracesVector(:)=nan;

AllTracesAP=zeros(length(Nuclei),1);

for i=1:length(Nuclei)
    for j=1:length(Nuclei(i).Frames)
        AllTracesVector(Nuclei(i).Frames(j),i)=...
     Nuclei(i).FluoMax(Nuclei(i).Frames==Nuclei(i).Frames(j));
     AllTracesAP(i)= Nuclei(i).MeanAP;
    end
   
end

