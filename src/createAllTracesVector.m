function [AllTracesVector,AllTracesAP,AllTracesDV]= createAllTracesVector(FrameInfo,Particles,varargin)

%Order all traces in a time array. Also create a vector with the
%corresponding AP positions

AllTracesVector=zeros(length(FrameInfo),length(Particles));
AllTracesVector(:)=nan;

AllTracesAP=zeros(length(Particles),1);
AllTracesDV=zeros(length(Particles),1); %Added DV compatibility

for i=1:length(Particles)
    for j=1:length(Particles(i).Frame)
        AllTracesVector(Particles(i).Frame(j),i)=...
            Particles(i).Fluo(Particles(i).Frame==Particles(i).Frame(j));
    end
    
    if length(varargin)==1
        if strcmp(varargin{1},'NoAP')
            AllTracesAP(i)=nan;
        end
    else
        AllTracesAP(i)=Particles(i).MeanAP;
        AllTracesDV(i)=Particles(i).MeanDV; %Added DV compatibility
    end
end

