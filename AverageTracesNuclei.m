function [MeanVector,SDVector,NParticles]=AverageTracesNuclei(FrameInfo,Nuclei)

%Average and SD over each time point. In order to do this we'll generate a
%cell array with all the values for a given time point

TraceCell=cell(length(FrameInfo),1);


for i=1:length(Nuclei)
    for j=1:length(Nuclei(i).Frames)
        TraceCell{Nuclei(i).Frames(j)}=[TraceCell{Nuclei(i).Frames(j)},...
            Nuclei(i).FluoMax(j)];
    end
end


MeanTrace=cellfun(@mean,TraceCell,'UniformOutput',false);
SDTrace=cellfun(@std,TraceCell,'UniformOutput',false);
NParticlesTrace=cellfun(@length,TraceCell,'UniformOutput',false);


MeanVector=[MeanTrace{:}];
SDVector=[SDTrace{:}];
NParticles=[NParticlesTrace{:}];
