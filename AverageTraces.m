function [MeanVector,SDVector,NParticles,ErrorVarVector]=AverageTraces(FrameInfo,Particles)

%Average and SD over each time point. In order to do this we'll generate a
%cell array with all the values for a given time point

%2013/09/08: Modified to add 'ErrorSDVector', which uses bootstrapping to
%assign an error to the variance.

TraceCell=cell(length(FrameInfo),1);


for i=1:length(Particles)
    for j=1:length(Particles(i).Frame)
        TraceCell{Particles(i).Frame(j)}=[TraceCell{Particles(i).Frame(j)},...
            Particles(i).Fluo(j)];
    end
end


MeanTrace=cellfun(@mean,TraceCell,'UniformOutput',false);
SDTrace=cellfun(@std,TraceCell,'UniformOutput',false);
NParticlesTrace=cellfun(@length,TraceCell,'UniformOutput',false);

MeanVector=[MeanTrace{:}];
SDVector=[SDTrace{:}];
NParticles=[NParticlesTrace{:}];
ErrorVarVector=nan(size(MeanTrace));

%Estimate the error in the SD by bootstrapping
for i=1:length(TraceCell)
    if length(TraceCell{i})>1
        ErrorVarVector(i)=std(bootstrp(100,@var,TraceCell{i}));
    end
end



