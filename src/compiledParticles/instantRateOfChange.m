function [CompiledParticles, MeanSlopeVectorAP, SDSlopeVectorAP, NSlopeAP]...
    = instantRateOfChange(CompiledParticles, ElapsedTime, ExperimentAxis, APFilter)
%INSTANTRATEOFCHANGE Summary of this function goes here
%   Detailed explanation goes here

FrameWindow=5; %AR 1/12/18 where did this number come from?

for ChN=1:NChannels
    if ~isempty(CompiledParticles{ChN})
        
        %Calculate the derivative as a function of time for each particle
        for j=1:length(CompiledParticles{ChN})
            [CompiledParticles{ChN}(j).SlopeTrace,...
                CompiledParticles{ChN}(j).SDSlopeTrace]=...
                TraceDerivative(CompiledParticles{ChN}(j),...
                ElapsedTime,FrameWindow);
        end
        
        
        if strcmpi(ExperimentAxis,'AP') || strcmpi(ExperimentAxis,'DV')
            %Calculate the average slope over an AP window
            MeanSlopeVectorAP{ChN}=nan(size(MeanVectorAP{ChN}));
            SDSlopeVectorAP{ChN}=nan(size(MeanVectorAP{ChN}));
            NSlopeAP{ChN}=nan(size(MeanVectorAP{ChN}));
            
            
            APBins=find(sum(APFilter{ChN}));
            
            
            
            for j=1:length(APBins)
                TraceCell=cell(length(ElapsedTime),1);
                
                ParticlesToAverage=find(APFilter{ChN}(:,APBins(j)));
                
                for k=1:length(ParticlesToAverage)
                    
                    for m=1:length(CompiledParticles{ChN}(ParticlesToAverage(k)).SlopeTrace)
                        
                        TraceCell{CompiledParticles{ChN}(ParticlesToAverage(k)).Frame(m)}=...
                            [TraceCell{CompiledParticles{ChN}(ParticlesToAverage(k)).Frame(m)},...
                            CompiledParticles{ChN}(ParticlesToAverage(k)).SlopeTrace(m)];
                    end
                end
                
                %Get rid of the nan in certain time points
                TraceCell=cellfun(@(x) x(~isnan(x)),TraceCell,'UniformOutput',false);
                
                MeanTrace=cellfun(@mean,TraceCell,'UniformOutput',false);
                SDTrace=cellfun(@std,TraceCell,'UniformOutput',false);
                NParticlesTrace=cellfun(@length,TraceCell,'UniformOutput',false);
                
                MeanSlopeVectorAP{ChN}(:,APBins(j))=[MeanTrace{:}];
                SDSlopeVectorAP{ChN}(:,APBins(j))=[SDTrace{:}];
                NSlopeAP{ChN}(:,APBins(j))=[NParticlesTrace{:}];
            end
        end
    end
end
end

