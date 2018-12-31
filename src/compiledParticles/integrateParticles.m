function CompiledParticles = integrateParticles(NChannels, ...
    ElapsedTime, CompiledParticles)
%INTEGRATEPARTICLES Summary of this function goes here
%   Detailed explanation goes here
for ChN=1:NChannels
    
    %In order to take this seriously I need to come up with a way to deal with
    %the particles that survive during mitosis.
    
    for i=1:length(CompiledParticles{ChN})
        if length(ElapsedTime(CompiledParticles{ChN}(i).Frame))>1
            CompiledParticles{ChN}(i).TotalmRNA=trapz(ElapsedTime(CompiledParticles{ChN}(i).Frame),CompiledParticles{ChN}(i).Fluo);
            
            %Estimate the error
            if length(CompiledParticles{ChN}(i).Frame)==2
                CompiledParticles{ChN}(i).TotalmRNAError=(ElapsedTime(CompiledParticles{ChN}(i).Frame(2))-...
                    ElapsedTime(CompiledParticles{ChN}(i).Frame(1)))/2*...
                    CompiledParticles{ChN}(i).FluoError*sqrt(2);
            else
                
                ErrorTemp=[];
                %Calculate the error of the inner points
                for j=2:(length(CompiledParticles{ChN}(i).Frame)-1)
                    ErrorTemp(j)=(ElapsedTime(CompiledParticles{ChN}(i).Frame(j+1))-...
                        ElapsedTime(CompiledParticles{ChN}(i).Frame(j-1)))*...
                        CompiledParticles{ChN}(i).FluoError;
                end
                
                %Calculate the error of the outer points
                ErrorTemp(1)=(ElapsedTime(CompiledParticles{ChN}(i).Frame(2))-...
                    ElapsedTime(CompiledParticles{ChN}(i).Frame(1)))/2*...
                    CompiledParticles{ChN}(i).FluoError;
                
                ErrorTemp(length(CompiledParticles{ChN}(i).Frame))=...
                    (ElapsedTime(CompiledParticles{ChN}(i).Frame(end))-...
                    ElapsedTime(CompiledParticles{ChN}(i).Frame(end-1)))/2*...
                    CompiledParticles{ChN}(i).FluoError;
                
                %Now, add it all up
                CompiledParticles{ChN}(i).TotalmRNAError=sqrt(sum(ErrorTemp.^2));
                
                
            end
        else
            CompiledParticles{ChN}(i).TotalmRNA=[];
        end
    end
end
end

