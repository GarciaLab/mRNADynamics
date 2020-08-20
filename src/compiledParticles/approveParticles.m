function Particles = approveParticles(Particles,...
    shouldApproveAll, nSpotChannels, haveHistoneChannel)

%Approve all particles if the mode has been selected
if shouldApproveAll
    for ChN=1:nSpotChannels
        %Check that the approved field is present. If not include
        %it. This can occur if CheckParticleTracking is not run
        %first.
        if ~isfield(Particles{ChN},'Approved')
            for i=1:length(Particles{ChN})
                Particles{ChN}(i).Approved=0;
            end
        end
        
        %Make sure the particle has an associated nucleus if we are in
        %HistoneChannel mode
        if haveHistoneChannel
            for i=1:length(Particles{ChN})
                if ~isempty(Particles{ChN}(i).Nucleus)
                    %If a particle has been explicitly rejected then don't
                    %approve it!
                    if Particles{ChN}(i).Approved~=-1
                        Particles{ChN}(i).Approved=1;
                    end
                end
            end
        else % in case there's no histone channel
            for i=1:length(Particles{ChN})
                %If a particle has been explicitly rejected then don't
                %approve it!
                if Particles{ChN}(i).Approved~=-1
                    Particles{ChN}(i).Approved=1;
                end
            end
        end
    end
end

end