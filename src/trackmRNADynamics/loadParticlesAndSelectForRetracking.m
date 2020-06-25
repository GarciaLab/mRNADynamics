function Particles = loadParticlesAndSelectForRetracking(OutputFolder, NCh,retrack)
% Check if particle tracking has already been done on this dataset

%initialize Particles. replace if it already exists.
for ch = 1:NCh
    Particles{ch} = []; % This is the structure where we'll be tracking all particles.
end

if exist([OutputFolder, filesep, 'Particles.mat'], 'file')
    
    load([OutputFolder, filesep, 'Particles.mat'], 'Particles')
    
    % If there's only one channel, Particles, Spots and other structures are
    % not saved as cells. We turn them into a cell to ensure
    % compatibility.
    if NCh == 1
        Particles = {Particles};
    end
    
    for ch = 1:NCh
        
        if isfield(Particles{ch}, 'Approved') && retrack
            
            display(['Performing retracking on channel ', num2str(ch)])
            
            %Only keep the approved particles and start the tracking from there
            k = 1;
            
            for p = 1:length(Particles{ch})
                
                if Particles{ch}(p).Approved
                    NewParticles{ch}(k) = Particles{ch}(p);
                    k = k + 1;
                end
                
            end
            
            if exist('NewParticles', 'var')
                Particles{ch} = NewParticles{ch};
            else
                Particles{ch} = [];
            end
            
        else
            Particles{ch} = [];
        end
        
    end
    
end

end