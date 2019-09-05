function Particles = loadParticlesAndSelectForRetracking(OutputFolder, NCh,retrack)
% Check if particle tracking has already been done on this dataset

%initialize Particles. replace if it already exists.
for Channel = 1:NCh
    Particles{Channel} = []; % This is the structure where we'll be tracking all particles.
end

if exist([OutputFolder, filesep, 'Particles.mat'], 'file')
    
    load([OutputFolder, filesep, 'Particles.mat'], 'Particles')
    
    % If there's only one channel, Particles, Spots and other structures are
    % not saved as cells. We turn them into a cell to ensure
    % compatibility.
    if NCh == 1
        Particles = {Particles};
    end
    
    for Channel = 1:NCh
        
        if isfield(Particles{1}, 'Approved') && retrack
            display(['Performing retracking on channel ', num2str(Channel)])
            
            %Only keep the approved particles and start the tracking from there
            k = 1;
            
            for i = 1:length(Particles{Channel})
                
                if Particles{Channel}(i).Approved ~= 0
                    NewParticles{Channel}(k) = Particles{Channel}(i);
                    k = k + 1;
                end
                
            end
            
            if exist('NewParticles', 'var')
                Particles{Channel} = NewParticles{Channel};
            else
                Particles{Channel} = [];
            end
            
        else
            Particles{Channel} = [];
        end
        
    end
    
end

end