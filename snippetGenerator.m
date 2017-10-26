function generateElongationSnippets(PrefixStrong, PrefixWeak, varargin)

    % extractSnippetsFor2ColorElongation
    %
    % DESCRIPTION
    % Isolate spots from different channels in 2 color experiments to compare
    % green and red time traces.
    %
    % ARGUMENTS
    % Prefixes for the elongation data (one prefix per channel)
    %
    % OPTIONS
    % None.
    %               
    % OUTPUT
    % Plots.
    %
    % Author (contact): Armando Reimer (areimer@berkeley.edu)
    % Created: 6/27/17
    % Last Updated: 8/31/2017
    %
    % Documented by: Armando Reimer (areimer@berkeley.edu)

    %Data for elongation experiment

    [SourcePath,FISHPath,DefaultDropboxFolder,MS2CodePath,PreProcPath]=...
        DetermineLocalFolders;


    pathWeak = [DefaultDropboxFolder, filesep,PrefixWeak, filesep];
    pathWeak = [DefaultDropboxFolder, filesep,PrefixStrong, filesep];
    weakspots = load([pathWeak, 'Spots.mat']);
    weakspots = Weakspots.Spots;
    weakframeinfo = load([pathWeak, 'FrameInfo.mat']);
    weakframeinfo = Weakframeinfo.FrameInfo;
    weakparticles = load([pathWeak, 'Particles.mat']);
    weakparticles = Weakparticles.Particles;
    strongspots = load([pathWeak, 'Spots.mat']);
    strongspots = Strongspots.Spots;
    strongparticles = load([pathWeak, 'Particles.mat']);
    strongparticles = Strongparticles.Particles;

    [SourcePathWeak,FISHPathWeak,DropboxFolderWeak,MS2CodePathWeak,PreProcPathWeak]=...
        DetermineLocalFolders(PrefixWeak);
    [SourcePathStrong,FISHPathStrong,DropboxFolderStrong,MS2CodePathStrong,PreProcPathStrong]=...
        DetermineLocalFolders(PrefixStrong);

    firstFrame = 1;

    for i=1:length(varargin)
        if strcmp(varargin{i},'firstFrame') && isnumeric(varargin{i+1})
            firstFrame = varargin{i+1};
        end
    end


    % loop over strong particles 

    for sp = 1:length(Strongparticles)
        % get particle frames
        strongFrames = Strongparticles(sp).Frames;

        %add 20 frames before the particle is detected, or as many as possible
        firstStrongFrame = strongFrames(1);
        if firstStrongFrame>20
            strongFrames = strongFrames-20:strongFrames(end);
        else
            strongFrames = 1:strongFrames(end);
        end
    end

end
    
    
    
    
    


