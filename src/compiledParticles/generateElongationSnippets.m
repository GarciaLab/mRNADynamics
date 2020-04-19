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
    pathStrong = [DefaultDropboxFolder, filesep,PrefixStrong, filesep];
    strongFrameInfo = load([pathStrong, 'FrameInfo.mat']);
    strongFrameInfo = strongFrameInfo.FrameInfo;
    strongParticles = load([pathStrong, 'CompiledParticles.mat'], 'CompiledParticles'); 
    strongParticles = strongParticles.CompiledParticles;
    zSize= strongFrameInfo(1).NumberSlices + 2;

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
    AllSnippets = struct;
    for sp = 1:length(strongParticles)
        % get particle frames
        strongFrames = strongParticles(sp).Frame;
        strongX = strongParticles(sp).xPos;
        strongY = strongParticles(sp).yPos;
        %add 20 frames before the particle is detected, or as many as possible
        firstStrongFrame = strongFrames(1);
        firstStrongX = strongX(1);
        firstStrongY = strongY(1);
        numBackTrackFrames = 20;
        snippetSize = 16; %pixels width

        if firstStrongFrame>numBackTrackFrames
            strongFrames = (firstStrongFrame-numBackTrackFrames):strongFrames(end);
            strongX =[ones(numBackTrackFrames,1)'*firstStrongX,strongX];
            strongY =[ones(numBackTrackFrames,1)'*firstStrongY,strongY];
        else
            strongFrames = firstFrame:strongFrames(end);
            strongX = [ones(firstStrongFrame-1,1)'*firstStrongX,strongX];
            strongY = [ones(firstStrongFrame-1,1)'*firstStrongY,strongY];
        end
        %
        
        strongMaskMovie = [];
        weakMaskMovie = [];
        
        for sframe = 1:length(strongFrames) %loop over z to max project
            %initialize max projected frame
            projTempStrong = zeros(strongFrameInfo(1).LinesPerFrame,strongFrameInfo(1).PixelsPerLine);
            projTempWeak = projTempStrong;
            
            
            %project frame
            for zslice = 1:zSize
                strongSlice = double(imread([PreProcPathStrong,filesep,PrefixStrong,filesep,PrefixStrong,'_',iIndex(sframe,3),'_z',iIndex(zslice,2),'.tif']));      
                projTempStrong = max(strongSlice,projTempStrong);                
                weakSlice = double(imread([PreProcPathWeak,filesep,PrefixWeak,filesep,PrefixWeak,'_',iIndex(sframe,3),'_z',iIndex(zslice,2),'.tif']));      
                projTempWeak= max(weakSlice,projTempWeak);
                
            end
            
            %generate the mask. 
            
            maskStrong = zeros(size(projTempStrong)); 
            maskWeak = maskStrong; %if these aren't the same dimensions, this algorithm obviously won't work.              
            
            % check if the particle is too close to the edges
            % first get the position of the snippet in this frame
            yo = strongY(sframe) - snippetSize/2;
            y1 = strongY(sframe) + snippetSize/2;
            xo = strongX(sframe) - snippetSize/2;
            x1 = strongX(sframe) + snippetSize/2;
            
            if yo > 0 && xo > 0 && y1 < size(maskStrong,1) && ...
               x1 < size(maskStrong,2)
                maskStrong(yo:y1, xo:x1) = projTempStrong(yo:y1, xo:x1);
                snipStrong = projTempStrong(yo:y1,xo:x1);
                bgStrong = mode(snipStrong(:));
%                 maskStrong(yo:y1, xo:x1) = maskStrong(yo:y1, xo:x1) - bgStrong;
                maskWeak(yo:y1, xo:x1) = projTempWeak(yo:y1, xo:x1);
                snipWeak = projTempWeak(yo:y1,xo:x1);
                bgWeak = mode(snipWeak(:));
%                 maskWeak(yo:y1, xo:x1) = maskWeak(yo:y1, xo:x1) - bgWeak;

                % make a 3D array (rows,columns,frames, i.e a movie) of
                % snippets of this nucleus                
            else %if it was too close make a mask of NaNs
                maskStrong = nan(size(maskStrong));
                maskWeak = maskStrong;
            end
            strongMaskMovie = cat(3, strongMaskMovie, maskStrong);
            weakMaskMovie = cat(3, weakMaskMovie, maskWeak);
        end
        % save the snippet movie of this particle into a format containing
        % all snippet movies of all particles
        AllSnippets(sp).strongSnippetsMovie = strongMaskMovie;
        AllSnippets(sp).weakSnippetsMovie = weakMaskMovie;
    end
    save([DefaultDropboxFolder,filesep,PrefixStrong,filesep,'AllSnippets.mat'],'AllSnippets')
    'debug here';
    
%     %% Look at the results
%     for nucleus = 1:length(AllSnippets)
%         frameNumber = size(AllSnippets(nucleus).strongSnippetsMovie,3);
%         strongMovie = AllSnippets(nucleus).strongSnippetsMovie;
%         weakMovie = AllSnippets(nucleus).weakSnippetsMovie;
%         for frame = 1:frameNumber
%             strongImage = strongMovie(:,:,frame);
%             weakImage = weakMovie(:,:,frame);
%             if sum(~isnan(strongImage))
%                 figure
%                 subplot(2,1,1)
%                 imshow(strongImage,[]);
%                 title('strong');
%                 subplot(2,1,2)
%                 imshow(weakImage,[]);
%                 title('weak');
%                 waitforbuttonpress
%                 close all force
%             end
%         end
end 