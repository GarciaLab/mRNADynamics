function compileNuclearTiffStacks(Prefix)

% get basic project info 
liveExperiment = LiveExperiment(Prefix);

% set input director
nucleusProbDir = [liveExperiment.procFolder 'nucleusProbabilityMaps' filesep];
nucleusProbDirFinal = [liveExperiment.procFolder 'nucleusProbabilityMapsFull' filesep];
mkdir(nucleusProbDirFinal);

% load info file
load([nucleusProbDir 'nucleusInfo.mat'], 'nucleusInfo')

% get list of files to load
nucleusProbFiles = nucleusInfo.originalFileNames;
frameVecOrig = nucleusInfo.originalFrames;

   
% load in probability maps and interpolate to obtain intermediate frame
% probabilities
% wb = waitbar(0,'Interpolating nucleus mask files...');
disp('Interpolating nucleus mask files...');
parfor frame_i = 1:length(frameVecOrig)-1
  
    refFrames = frameVecOrig(frame_i):frameVecOrig(frame_i+1);
    
    if true%length(refFrames)>2
        % load stacks
        refStack1 = imreadStack([nucleusProbDir  'prob' nucleusProbFiles{frame_i}]);
        refStack2 = imreadStack([nucleusProbDir  'prob' nucleusProbFiles{frame_i+1}]);        
        % linearize
        refStack1Lin = refStack1(:)';
        refStack2Lin = refStack2(:)';

        % interpolate    
        interpLin = interp1(refFrames([1 end]),vertcat(refStack1Lin,refStack2Lin),refFrames);

        % iterate through intermediate frames and write to file
        for frame_j = 1:length(refFrames)
            % shape back into image stacks
            tempStack = cast(reshape(interpLin(frame_j,:),size(refStack1)),'uint16');
            % write to file
            outName = ['prob' Prefix '_' sprintf('%03d',refFrames(frame_j)) '_ch00.tif'];
            imwrite(tempStack(:,:,1),[nucleusProbDirFinal outName],'WriteMode','append');   
            for i = 2:size(refStack1,3)
                imwrite(tempStack(:,:,i),[nucleusProbDirFinal outName],'WriteMode','append');   
            end
        end
    end
%     waitbar(frame_i/(length(frameVecOrig)-1),wb);
end
disp('done.')
% delete(wb);