function compileNuclearTiffStacks(Prefix)

% get basic project info 
liveExperiment = LiveExperiment(Prefix);

% set input director
nucleusProbDir = [liveExperiment.procFolder filesep 'nucleusProbabilityMaps' filesep];

% get list of files in directory
nucleusProbFiles = dir([nucleusProbDir '*.tif']);

% iterate through list to obtain 
frameVecOrig = [];
for fileNum = 1:length(nucleusProbFiles)
    fName = nucleusProbFiles(fileNum).name;
    frameVecOrig(fileNum) = str2num(fName(length(fName)-11:length(fName)-9));
end

   
% load in probability maps and interpolate to obtain intermediate frame
% probabilities
wb = waitbar(0,'Interpolating nucleus mask files...');
for frame_i = 1:length(frameVecOrig)-1
  
    refFrames = frameVecOrig(frame_i):frameVecOrig(frame_i+1);
    
    if length(refFrames)>2
        % load stacks
        refStack1 = imread([nucleusProbFiles(frame_i).folder filesep nucleusProbFiles(frame_i).name]);
        refStack2 = imread([nucleusProbFiles(frame_i+1).folder filesep nucleusProbFiles(frame_i+1).name]);

        % linearize
        refStack1Lin = refStack1(:)';
        refStack2Lin = refStack2(:)';

        % interpolate    
        interpLin = interp1(refFrames([1 end]),vertcat(refStack1Lin,refStack2Lin),refFrames);

        % iterate through intermediate frames and write to file
        for frame_j = 2:length(refFrames)-1
            % shape back into image stacks
            tempStack = double(reshape(interpLin(frame_j,:),size(refStack1)));
            % write to file
            outName = ['prob' Prefix '_' sprintf('%03d',refFrames(frame_j)) '_ch00.tif'];
            imwrite(tempStack(:,:,1),[nucleusProbDir outName]);   
            for i = 2:size(refStack1,3)
                imwrite(tempStack(:,:,i),[nucleusProbDir outName],'WriteMode','append');   
            end
        end
    end
    waitbar(frame_i/(length(frameVecOrig)-1),wb);
end
disp('done.')
delete(wb);