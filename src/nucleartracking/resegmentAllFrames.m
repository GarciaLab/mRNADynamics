function resegmentAllFrames(Prefix, varargin)

hisMat = [];
min_rad_um = 2; %default. good for fly embryos. 
max_rad_um = 6;

%options must be specified as name, value pairs. unpredictable errors will
%occur, otherwise.
for k = 1:2:(numel(varargin)-1)
    if k ~= numel(varargin)
        eval([varargin{k} '=varargin{k+1};']);
    end
end


thisExperiment = liveExperiment(Prefix);

if isempty(hisMat)    
    hisMat = getHisMat(thisExperiment);
end
% if exist([thisExperiment.resultsFolder, filesep, 'Ellipses.mat'], 'file')
%     Ellipses = getEllipses(thisExperiment);
% else
Ellipses = cell(thisExperiment.nFrames, 1);
% end

parfor CurrentFrame = 1:thisExperiment.nFrames
        
        HisImage = hisMat(:, :, CurrentFrame);
        [~, circles] = kSnakeCircles(HisImage,...
            thisExperiment.pixelSize_nm/1000, 'min_rad_um', min_rad_um,...
            'max_rad_um', max_rad_um);    
%         circles(:, 4) = circles(:, 3);
        circles(:, 6:9) = zeros(size(circles, 1), 4);
        Ellipses{CurrentFrame} = circles;
                
        %report progress every tenth frame
        if ~mod(CurrentFrame, 10), disp(['Segmenting frame: ' num2str(CurrentFrame), '...']); end
        
end

   %to prevent errors later on, let's fill in empty frames. 
    if sum(cellfun(@(x) size(x,1),Ellipses)==0)
        %Find the frames where we have issues
        FramesToFix = find(cellfun(@(x) size(x,1),Ellipses)==0);
        for k=1:length(FramesToFix)
            if FramesToFix(k)==1
                FrameToCopy=1;
                while sum(FramesToFix==NextFrameToCopy)
                    FrameToCopy=FrameToCopy+1;
                end
            else
                FrameToCopy=FramesToFix(k)-1;
            end
            Ellipses{FramesToFix(k)}=Ellipses{FrameToCopy};
        end
    end


save([thisExperiment.resultsFolder, 'Ellipses.mat'],'Ellipses', '-v6');

TrackNuclei(Prefix,'retrack');


end