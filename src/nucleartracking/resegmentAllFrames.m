function Ellipses = resegmentAllFrames(Prefix, varargin)

hisMat = [];
%defaults. good for fly embryos (Dorsal synthetics settings)
% to change use them as options in TrackNuclei
sigmaK_um = 0.2;
min_rad_um = 2; 
max_rad_um = 6;
mu = 0.1;
nIterSnakes = 100;

%options must be specified as name, value pairs. unpredictable errors will
%occur, otherwise.
for k = 1:2:(numel(varargin)-1)
    if k ~= numel(varargin)
        eval([varargin{k} '=varargin{k+1};']);
    end
end

liveExperiment = LiveExperiment(Prefix);
pixelSize_um = liveExperiment.pixelSize_um;
hisMat = getHisMat(liveExperiment);


Ellipses = cell(liveExperiment.nFrames, 1);

parfor CurrentFrame = 1:liveExperiment.nFrames
        
        
        HisImage = hisMat(:, :, CurrentFrame);
        
        %segment current frame if nonempty
        %and add ellipse information to Ellipses
        if max(max(HisImage))~=0
            [~, ellipseFrame] = kSnakeCircles(HisImage,...
                pixelSize_um, 'min_rad_um', min_rad_um,...
                'max_rad_um', max_rad_um,'sigmaK_um',sigmaK_um,'mu', mu,...
            'nIterSnakes',nIterSnakes);    
            ellipseFrame(:, 6:9) = zeros(size(ellipseFrame, 1), 4);
            Ellipses{CurrentFrame} = ellipseFrame;
        else
            Ellipses{CurrentFrame} = [];
        end
                
        %report progress every tenth frame
        if ~mod(CurrentFrame, 10),...
                disp(['Segmenting frame: ' num2str(CurrentFrame), '...']); end
        
end

%TrackNuclei handles empty frames poorly, so let's fill them in. 
Ellipses = fillEmptyXYFrames(Ellipses);

save2([liveExperiment.resultsFolder, 'Ellipses.mat'], Ellipses);

ellipsesStats = getEllipsesStatistics(Ellipses);

% 
% figure(2); tiledlayout('flow');
% nexttile;
% hist(ellipsesStats.cenxList);
% nexttile;
% hist(ellipsesStats.cenyList);
% nexttile;
% hist(ellipsesStats.semiMajList);
% nexttile;
% hist(ellipsesStats.semiMinList);
% nexttile;
% hist(ellipsesStats.orientationAngleList);
% 
end