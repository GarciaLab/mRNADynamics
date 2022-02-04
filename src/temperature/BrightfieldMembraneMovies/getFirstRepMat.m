

function outmat = getFirstRepMat(liveExperiment)



load([liveExperiment.resultsFolder,filesep,'MarkAndFindInfo.mat'], 'MarkAndFindInfo');



preTifDir = dir([liveExperiment.preFolder, '*_ch0*.tif']);

%just return an empty array if we can't load the movie.
%leave the handling to the caller, presumably by enabling
%sequential file loading.
if ~haveSufficientMemory(preTifDir)
    outmat = [];
    return;
end


exportedChannels = [];
% find what channels were exported
for k = 1:6  %i don't see channel number going beyond 6 any time soon.
    exportedChannels(k) =  any(contains(...
        string({preTifDir.name}), ['_ch0',num2str(k)]));
end
channelsToRead = find(exportedChannels);

% this is for backwards compatibility,
%exported tiffs used to be one per z slice.
haveTifStacks = any(~contains(...
    string({preTifDir.name}), '_z'));



if haveTifStacks
    
    
    moviePrecision = 'uint16';
    movieMat = zeros(liveExperiment.yDim, liveExperiment.xDim,...
        liveExperiment.zDim, MarkAndFindInfo.NSeries,...
        length(channelsToRead), moviePrecision); % y x z t ch
    
    chIndex = 0;
    
    for ch = channelsToRead
        
        chIndex = chIndex + 1;
        
        preChDir = preTifDir( ...
            contains(...
            string({preTifDir.name}), ['001_ch0', num2str(ch)]) &...
            ~contains(string({preTifDir.name}), '_z') );
        
        %making these temporary variables to avoid passing all
        %of
        %liveExperiment to the parforloop
        this_nEmbryos = MarkAndFindInfo.NSeries;
        this_yDim = liveExperiment.yDim;
        this_preFolder = liveExperiment.preFolder;
        this_xDim = liveExperiment.xDim;
        this_zDim = liveExperiment.zDim;
        
        for f = 1:this_nEmbryos
            movieMat(:, :, :, f, chIndex) =...
                imreadStack2([this_preFolder, filesep, preChDir(f).name],...
                this_yDim, this_xDim, this_zDim);
        end
        
        
    end
    
else
    error('can''t load movie.')
end


outmat = movieMat;

%let's reduce the memory footprint of the movie if we can
%             if max(movieMat(:)) < 255
%                 movieMat = uint8(movieMat);
%             end
outmat = double(outmat);

end