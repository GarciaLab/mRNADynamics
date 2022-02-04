function membraneMat = getMembraneMat(liveExperiment)

membraneFile = [liveExperiment.preFolder, filesep,liveExperiment.Prefix, '-Membrane.tif'];
if exist(membraneFile, 'file')
    haveMemTifStack = true;
else
    haveMemTifStack = false;
end





if haveMemTifStack
    %load in sequential tif stacks
    membraneMat = imreadStack2([liveExperiment.preFolder, filesep,...
        liveExperiment.Prefix, '-Membrane.tif'], liveExperiment.yDim, liveExperiment.xDim, liveExperiment.nFrames);
    
    
else
    %load in individual tif slices
    error('Membrane Tif Stack does not exist. Run ExportDataForLivemRNA.')
    
end



%             end
membraneMat = double(membraneMat);