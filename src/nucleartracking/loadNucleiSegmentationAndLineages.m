function [Ellipses, schnitzcells] = loadNucleiSegmentationAndLineages(outputFolder, prefix)
% Check if particle tracking has already been done on this dataset

if exist([outputFolder, filesep, 'Ellipses.mat'], 'file') && ...
   exist([outputFolder, filesep, prefix, '_lin.mat'], 'file')
    
    load([outputFolder, filesep, 'Ellipses.mat'], 'Ellipses')
    load([outputFolder, filesep, prefix, '_lin.mat'], 'schnitzcells')
else
    error('Ellipses and/or _lin files missing from Results folder.')
end