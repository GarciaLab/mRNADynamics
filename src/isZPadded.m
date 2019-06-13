function [zPadded, saveType, dogProb] = isZPadded(ProcPath, Prefix, spotChannel, saveType, dogs)

zPadded = false;

nameSuffix = ['_ch', iIndex(spotChannel, 2)];

dogProb = 'DOG_';

firstdogname = [dogProb, Prefix, '_', iIndex(1, 3), '_z', iIndex(1, 2), nameSuffix];

firstdogpath = [ProcPath, filesep,Prefix,'_',filesep,'dogs',filesep,firstdogname];

matsPresent = exist([firstdogpath, '.mat'], 'file');
tifsPresent = exist([firstdogpath, '.tif'], 'file');

if ~strcmpi(saveType, 'none') | isempty(saveType)
    
    if tifsPresent & ~matsPresent
        saveType = '.tif';
    elseif matsPresent & ~tifsPresent
        saveType = '.mat';
    elseif matsPresent & tifsPresent
        error('not sure which files to pick. check your processed folder.');
    end
    
end

firstdogpath = [firstdogpath, saveType];

if strcmpi(saveType, '.tif')
    try
        firstDoG = imread(firstdogpath);
    catch
        dogProb = 'prob';
        firstdogname = [dogProb, Prefix, '_', iIndex(1, 3), '_z', iIndex(1, 2), nameSuffix];
        firstdogpath = [ProcPath, filesep,Prefix,'_',filesep,'dogs',filesep,firstdogname,saveType];
        firstDoG = imread(firstdogpath);
    end
elseif strcmpi(saveType, '.mat')
    try
        load(firstdogpath, 'plane');
    catch
        dogProb = 'prob';
        firstdogname = [dogProb, Prefix, '_', iIndex(1, 3), '_z', iIndex(1, 2), nameSuffix];
        firstdogpath = [ProcPath, filesep,Prefix,'_',filesep,'dogs',filesep,firstdogname, saveType];
        load(firstdogpath, 'plane');
    end
    firstDoG = plane;
elseif strcmpi(saveType, 'none')
    firstDoG = dogs(:, :, 1, 1);
end

if sum(firstDoG(:)) == 0
    zPadded = true;
end

end