function cleanDir = getImageFiletypes(imageDir, isDogFolder)
%{
function for finding valid filetypes 
in a given directory for downstream
processing. if multiple valid filetypes are
found, the preference is:
1. stacks over planes
2. tifs over mats
3. probability maps over DoGs
%}


%isdogfolder is optional. default is false;
arguments
    imageDir struct
    isDogFolder logical = false;
end


isStack = contains(...
    string({imageDir.name}), '_z');

haveStacks = any(~isStack);
nStacks = sum(~isStack);

havePlanes = any(isStack);
nPlanes = sum(isStack);

isTif = contains(...
    string({imageDir.name}), '.tif');
haveTifs = any(isTif);
nTifs = sum(isTif);

isMat = contains(...
    string({imageDir.name}), '.mat');
haveMats = any(isMat);
nMats = sum(isMat);


if isDogFolder
    
    isProb =  startsWith(...
    string({imageDir.name}), 'prob');
    
    haveProbs = any(isProb);
    
    if haveProbs
        dogStr = 'prob';
    elseif haveStacks
        dogStr = 'dogStack_';
    elseif havePlanes
        dogStr = 'DOG_';
    else, error('No supported filetypes found.')
    end
    
end

%let's prefer stacks and tifs. 
%everything else is just to 
%ensure backwards compatibility
if haveStacks
    %we'll not include planes (if they exist)
    cleanDir = imageDir(isStack);
end

if haveTifs
    %we'll not include mats (if they exist)
    cleanDir = imageDir(isTif);
end

%if we have a dog folder, let's
%prefer probability maps
if haveProbs
    cleanDir = imageDir(isProb);
end


%By now, the cleaned directory should
%only list one set of file types. 

end