function [avgNuclearFluo, totalNumNuclei] = extractNuclearFluoDosage(prefix)

close all

%Load the folder information
[~,~,~,~,preProcPath]=...
    DetermineLocalFolders(prefix);
fileProperties = dir([preProcPath,filesep,prefix,filesep,'*_ch01.tif']);

totalNuclearFluo = 0;
totalNumNuclei = 0;
for i = 1:length(fileProperties)
    preProsImage=imread([preProcPath,filesep,prefix,filesep,fileProperties(i).name]);
%     preProsSharpenImage = imsharpen(preProsImage);
    [nuclearFluo,numNuclei] = integrateNuclearFluo(preProsImage);
    totalNuclearFluo = nuclearFluo + totalNuclearFluo;
    totalNumNuclei = numNuclei + totalNumNuclei;
end


avgNuclearFluo = totalNuclearFluo/totalNumNuclei;
end