scriptPath = fileparts(mfilename('fullpath'));
jar = fullfile(scriptPath, 'tr2d.jar');
projectPath = 'C:\Users\jug\LivemRNA\Data\PreProcessedData\2015-05-31-89B8-3-P2P\tr2dProject';
exportPath = 'C:\Users\jug\LivemRNA\Data\PreProcessedData\2015-05-31-89B8-3-P2P\tr2dProject\mRNADynamicsExport';

%In case you want to put a certain classifier in place before Tr2d starts
%Copy classifier file to $projectPath/segmentation/weka
%Create a savedState.csv file in there, containing a line per model file
% e.g.: <filename.model>, 0.5,

system(['java -jar ', jar, ' -p ', projectPath, ' -e ', exportPath]);