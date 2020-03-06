function plotDorsalResultsLoop(varargin)

dataType = '1Dg-5_FFF';
activity = 'fraction';

%options must be specified as name, value pairs. unpredictable errors will
%occur, otherwise.
for i = 1:2:(numel(varargin)-1)
    if i ~= numel(varargin)
        eval([varargin{i} '=varargin{i+1};']);
    end
end


dorsalResults = plotFracByDlFluo2(dataType, activity); %dorsalResults is a struct for each nc 12, 13, 14
nc = 12;
dorsalResults = dorsalResults{nc-11};
x = dorsalResults.dorsalFluoBins;
y = dorsalResults.fracFluoEmbryo;
ymean = dorsalResults.meanFracFluoEmbryo;
se = dorsalResults.seFracFluoEmbryo;
n = 1;
R = 1;
plotDorsalActivity(x, y,activity, nc, dataType, ymean, se, 'fixRate', R, 'fixN', n)


end
