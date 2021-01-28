function [PlottedParams, PlottedParamSEs,ylab,OutputString,GlobalPlotYmax,GlobalPlotYmin,LogPlotYmin, includeLogPlots] = ...
    getEmbryoStatsPlottingVariables(this, parameter)
%% Load relevant parameters into memory
SchnitzCount = getEmbryoStatsParameters(this, 'SchnitzCount');
FractionSickNuclei = getEmbryoStatsParameters(this, 'FractionSickNuclei');
FractionRejectedNuclei = getEmbryoStatsParameters(this, 'FractionRejectedNuclei');
FractionCompleteNuclei = getEmbryoStatsParameters(this, 'FractionCompleteNuclei');
FractionFirstLastFrameNuclei = getEmbryoStatsParameters(this, 'FractionFirstLastFrameNuclei');
[TotalXDistanceTraveled, TotalXDistanceTraveledStd] = getEmbryoStatsParameters(this, 'MeanTotalXDistanceTraveled');
[TotalYDistanceTraveled, TotalYDistanceTraveledStd] = getEmbryoStatsParameters(this, 'MeanTotalYDistanceTraveled');
[TotalDistanceTraveled, TotalDistanceTraveledStd] = getEmbryoStatsParameters(this, 'MeanTotalDistanceTraveled');
[DistanceTraveledPerSecond, DistanceTraveledPerSecondStd] = getEmbryoStatsParameters(this, 'MeanDistanceTraveledPerSecond');
[TotalXDisplacement, TotalXDisplacementStd] = getEmbryoStatsParameters(this, 'MeanTotalXDisplacement');
[TotalYDisplacement, TotalYDisplacementStd] = getEmbryoStatsParameters(this, 'MeanTotalYDisplacement');
[TotalDisplacement, TotalDisplacementStd] = getEmbryoStatsParameters(this, 'MeanTotalDisplacement');
[DisplacementPerSecond, DisplacementPerSecondStd] = getEmbryoStatsParameters(this, 'MeanDisplacementPerSecond');
[NCDivisionTimes, NCDivisionTimesStd] = getEmbryoStatsParameters(this, 'NCDivisionInfo');

%%
if strcmpi(parameter, 'FractionSickNuclei')
    OutputString = 'FractionSickNuclei';
    PlottedParams = FractionSickNuclei;
    PlottedParamSEs = NaN(size(PlottedParams));
    ylab = 'sick nuclear fraction';
    GlobalPlotYmax = 1.05;
    GlobalPlotYmin = -0.05;
    LogPlotYmin = 0.05;
    includeLogPlots = false;
elseif strcmpi(parameter, 'FractionHealthyNuclei')
    OutputString = 'FractionHealthyNuclei';
    PlottedParams = 1-FractionSickNuclei;
    PlottedParamSEs = NaN(size(PlottedParams));
    ylab = 'healthy nuclear fraction';
    GlobalPlotYmax = 1.05;
    GlobalPlotYmin = -0.05;
    LogPlotYmin = 0.05;
    includeLogPlots = false;
elseif strcmpi(parameter, 'FractionRejectedNuclei')
    OutputString = 'FractionRejectedNuclei';
    PlottedParams = FractionRejectedNuclei;
    PlottedParamSEs = NaN(size(PlottedParams));
    ylab = 'rejected nuclear fraction';
    GlobalPlotYmax = 1.05;
    GlobalPlotYmin = -0.05;
    LogPlotYmin = 0.05;
    includeLogPlots = false;
elseif strcmpi(parameter, 'FractionApprovedNuclei')
    OutputString = 'FractionApprovedNuclei';
    PlottedParams = 1-FractionRejectedNuclei;
    PlottedParamSEs = NaN(size(PlottedParams));
    ylab = 'approved nuclear fraction';
    GlobalPlotYmax = 1.05;
    GlobalPlotYmin = -0.05;
    LogPlotYmin = 0.05;
    includeLogPlots = false;
elseif strcmpi(parameter, 'FractionCompleteNuclei')
    OutputString = 'FractionCompleteNuclei';
    PlottedParams = FractionCompleteNuclei;
    PlottedParamSEs = NaN(size(PlottedParams));
    ylab = 'complete cycle nuclear fraction';
    GlobalPlotYmax = 1.05;
    GlobalPlotYmin = -0.05;
    LogPlotYmin = 0.05;
    includeLogPlots = false;
elseif strcmpi(parameter, 'FractionFirstLastFrameNuclei')
    OutputString = 'FractionFirstLastFrameNuclei';
    PlottedParams = FractionFirstLastFrameNuclei;
    PlottedParamSEs = NaN(size(PlottedParams));
    ylab = 'first/last frame nuclear fraction';
    GlobalPlotYmax = 1.05;
    GlobalPlotYmin = -0.05;
    LogPlotYmin = 0.05;
    includeLogPlots = false;
elseif strcmpi(parameter, 'TotalXDistanceTraveled')
    OutputString = 'TotalXDistanceTraveled';
    PlottedParams = TotalXDistanceTraveled;
    PlottedParamSEs = TotalXDistanceTraveledStd;
    ylab = 'cumulative distance traveled along x axis per nucleus (\mum)';
    [GlobalPlotYmax, GlobalPlotYmin] = getEmbryoStatsGlobalYlims(PlottedParams, PlottedParamSEs);GlobalPlotYmax = min(500, GlobalPlotYmax*1.05);
    GlobalPlotYmin = max(-5, GlobalPlotYmin*0.95);
    LogPlotYmin = max(1, GlobalPlotYmin);
    includeLogPlots = false;
elseif strcmpi(parameter, 'TotalYDistanceTraveled')
    OutputString = 'TotalYDistanceTraveled';
    PlottedParams = TotalYDistanceTraveled;
    PlottedParamSEs = TotalYDistanceTraveledStd;
    ylab = 'cumulative distance traveled along y axis per nucleus (\mum)';
    [GlobalPlotYmax, GlobalPlotYmin] = getEmbryoStatsGlobalYlims(PlottedParams, PlottedParamSEs);GlobalPlotYmax = min(500, GlobalPlotYmax*1.05);
    GlobalPlotYmin = max(-5, GlobalPlotYmin*0.95);
    LogPlotYmin = max(1, GlobalPlotYmin);
    includeLogPlots = false;
elseif strcmpi(parameter, 'TotalDistanceTraveled')
    OutputString = 'TotalDistanceTraveled';
    PlottedParams = TotalDistanceTraveled;
    PlottedParamSEs = TotalDistanceTraveledStd;
    ylab = 'cumulative distance traveled per nucleus (\mum)';
    [GlobalPlotYmax, GlobalPlotYmin] = getEmbryoStatsGlobalYlims(PlottedParams, PlottedParamSEs);GlobalPlotYmax = min(500, GlobalPlotYmax*1.05);
    GlobalPlotYmin = max(-5, GlobalPlotYmin*0.95);
    LogPlotYmin = max(1, GlobalPlotYmin);
    includeLogPlots = false;    
elseif strcmpi(parameter, 'DistanceTraveledPerSecond')
    OutputString = 'DistanceTraveledPerSecond';
    PlottedParams = DistanceTraveledPerSecond;
    PlottedParamSEs = DistanceTraveledPerSecondStd;
    ylab = 'mean distance traveled per nucleus per second (\mum/s)';
    [GlobalPlotYmax, GlobalPlotYmin] = getEmbryoStatsGlobalYlims(PlottedParams, PlottedParamSEs);GlobalPlotYmax = min(5, GlobalPlotYmax*1.05);
    GlobalPlotYmin = max(-1, GlobalPlotYmin*0.95);
    LogPlotYmin = max(.01, GlobalPlotYmin);
    includeLogPlots = false;    
elseif strcmpi(parameter, 'TotalXDisplacement')
    OutputString = 'TotalXDisplacement';
    PlottedParams = TotalXDisplacement;
    PlottedParamSEs = TotalXDisplacementStd;
    ylab = 'mean x displacement per nucleus (\mum)';
   [GlobalPlotYmax, GlobalPlotYmin] = getEmbryoStatsGlobalYlims(PlottedParams, PlottedParamSEs);GlobalPlotYmax = min(100, GlobalPlotYmax*1.05);
    GlobalPlotYmin = max(-5, GlobalPlotYmin*0.95);
    LogPlotYmin = max(1, GlobalPlotYmin);
    includeLogPlots = false;         
elseif strcmpi(parameter, 'TotalYDisplacement')
    OutputString = 'TotalYDisplacement';
    PlottedParams = TotalYDisplacement;
    PlottedParamSEs = TotalYDisplacementStd;
    ylab = 'mean y displacement per nucleus (\mum)';
    [GlobalPlotYmax, GlobalPlotYmin] = getEmbryoStatsGlobalYlims(PlottedParams, PlottedParamSEs);GlobalPlotYmax = min(100, GlobalPlotYmax*1.05);
    GlobalPlotYmin = max(-5, GlobalPlotYmin*0.95);
    LogPlotYmin = max(1, GlobalPlotYmin);
    includeLogPlots = false;
elseif strcmpi(parameter, 'TotalDisplacement')
    OutputString = 'TotalDisplacement';
    PlottedParams = TotalDisplacement;
    PlottedParamSEs = TotalDisplacementStd;
    ylab = 'mean displacement per nucleus (\mum)';
    [GlobalPlotYmax, GlobalPlotYmin] = getEmbryoStatsGlobalYlims(PlottedParams, PlottedParamSEs);GlobalPlotYmax = min(100, GlobalPlotYmax*1.05);
    GlobalPlotYmin = max(-5, GlobalPlotYmin*0.95);
    LogPlotYmin = max(1, GlobalPlotYmin);
    includeLogPlots = false;
elseif strcmpi(parameter, 'DisplacementPerSecond')
    OutputString = 'DisplacementPerSecond';
    PlottedParams = DisplacementPerSecond;
    PlottedParamSEs = DisplacementPerSecondStd;
    ylab = 'mean displacement per nucleus per second (\mum/s)';
    [GlobalPlotYmax, GlobalPlotYmin] = getEmbryoStatsGlobalYlims(PlottedParams, PlottedParamSEs);
    GlobalPlotYmax = min(5, GlobalPlotYmax*1.05);
    GlobalPlotYmin = max(-1, GlobalPlotYmin*0.95);
    LogPlotYmin = max(.001, GlobalPlotYmin);
    includeLogPlots = false;
elseif strcmpi(parameter, 'NCDivisionTimes') | strcmpi(parameter, 'NCDivisionDurations') | strcmpi(parameter, 'NCDivisions') 
    OutputString = 'NCDivisions';
    PlottedParams = NCDivisionTimes;
    PlottedParamSEs = NCDivisionTimesStd;
    ylab = 'division cycle duration (min)';
    [GlobalPlotYmax, GlobalPlotYmin] = getEmbryoStatsGlobalYlims(PlottedParams, PlottedParamSEs);
    GlobalPlotYmax = min(100, GlobalPlotYmax*1.05);
    GlobalPlotYmin = max(0, GlobalPlotYmin*0.95);
    LogPlotYmin = max(1, GlobalPlotYmin);
    includeLogPlots = true; 
else
    error('Invalid choice of parameter.')
end

