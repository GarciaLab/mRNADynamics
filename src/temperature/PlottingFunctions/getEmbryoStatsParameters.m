function [params, paramSEs] = getEmbryoStatsParameters(this, ParamType)

% First make sure ParamType is a valid argument
if strcmpi(ParamType, 'SchnitzCount') 
    paramName = 'SchnitzCount';
    seName = '';
elseif strcmpi(ParamType, 'FractionSickNuclei') 
    paramName = 'FractionSickNuclei';
    seName = '';
elseif strcmpi(ParamType, 'FractionRejectedNuclei') 
    paramName = 'FractionRejectedNuclei';
    seName = '';
elseif strcmpi(ParamType, 'FractionCompleteNuclei') 
    paramName = 'FractionCompleteNuclei';
    seName = '';
elseif strcmpi(ParamType, 'FractionFirstLastFrameNuclei') 
    paramName = 'FractionFirstLastFrameNuclei';
    seName = '';
elseif strcmpi(ParamType, 'MeanTotalXDistanceTraveled') 
    paramName = 'MeanTotalXDistanceTraveled';
    seName = 'StdTotalXDistanceTraveled'; 
elseif strcmpi(ParamType, 'MeanTotalYDistanceTraveled') 
    paramName = 'MeanTotalYDistanceTraveled';
    seName = 'StdTotalYDistanceTraveled';  
elseif strcmpi(ParamType, 'MeanTotalDistanceTraveled') 
    paramName = 'MeanTotalDistanceTraveled';
    seName = 'StdTotalDistanceTraveled';  
elseif strcmpi(ParamType, 'MeanDistanceTraveledPerSecond') 
    paramName = 'MeanDistanceTraveledPerSecond';
    seName = 'StdDistanceTraveledPerSecond';
elseif strcmpi(ParamType, 'MeanTotalXDisplacement') 
    paramName = 'MeanTotalXDisplacement';
    seName = 'StdTotalXDisplacement';
elseif strcmpi(ParamType, 'MeanTotalYDisplacement') 
    paramName = 'MeanTotalYDisplacement';
    seName = 'StdTotalYDisplacement';
elseif strcmpi(ParamType, 'MeanTotalDisplacement') 
    paramName = 'MeanTotalDisplacement';
    seName = 'StdTotalDisplacement';
elseif strcmpi(ParamType, 'MeanDisplacementPerSecond') 
    paramName = 'MeanDisplacementPerSecond';
    seName = 'StdDisplacementPerSecond';
elseif strcmpi(ParamType, 'NCDivisionInfo') 
    paramName = 'NCDivisionInfo';
    seName = 'DivisionStdInfo';
else
    error(['Invalid choice of ParamType: ', ParamType])
end


params = this.EmbryoStats.(paramName);
if ~isempty(seName)
    paramSEs = this.EmbryoStats.(seName);
else
    paramSEs = NaN(size(params));
end
        