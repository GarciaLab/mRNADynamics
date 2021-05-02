function this = AddFluoAdjustedActivationEnergies(this)

NumExperiments = length(this.ExperimentPrefixes);
APResolution = this.Experiments{1}.APResolution;
NumAPbins = uint16(1/APResolution)+1;

parameters = {'LoadingRates', 'MaxFluos', 'PlateauHeights'};
traceNames = {'AnaphaseAligned', 'AnaphaseAligned3D', 'Tbinned', 'Tbinned3D',...
    'Unaligned', 'Unaligned3D'};

for paramIndex = 1:length(parameters)
    this.FluoAdjustedActivationEnergyFits.(parameters{paramIndex}) = {};
    this.FluoAdjustedActivationEnergyFits.(parameters{paramIndex}).ActivationEnergies = {};
    this.FluoAdjustedActivationEnergyFits.(parameters{paramIndex}).LogScalingCoefficients = {};
    this.FluoAdjustedActivationEnergyFits.(parameters{paramIndex}).Fits = {};
    if ~strcmpi(parameters{paramIndex}, 'CycleDurations')
        
        for traceIndex = 1:length(traceNames)
            
            this.FluoAdjustedActivationEnergyFits.(parameters{paramIndex}).ActivationEnergies.(traceNames{traceIndex}) =...
                NaN(NumAPbins, 6);
            this.FluoAdjustedActivationEnergyFits.(parameters{paramIndex}).ActivationEnergies.([traceNames{traceIndex}, 'StdError']) =...
                NaN(NumAPbins, 6);
            this.FluoAdjustedActivationEnergyFits.(parameters{paramIndex}).ActivationEnergies.([traceNames{traceIndex}, 'CIhigh']) =...
                NaN(NumAPbins, 6);
            this.FluoAdjustedActivationEnergyFits.(parameters{paramIndex}).ActivationEnergies.([traceNames{traceIndex}, 'CIlow']) =...
                NaN(NumAPbins, 6);
            
            
            this.FluoAdjustedActivationEnergyFits.(parameters{paramIndex}).LogScalingCoefficients.(traceNames{traceIndex}) =...
                NaN(NumAPbins, 6);
            this.FluoAdjustedActivationEnergyFits.(parameters{paramIndex}).LogScalingCoefficients.([traceNames{traceIndex}, 'StdError']) =...
                NaN(NumAPbins, 6);
            this.FluoAdjustedActivationEnergyFits.(parameters{paramIndex}).LogScalingCoefficients.([traceNames{traceIndex}, 'CIhigh']) =...
                NaN(NumAPbins, 6);
            this.FluoAdjustedActivationEnergyFits.(parameters{paramIndex}).LogScalingCoefficients.([traceNames{traceIndex}, 'CIlow']) =...
                NaN(NumAPbins, 6);
            
            this.FluoAdjustedActivationEnergyFits.(parameters{paramIndex}).Fits.(traceNames{traceIndex}) =...
                cell(NumAPbins, 6);
            
            this.FluoAdjustedActivationEnergyFits.(parameters{paramIndex}).R2s.([traceNames{traceIndex}]) =...
                NaN(NumAPbins, 6);
            
            
        end
        
    else
        
        this.FluoAdjustedActivationEnergyFits.(parameters{paramIndex}).ActivationEnergies.General =...
            NaN(1, 6);
        this.FluoAdjustedActivationEnergyFits.(parameters{paramIndex}).ActivationEnergies.GeneralStdError =...
            NaN(1, 6);
        this.FluoAdjustedActivationEnergyFits.(parameters{paramIndex}).ActivationEnergies.GeneralCIhigh =...
            NaN(1, 6);
        this.FluoAdjustedActivationEnergyFits.(parameters{paramIndex}).ActivationEnergies.GeneralCIlow = ...
            NaN(1, 6);
        
        this.FluoAdjustedActivationEnergyFits.(parameters{paramIndex}).LogScalingCoefficients.General =...
            NaN(1, 6);
        this.FluoAdjustedActivationEnergyFits.(parameters{paramIndex}).LogScalingCoefficients.GeneralStdError =...
            NaN(1, 6);
        this.FluoAdjustedActivationEnergyFits.(parameters{paramIndex}).LogScalingCoefficients.GeneralCIhigh =...
            NaN(1, 6);
        this.FluoAdjustedActivationEnergyFits.(parameters{paramIndex}).LogScalingCoefficients.GeneralCIlow =...
            NaN(1, 6);
        
        this.FluoAdjustedActivationEnergyFits.(parameters{paramIndex}).Fits.General = cell(1, 6);
        
        this.FluoAdjustedActivationEnergyFits.(parameters{paramIndex}).R2s.General = NaN(1, 6);
    end
end


for paramIndex = 1:length(parameters)
    parameter = parameters{paramIndex};
    for NC=9:14
        
        if strcmpi(parameter, 'CycleDurations')
            %disp([parameter,', NC: ', num2str(NC)])
            this = FitLTMActivationEnergy(this, parameter,  NC);
        else
            for APindex=1:NumAPbins
                for traceIndex = 1:length(traceNames)
                    TraceType = traceNames{traceIndex};
                    %disp([parameter,', ', TraceType,', NC: ', num2str(NC),', APbin: ', num2str(APindex)])
                    this = FitFluoAdjustedActivationEnergy(this, parameter,  NC, APindex, TraceType);
                end
                
            end
        end
    end
end