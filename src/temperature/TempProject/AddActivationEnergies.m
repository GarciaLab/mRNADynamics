function this = AddActivationEnergies(this)

NumExperiments = length(this.ExperimentPrefixes);
APResolution = this.Experiments{1}.APResolution;
NumAPbins = uint16(1/APResolution)+1;

parameters = {'CycleDurations', 'TimeOns', 'ElongationTimes','ElongationRates',...
    'TranscriptionWindows', 'PostTranscriptionDurations', 'LoadingRates'};
traceNames = {'AnaphaseAligned', 'AnaphaseAligned3D', 'Tbinned', 'Tbinned3D',...
    'Unaligned', 'Unaligned3D'};

for paramIndex = 1:length(parameters)
    this.ActivationEnergyFits.(parameters{paramIndex}) = {};
    this.ActivationEnergyFits.(parameters{paramIndex}).ActivationEnergies = {};
    this.ActivationEnergyFits.(parameters{paramIndex}).LogScalingCoefficients = {};
    this.ActivationEnergyFits.(parameters{paramIndex}).Fits = {};
    if ~strcmpi(parameters{paramIndex}, 'CycleDurations')
        
        for traceIndex = 1:length(traceNames)
            
            this.ActivationEnergyFits.(parameters{paramIndex}).ActivationEnergies.(traceNames{traceIndex}) =...
                NaN(NumAPbins, 6);
            this.ActivationEnergyFits.(parameters{paramIndex}).ActivationEnergies.([traceNames{traceIndex}, 'StdError']) =...
                NaN(NumAPbins, 6);
            this.ActivationEnergyFits.(parameters{paramIndex}).ActivationEnergies.([traceNames{traceIndex}, 'CIhigh']) =...
                NaN(NumAPbins, 6);
            this.ActivationEnergyFits.(parameters{paramIndex}).ActivationEnergies.([traceNames{traceIndex}, 'CIlow']) =...
                NaN(NumAPbins, 6);
            
            
            this.ActivationEnergyFits.(parameters{paramIndex}).LogScalingCoefficients.(traceNames{traceIndex}) =...
                NaN(NumAPbins, 6);
            this.ActivationEnergyFits.(parameters{paramIndex}).LogScalingCoefficients.([traceNames{traceIndex}, 'StdError']) =...
                NaN(NumAPbins, 6);
            this.ActivationEnergyFits.(parameters{paramIndex}).LogScalingCoefficients.([traceNames{traceIndex}, 'CIhigh']) =...
                NaN(NumAPbins, 6);
            this.ActivationEnergyFits.(parameters{paramIndex}).LogScalingCoefficients.([traceNames{traceIndex}, 'CIlow']) =...
                NaN(NumAPbins, 6);
            
            this.ActivationEnergyFits.(parameters{paramIndex}).Fits.(traceNames{traceIndex}) =...
                cell(NumAPbins, 6);
            
            this.ActivationEnergyFits.(parameters{paramIndex}).R2s.([traceNames{traceIndex}]) =...
                NaN(NumAPbins, 6);
            
            
        end
        
    else
        
        this.ActivationEnergyFits.(parameters{paramIndex}).ActivationEnergies.General =...
            NaN(1, 6);
        this.ActivationEnergyFits.(parameters{paramIndex}).ActivationEnergies.GeneralStdError =...
            NaN(1, 6);
        this.ActivationEnergyFits.(parameters{paramIndex}).ActivationEnergies.GeneralCIhigh =...
            NaN(1, 6);
        this.ActivationEnergyFits.(parameters{paramIndex}).ActivationEnergies.GeneralCIlow = ...
            NaN(1, 6);
        
        this.ActivationEnergyFits.(parameters{paramIndex}).LogScalingCoefficients.General =...
            NaN(1, 6);
        this.ActivationEnergyFits.(parameters{paramIndex}).LogScalingCoefficients.GeneralStdError =...
            NaN(1, 6);
        this.ActivationEnergyFits.(parameters{paramIndex}).LogScalingCoefficients.GeneralCIhigh =...
            NaN(1, 6);
        this.ActivationEnergyFits.(parameters{paramIndex}).LogScalingCoefficients.GeneralCIlow =...
            NaN(1, 6);
        
        this.ActivationEnergyFits.(parameters{paramIndex}).Fits.General = cell(1, 6);
        
        this.ActivationEnergyFits.(parameters{paramIndex}).R2s.General = NaN(1, 6);
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
                    this = FitLTMActivationEnergy(this, parameter,  NC, APindex, TraceType);
                end
                
            end
        end
    end
end