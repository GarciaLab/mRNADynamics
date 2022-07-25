function this = AddPerNucleusActivationEnergies(this)

NumExperiments = length(this.ExperimentPrefixes);
APResolution = this.Experiments{1}.APResolution;
NumAPbins = uint16(1/APResolution)+1;

parameters = {'CycleDurations', 'TimeOns', 'ElongationTimes','ElongationRates',...
    'TranscriptionWindows', 'PostTranscriptionDurations', 'LoadingRates', 'MaxFluos', 'PlateauHeights'};
traceNames = {'AnaphaseAligned', 'AnaphaseAligned3D', 'Tbinned', 'Tbinned3D',...
    'Unaligned', 'Unaligned3D'};

for paramIndex = 1:length(parameters)
    this.PerNucleusActivationEnergyFits.(parameters{paramIndex}) = {};
    this.PerNucleusActivationEnergyFits.(parameters{paramIndex}).ActivationEnergies = {};
    this.PerNucleusActivationEnergyFits.(parameters{paramIndex}).LogScalingCoefficients = {};
    this.PerNucleusActivationEnergyFits.(parameters{paramIndex}).Fits = {};
    if ~strcmpi(parameters{paramIndex}, 'CycleDurations')
        
        for traceIndex = 1:length(traceNames)
            
            this.PerNucleusActivationEnergyFits.(parameters{paramIndex}).ActivationEnergies.(traceNames{traceIndex}) =...
                NaN(NumAPbins, 6);
            this.PerNucleusActivationEnergyFits.(parameters{paramIndex}).ActivationEnergies.([traceNames{traceIndex}, 'StdError']) =...
                NaN(NumAPbins, 6);
            this.PerNucleusActivationEnergyFits.(parameters{paramIndex}).ActivationEnergies.([traceNames{traceIndex}, 'CIhigh']) =...
                NaN(NumAPbins, 6);
            this.PerNucleusActivationEnergyFits.(parameters{paramIndex}).ActivationEnergies.([traceNames{traceIndex}, 'CIlow']) =...
                NaN(NumAPbins, 6);
            
            
            this.PerNucleusActivationEnergyFits.(parameters{paramIndex}).LogScalingCoefficients.(traceNames{traceIndex}) =...
                NaN(NumAPbins, 6);
            this.PerNucleusActivationEnergyFits.(parameters{paramIndex}).LogScalingCoefficients.([traceNames{traceIndex}, 'StdError']) =...
                NaN(NumAPbins, 6);
            this.PerNucleusActivationEnergyFits.(parameters{paramIndex}).LogScalingCoefficients.([traceNames{traceIndex}, 'CIhigh']) =...
                NaN(NumAPbins, 6);
            this.PerNucleusActivationEnergyFits.(parameters{paramIndex}).LogScalingCoefficients.([traceNames{traceIndex}, 'CIlow']) =...
                NaN(NumAPbins, 6);
            
            this.PerNucleusActivationEnergyFits.(parameters{paramIndex}).Fits.(traceNames{traceIndex}) =...
                cell(NumAPbins, 6);
            
            this.PerNucleusActivationEnergyFits.(parameters{paramIndex}).R2s.([traceNames{traceIndex}]) =...
                NaN(NumAPbins, 6);
            
            
        end
        
    else
        
        this.PerNucleusActivationEnergyFits.(parameters{paramIndex}).ActivationEnergies.General =...
            NaN(1, 6);
        this.PerNucleusActivationEnergyFits.(parameters{paramIndex}).ActivationEnergies.GeneralStdError =...
            NaN(1, 6);
        this.PerNucleusActivationEnergyFits.(parameters{paramIndex}).ActivationEnergies.GeneralCIhigh =...
            NaN(1, 6);
        this.PerNucleusActivationEnergyFits.(parameters{paramIndex}).ActivationEnergies.GeneralCIlow = ...
            NaN(1, 6);
        
        this.PerNucleusActivationEnergyFits.(parameters{paramIndex}).LogScalingCoefficients.General =...
            NaN(1, 6);
        this.PerNucleusActivationEnergyFits.(parameters{paramIndex}).LogScalingCoefficients.GeneralStdError =...
            NaN(1, 6);
        this.PerNucleusActivationEnergyFits.(parameters{paramIndex}).LogScalingCoefficients.GeneralCIhigh =...
            NaN(1, 6);
        this.PerNucleusActivationEnergyFits.(parameters{paramIndex}).LogScalingCoefficients.GeneralCIlow =...
            NaN(1, 6);
        
        this.PerNucleusActivationEnergyFits.(parameters{paramIndex}).Fits.General = cell(1, 6);
        
        this.PerNucleusActivationEnergyFits.(parameters{paramIndex}).R2s.General = NaN(1, 6);
    end
end


for paramIndex = 1:length(parameters)
    parameter = parameters{paramIndex};
    for NC=9:14
        
        if strcmpi(parameter, 'CycleDurations')
            %disp([parameter,', NC: ', num2str(NC)])
            this = FitLTMPerNucleusActivationEnergy(this, parameter,  NC);
        else
            for APindex=1:NumAPbins
                for traceIndex = 1:length(traceNames)
                    TraceType = traceNames{traceIndex};
                    %disp([parameter,', ', TraceType,', NC: ', num2str(NC),', APbin: ', num2str(APindex)])
                    this = FitLTMPerNucleusActivationEnergy(this, parameter,  NC, APindex, TraceType);
                end
                
            end
        end
    end
end