function this = AddHmmActivationEnergies(this, IncludeRescaling)
if ~exist('IncludeRescaling', 'var')
    IncludeRescaling   = true;
end


parameters = {'InitiationRates', 'Durations', 'Frequencies'};
NumAPbins = size(this.InitiationRates, 3);
NumTimeBins = size(this.InitiationRates, 2);
NumScaledTimeBins = size(this.ScaledInitiationRates, 2);
for paramIndex = 1:length(parameters)
    
    this.ActivationEnergyFits.(parameters{paramIndex}) = {};
    this.ActivationEnergyFits.(parameters{paramIndex}).ActivationEnergies = NaN(NumTimeBins, NumAPbins);
    this.ActivationEnergyFits.(parameters{paramIndex}).SEActivationEnergies = NaN(NumTimeBins, NumAPbins);
    this.ActivationEnergyFits.(parameters{paramIndex}).LogScalingCoefficients = NaN(NumTimeBins, NumAPbins);
    this.ActivationEnergyFits.(parameters{paramIndex}).SELogScalingCoefficients = NaN(NumTimeBins, NumAPbins);
    this.ActivationEnergyFits.(parameters{paramIndex}).R2s = NaN(NumTimeBins, NumAPbins);
    this.ActivationEnergyFits.(parameters{paramIndex}).Fits = cell(NumTimeBins, NumAPbins);
    
    if IncludeRescaling
        this.FluoScaledActivationEnergyFits.(parameters{paramIndex}) = {};
        this.FluoScaledActivationEnergyFits.(parameters{paramIndex}).ActivationEnergies = NaN(NumTimeBins, NumAPbins);
        this.FluoScaledActivationEnergyFits.(parameters{paramIndex}).SEActivationEnergies = NaN(NumTimeBins, NumAPbins);
        this.FluoScaledActivationEnergyFits.(parameters{paramIndex}).LogScalingCoefficients = NaN(NumTimeBins, NumAPbins);
        this.FluoScaledActivationEnergyFits.(parameters{paramIndex}).SELogScalingCoefficients = NaN(NumTimeBins, NumAPbins);
        this.FluoScaledActivationEnergyFits.(parameters{paramIndex}).R2s = NaN(NumTimeBins, NumAPbins);
        this.FluoScaledActivationEnergyFits.(parameters{paramIndex}).Fits = cell(NumTimeBins, NumAPbins);
        
        this.TimeScaledActivationEnergyFits.(parameters{paramIndex}) = {};
        this.TimeScaledActivationEnergyFits.(parameters{paramIndex}).ActivationEnergies = NaN(NumScaledTimeBins, NumAPbins);
        this.TimeScaledActivationEnergyFits.(parameters{paramIndex}).SEActivationEnergies = NaN(NumScaledTimeBins, NumAPbins);
        this.TimeScaledActivationEnergyFits.(parameters{paramIndex}).LogScalingCoefficients = NaN(NumScaledTimeBins, NumAPbins);
        this.TimeScaledActivationEnergyFits.(parameters{paramIndex}).SELogScalingCoefficients = NaN(NumScaledTimeBins, NumAPbins);
        this.TimeScaledActivationEnergyFits.(parameters{paramIndex}).R2s = NaN(NumScaledTimeBins, NumAPbins);
        this.TimeScaledActivationEnergyFits.(parameters{paramIndex}).Fits = cell(NumScaledTimeBins, NumAPbins);
        
        this.DoubleScaledActivationEnergyFits.(parameters{paramIndex}) = {};
        this.DoubleScaledActivationEnergyFits.(parameters{paramIndex}).ActivationEnergies = NaN(NumScaledTimeBins, NumAPbins);
        this.DoubleScaledActivationEnergyFits.(parameters{paramIndex}).SEActivationEnergies = NaN(NumScaledTimeBins, NumAPbins);
        this.DoubleScaledActivationEnergyFits.(parameters{paramIndex}).LogScalingCoefficients = NaN(NumScaledTimeBins, NumAPbins);
        this.DoubleScaledActivationEnergyFits.(parameters{paramIndex}).SELogScalingCoefficients = NaN(NumScaledTimeBins, NumAPbins);
        this.DoubleScaledActivationEnergyFits.(parameters{paramIndex}).R2s = NaN(NumScaledTimeBins, NumAPbins);
        this.DoubleScaledActivationEnergyFits.(parameters{paramIndex}).Fits = cell(NumScaledTimeBins, NumAPbins);
    end
end


for paramIndex = 1:length(parameters)
    parameter = parameters{paramIndex};
    
    for TimeIndex = 1:NumTimeBins
        
        for APindex=1:NumAPbins
            
            this = FitHMMActivationEnergy(this, parameter, false, false,  TimeIndex, APindex);
            if IncludeRescaling
                this = FitHMMActivationEnergy(this, parameter, true, false,  TimeIndex, APindex);
            end
        end
        
    end
    if IncludeRescaling
        for TimeIndex = 1:NumScaledTimeBins
            
            for APindex=1:NumAPbins
                
                this = FitHMMActivationEnergy(this, parameter, false, true,  TimeIndex, APindex);
                this = FitHMMActivationEnergy(this, parameter, true, true,  TimeIndex, APindex);
            end
            
        end
    end
end