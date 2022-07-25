classdef LTM
    %liveProject Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        ProjectName = '';
        ProjectType = '';
        ExperimentTypes = {};
        ExperimentPrefixes = {};
        Experiments = {};
        ExperimentStatuses = {};
        IncludedExperiments = [];
        ProcessedExperiments = [];
        
        IncludedNCs = [];
        Regions = {};
        Temp_sps = [];
        Temp_obs = [];
        MinAPs = [];
        MaxAPs = [];
        
        R = [];
        alpha = [];
        R2bound = [];
        time_delta = [];
        MinimumTraceCount = [];
        MinimumTimePoints = [];
        MinimumFittingPoints = [];
        MinimumBinCount = [];
        MinimumSchnitzCount = [];
        MinimumEmbryos = [];
        
        LegendLabels = {};
        
        APLengths = [];
        
        MeanProfiles = {};
        BinnedMeanProfiles = {};

        
        GeneLength = [];
        
        BinnedProfileParameters = {};
        PerNucleusProfileParameters = {};
        BinnedPerNucleusProfileParameters = {};
        MeanInitiationRates = {};
        TimeOns = {};
        TimeOffs = {};
        UnloadingRates = {};
        ElongationTimes = {};
        MeanFitR2s = {};
        
        SchnitzCounts = [];
        FractionOns = [];
        

        CompiledParameters = {};
        
        Fits = {};
        
        EmbryoStats = {};
        
        MeanSpotFluo = {};
        
        ActivationEnergyFits = {};
        BinnedActivationEnergyFits = {};
        BinnedPerNucleusActivationEnergyFits = {};
        PerNucleusActivationEnergyFits = {};
        FluoAdjustedBinnedActivationEnergyFits = {};
        FluoAdjustedBinnedPerNucleusActivationEnergyFits = {};
        FluoAdjustedPerNucleusActivationEnergyFits = {};
        
        UniqueTemperatures = [];
        AutocorrElongationTimes = [];
        TimingCoeffs = [];
        FluoCoeffs = [];
        FluoAdjustedActivationEnergyFits = {};
        
        AutocorrElongationTimeInfo = {};
        FluoScalingInfo = {};
        TimeScalingInfo = {};
        
        UseManualApproval = false;
        
        
    end
    
    
    
    methods
        %% Constructors
        
        
        function this = LTM(ProjectName, ProjectType, GeneLength, FluoParamPath, IncludedNCs, time_delta,...
                MinimumTraceCount, MinimumTimePoints, MinimumBinCount, MinimumFittingPoints, MinimumSchnitzCount,...
                MinimumEmbryos, UseManualApproval)
            %liveProject Construct an instance of this class
            %   Detailed explanation goes here
            this.ProjectName = ProjectName;
            
            if exist('ProjectType', 'var')
                this.ProjectType = ProjectType;
            end
            
            if exist('IncludedNCs', 'var') & ~isempty(IncludedNCs)
                this.IncludedNCs = IncludedNCs;
            else
                this.IncludedNCs =[10, 11, 12, 13, 14];
            end
            
            if exist('time_delta', 'var') & ~isempty(time_delta)
                this.time_delta = time_delta;
            else
                this.time_delta = 30; % unit: minutes
            end
            
            if exist('MinimumTraceCount', 'var') & ~isempty(MinimumTraceCount)
                this.MinimumTraceCount = MinimumTraceCount;
            else
                this.MinimumTraceCount = 5;
            end
            
            
            
            if exist('MinimumTimePoints', 'var') & ~isempty(MinimumTimePoints)
                this.MinimumTimePoints = MinimumTimePoints;
            else
                this.MinimumTimePoints = 4;
            end
            if exist('MinimumFittingPoints', 'var') & ~isempty(MinimumFittingPoints)
                this.MinimumFittingPoints = MinimumFittingPoints;
            else
                this.MinimumFittingPoints = 5;
            end
            if exist('MinimumBinCount', 'var')& ~isempty(MinimumBinCount)
                this.MinimumBinCount = MinimumBinCount;
            else
                this.MinimumBinCount = 1;
            end
            if exist('MinimumSchnitzCount', 'var') & ~isempty(MinimumSchnitzCount)
                this.MinimumSchnitzCount = MinimumSchnitzCount;
            else
                this.MinimumSchnitzCount = 5;
            end
            
            if exist('MinimumEmbryos', 'var') & ~isempty(MinimumEmbryos)
                this.MinimumEmbryos = MinimumEmbryos;
            else
                this.MinimumEmbryos = 2;
            end
            if exist('R2bound', 'var') & ~isempty(R2bound)
                this.R2bound = R2bound;
            else
                this.R2bound = 0.95;
            end
            if exist('alpha', 'var')  & ~isempty(alpha)
                this.alpha = alpha;
            else
                this.alpha = 0.95;
            end
            
            if exist('GeneLength', 'var') & ~isempty(GeneLength)
                this.GeneLength = GeneLength;
            else
                this.GeneLength = 5858;
            end
            
            if exist('UseManualApproval', 'var') & ~isempty(UseManualApproval)
                this.UseManualApproval = UseManualApproval;
            else
                this.UseManualApproval = false;
            end
            
            
            this.R = 8.314*10^(-3); % kJ * K^(-1)*mol^(-1)
            
            this.ExperimentPrefixes = getProjectPrefixes(ProjectName);
            
            for i = 1:length(this.ExperimentPrefixes)
                Prefix = this.ExperimentPrefixes{i};
                this.Experiments{i} = LiveExperiment(Prefix);
                this.ExperimentTypes{i} = this.Experiments{i}.experimentType;
                this.ExperimentStatuses{i} = TemperatureMS2Status(ProjectName, Prefix);
                include_set(i) = this.ExperimentStatuses{i}.include_set;
                if exist('ProjectType', 'var')
                    if lower(this.ExperimentTypes{i}) ~= lower(ProjectType)
                        include_set(i) = 0;
                    end
                end
                if this.ExperimentStatuses{i}.hasCompiledParticles
                    finished_set(i) = 1;
                else
                    finished_set(i) = 0;
                end
                this.Regions{i} = this.ExperimentStatuses{i}.Region;
                this.Temp_sps(i) = this.ExperimentStatuses{i}.Temp_sp;
                this.Temp_obs(i) = this.ExperimentStatuses{i}.Temp_obs;
                %                 if ~isempty(this.ExperimentStatuses{i}.MinAP)
                %                     this.MinAPs(i) = this.ExperimentStatuses{i}.MinAP;
                %                 else
                %                     this.MinAPs(i) = NaN;
                %                 end
                %                 if ~isempty(this.ExperimentStatuses{i}.MaxAP)
                %                     this.MaxAPs(i) = this.ExperimentStatuses{i}.MaxAP;
                %                 else
                %                     this.MaxAPs(i) = NaN;
                %                 end
            end
            
            this.IncludedExperiments = find(include_set > 0);
            this.ProcessedExperiments = find((include_set > 0) & (finished_set > 0));
            
            if exist('FluoParamPath', 'var')
                this = AddFluoCoeffs(this, FluoParamPath);
            end
            
            this.APLengths = AddAPLengths(this);
%             this = AddMeanProfiles(this);
%             
%             this = AddBinnedMeanProfiles(this);
%             this = AddMeanInitiationRates(this);
%             this = AddMeanPerNucleusInitiationRates(this);
%             this = AddTBinnedMeanInitiationRates(this);
%             this = AddTBinnedMeanPerNucleusInitiationRates(this);
%             this = AddHealthInfo(this);
%             this = CalculateMeanSpotFluo(this);
%             TimingParamPath = 'S:/Gabriella\Dropbox\TemperatureParameters\hbBAC-MS2-V3\';
%             this = AddTimingCoeffs(this, TimingParamPath);
%             this = AddActivationEnergies(this);
%             this = AddBinnedActivationEnergies(this);
%             this = AddPerNucleusActivationEnergies(this);
%             this = AddBinnedPerNucleusActivationEnergies(this);
        end
        
        
        
        %% Methods
        
        
        
        
        
        
    end
end