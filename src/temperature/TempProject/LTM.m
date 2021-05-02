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
        
        LegendLabels = {};
        
        APLengths = [];
        
        MeanProfiles = {};
        
        GeneLength = [];
        
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
        
        UniqueTemperatures = [];
        AutocorrElongationTimes = [];
        TimingCoeffs = [];
        FluoCoeffs = [];
        FluoAdjustedActivationEnergyFits = {};
        
        AutocorrElongationTimeInfo = {};
        FluoScalingInfo = {};
        TimeScalingInfo = {};
        
        
    end
    
    
    
    methods
        %% Constructors
        
        
        function this = LTM(ProjectName, ProjectType, GeneLength, FluoParamPath, IncludedNCs, time_delta,...
                MinimumTraceCount, MinimumTimePoints, MinimumBinCount, MinimumFittingPoints, MinimumSchnitzCount)
            %liveProject Construct an instance of this class
            %   Detailed explanation goes here
            this.ProjectName = ProjectName;
            
            if exist('ProjectType', 'var')
                this.ProjectType = ProjectType;
            end
            
            if exist('IncludedNCs', 'var')
                this.IncludedNCs = IncludedNCs;
            else
                this.IncludedNCs =[10, 11, 12, 13, 14];
            end
            
            if exist('time_delta', 'var')
                this.time_delta = time_delta;
            else
                this.time_delta = 30; % unit: minutes
            end
            
            if exist('MinimumTraceCount', 'var')
                this.MinimumTraceCount = MinimumTraceCount;
            else
                this.MinimumTraceCount = 5;
            end
            if exist('MinimumTimePoints', 'var')
                this.MinimumTimePoints = MinimumTimePoints;
            else
                this.MinimumTimePoints = 5;
            end
            if exist('MinimumFittingPoints', 'var')
                this.MinimumFittingPoints = MinimumFittingPoints;
            else
                this.MinimumFittingPoints = 5;
            end
            if exist('MinimumBinCount', 'var')
                this.MinimumBinCount = MinimumBinCount;
            else
                this.MinimumBinCount = 2;
            end
            if exist('MinimumSchnitzCount', 'var')
                this.MinimumSchnitzCount = MinimumSchnitzCount;
            else
                this.MinimumSchnitzCount = 5;
            end
            if exist('R2bound', 'var')
                this.R2bound = R2bound;
            else
                this.R2bound = 0.95;
            end
            if exist('alpha', 'var')
                this.alpha = alpha;
            else
                this.alpha = 0.95;
            end
            
            if exist('GeneLength', 'var')
                this.GeneLength = GeneLength;
            else
                this.GeneLength = 5858;
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
            this = AddMeanProfiles(this);
            this = AddMeanInitiationRates(this);
            this = AddFractionOns(this);
            this = AddHealthInfo(this);
            this = CalculateMeanSpotFluo(this);
            this = AddActivationEnergies(this);
            
            
        end
        
        
        
        %% Methods
        
        
        
        
        
        
    end
end