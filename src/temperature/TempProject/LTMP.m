classdef LTMP
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
        
        LegendLabels = {};
        APLengths = [];
        
        IncludedNCs = [];
        Regions = {};
        Temp_sps = [];
        Temp_obs = [];
        MinAPs = [];
        MaxAPs = [];
        
        time_delta = [];
        MinimumTimePoints = [];
        MinimumNuclearCount = []; % From LTP
        MinimumTraceCount = [];  % From LTM
        MinimumFittingPoints = []; % From LTM
        MinimumBinCount = []; % From LTM
        MinimumSchnitzCount = []; % From LTM
        
        FluoScalingFactors = [];
        FluoScalingFactorSEs = [];
        
        
        % LTP specific parameters for nuclear fluorescence offset 
        MinOffsetAP = [];
        MaxOffsetAP = [];
        
        % LTM specific parameters 
        R = [];
        alpha = [];
        R2bound = [];
        GeneLength = [];
        
         %% LTM parameters start
        MeanProfiles = {};
        MeanInitiationRates = {};
        TimeOns = {};
        TimeOffs = {};
        UnloadingRates = {};
        ElongationTimes = {};
        MeanFitR2s = {};

        SchnitzCounts = [];
        FractionOns = [];

        Fits = {};
        
        EmbryoStats = {};
        
        MeanSpotFluo = {};
        
        ActivationEnergyFits = {};
        
        % LTM Parameters end
        
        %% LTP parameters start 
        
        TFProfiles = {};
        TempFluoOffsets = {};
        SetFluoOffsets = {};

        
 
        
        
    end
    
    
    
    methods
        %% Constructors
        
        
        function this = LTMP(ProjectName, ProjectType, GeneLength, IncludedNCs, time_delta,...
                MinimumNuclearCount, MinimumTraceCount, MinimumTimePoints, MinimumBinCount, MinimumFittingPoints, MinimumSchnitzCount)
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
            
            if exist('MinimumNuclearCount', 'var')
                this.MinimumNuclearCount = MinimumNuclearCount;
            else
                this.MinimumNuclearCount = 5;
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
                this.GeneLength = 7675;
            end
                
            
            this.R = 8.314*10^(-3); % kJ * K^(-1)*mol^(-1)
            
            AllPrefixes = getProjectPrefixes(ProjectName);
            TempVector = zeros(1, length(AllPrefixes));
            for i = 1:length(AllPrefixes)
                TempVector(i) = pullTempSpFromDataStatus(ProjectName, AllPrefixes{i});
            end
            [~, indexorder] = sort(TempVector, 'descend');
            this.ExperimentPrefixes = AllPrefixes(indexorder);
            

            for i = 1:length(this.ExperimentPrefixes)
                Prefix = this.ExperimentPrefixes{i};
                this.Experiments{i} = LiveExperiment(Prefix);
                this.ExperimentTypes{i} = this.Experiments{i}.experimentType;
                disp(num2str(i))
                this.ExperimentStatuses{i} = TemperatureMS2JB3Status(ProjectName, Prefix);
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
            
            if isempty(this.MinOffsetAP)
                if contains(lower(this.ProjectName), 'hb')
                    this.MinOffsetAP = 0.575;
                end
            end
            if isempty(this.MaxOffsetAP)
                if contains(lower(this.ProjectName), 'hb')
                    this.MaxOffsetAP = 0.675;
                end
            end
            
            this = AddTFProfiles(this);
            this = AddNBFluoOffsets(this);
            
            
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