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
                this.IncludedNCs =[9, 10, 11, 12, 13, 14];
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
            this.MeanProfiles = AddMeanProfiles(this);
            this = AddMeanInitiationRates(this);
            this = AddFractionOns(this);
            this = AddHealthInfo(this);
            this = CalculateMeanSpotFluo(this, this.time_delta);
            this = AddActivationEnergies(this);
        end
        
        
        
        %% Methods
        

        
  
        
       
        function this = AddFractionOns(this)
            % author: G. Martini
            % date created: 1/19/21
            % date last modified: 1/19/21
            NumExperiments = length(this.ExperimentPrefixes);
            APResolution = this.Experiments{1}.APResolution;
            NumAPbins = uint16(1/APResolution)+1;
            APboundaries = 0:APResolution:1;
            APlowers = APboundaries(1:end-1);
            APuppers = APboundaries(2:end);
            
            
            this.SchnitzCounts = zeros(NumExperiments, NumAPbins, 6);
            this.FractionOns = NaN(NumExperiments, NumAPbins, 6);
            disp('Adding Fraction Ons.')
            for SetIndex = 1:length(this.ExperimentPrefixes)
                %disp(num2str(SetIndex))
                if ~ismember(SetIndex, this.ProcessedExperiments)
                    continue
                end
                schnitzcells = getSchnitzcells(this.Experiments{SetIndex});
                Particles = getParticles(this.Experiments{SetIndex});
                if iscell(Particles) & (length(Particles) == 1)
                    Particles = Particles{1};
                end
                load([this.Experiments{SetIndex}.resultsFolder, 'CompiledParticles.mat'], 'CompiledParticles')
                if iscell(CompiledParticles) & (length(CompiledParticles) == 1)
                    CompiledParticles = CompiledParticles{1};
                end
                
                CompiledParticleSchnitzMap = [CompiledParticles(:).Nucleus];
                if (length(unique(CompiledParticleSchnitzMap)) ~= length(CompiledParticleSchnitzMap)) | ...
                        (length(CompiledParticleSchnitzMap) ~= length(CompiledParticles))
                    %disp(num2str(SetIndex))
                    %warning('Compiled Particle mapping to schnitz cells either contains duplicates or unmapped particles.')
                    CompiledParticles = ...
                        correctParticleMappingsToSchnitzCells(this.ExperimentPrefixes{SetIndex});
                    if iscell(CompiledParticles) & (length(CompiledParticles) == 1)
                        CompiledParticles = CompiledParticles{1};
                    end
                    CompiledParticleSchnitzMap = [CompiledParticles(:).Nucleus];
                    
                end
                FrameInfo = getFrameInfo(this.Experiments{SetIndex});
                PixelSize = this.Experiments{SetIndex}.pixelSize_um;
                nucleusDiameters = zeros(1, 6);
                for nc=9:14
                    nucleusDiameters(nc-8) = getDefaultParameters(FrameInfo,[['d', num2str(nc)]])/PixelSize; % in pixels
                end
                xDim = this.Experiments{SetIndex}.xDim;
                yDim = this.Experiments{SetIndex}.yDim;
                zDim = this.Experiments{SetIndex}.zDim;
                
                SchnitzApproved = zeros(1, length(schnitzcells));
                HasAllFrames = zeros(1, length(schnitzcells));
                SchnitzInBoundsAllFrames = zeros(1, length(schnitzcells));
                SchnitzCycles = zeros(1, length(schnitzcells));
                SchnitzAPs = zeros(1, length(schnitzcells));
                SchnitzHasValidParticle = zeros(1, length(schnitzcells));
                
                
                
                
                for sc_idx=1:length(schnitzcells)
                    if (schnitzcells(sc_idx).Approved > 0)  & (schnitzcells(sc_idx).Flag ~= 6)
                        SchnitzApproved(sc_idx) = 1;
                    end
                    if schnitzcells(sc_idx).VelocityInfo.SchnitzHasAllFrames
                        HasAllFrames(sc_idx) =  1;
                        nc_idx = schnitzcells(sc_idx).cycle-8;
                        sc_d = nucleusDiameters(nc_idx);
                        IsInBounds = all((schnitzcells(sc_idx).cenx >= sc_d/2) & ...
                            (schnitzcells(sc_idx).cenx <= xDim-sc_d/2) & ...
                            (schnitzcells(sc_idx).ceny >= sc_d/2) & ...
                            (schnitzcells(sc_idx).ceny <= yDim-sc_d/2));
                        if IsInBounds
                            SchnitzInBoundsAllFrames(sc_idx) = 1;
                        end
                    end
                    
                    if ~isempty(schnitzcells(sc_idx).cycle)
                        SchnitzCycles(sc_idx) =  schnitzcells(sc_idx).cycle;
                    end
                    SchnitzAPs(sc_idx) = mean(schnitzcells(sc_idx).APpos);
                    cp_idx = find(CompiledParticleSchnitzMap == sc_idx);
                    if ~isempty(cp_idx)
                        CurrentCP = CompiledParticles(cp_idx);
                        NValidFrames = sum(CurrentCP.ManualFrameApproved(1:length(CurrentCP.Frame)));
                        if CurrentCP.ManualApproved & (NValidFrames >= this.MinimumTimePoints)
                            SchnitzHasValidParticle(sc_idx) = 1;
                        end
                    end
                end
                for NC = 9:14
                    for APidx = 1:NumAPbins-1
                        this.SchnitzCounts(SetIndex,APidx, NC-8)  = sum((SchnitzApproved == 1) &...
                            (SchnitzInBoundsAllFrames == 1) & (HasAllFrames == 1) & ...
                            (SchnitzCycles == NC) & (SchnitzAPs >= APlowers(APidx)) & (SchnitzAPs < APuppers(APidx)));
                        this.FractionOns(SetIndex,APidx, NC-8)  = sum((SchnitzApproved == 1) & ...
                            (SchnitzInBoundsAllFrames == 1) & (HasAllFrames == 1) & ...
                            (SchnitzCycles == NC) & (SchnitzAPs >= APlowers(APidx)) & (SchnitzAPs < APuppers(APidx)) & ...
                            (SchnitzHasValidParticle == 1))/this.SchnitzCounts(SetIndex,APidx, NC-8);
                    end
                end
                
            end
        end
        
        function this = AddHealthInfo(this)
            NumSets = length(this.Experiments);
            
            
            this.EmbryoStats.SchnitzCount = NaN(NumSets, 6);
            this.EmbryoStats.FractionSickNuclei = NaN(NumSets, 6);
            this.EmbryoStats.FractionRejectedNuclei = NaN(NumSets, 6);
            this.EmbryoStats.FractionCompleteNuclei = NaN(NumSets, 6);
            this.EmbryoStats.FractionFirstLastFrameNuclei = NaN(NumSets, 6);
            this.EmbryoStats.MeanTotalXDistanceTraveled = NaN(NumSets, 6);
            this.EmbryoStats.MeanTotalYDistanceTraveled = NaN(NumSets, 6);
            this.EmbryoStats.MeanTotalDistanceTraveled = NaN(NumSets, 6);
            this.EmbryoStats.MeanDistanceTraveledPerSecond = NaN(NumSets, 6);
            this.EmbryoStats.MeanTotalXDisplacement = NaN(NumSets, 6);
            this.EmbryoStats.MeanTotalYDisplacement = NaN(NumSets, 6);
            this.EmbryoStats.MeanTotalDisplacement = NaN(NumSets, 6);
            this.EmbryoStats.MeanDisplacementPerSecond = NaN(NumSets, 6);
            this.EmbryoStats.StdTotalXDistanceTraveled = NaN(NumSets, 6);
            this.EmbryoStats.StdTotalYDistanceTraveled = NaN(NumSets, 6);
            this.EmbryoStats.StdTotalDistanceTraveled = NaN(NumSets, 6);
            this.EmbryoStats.StdDistanceTraveledPerSecond = NaN(NumSets, 6);
            this.EmbryoStats.StdTotalXDisplacement = NaN(NumSets, 6);
            this.EmbryoStats.StdTotalYDisplacement = NaN(NumSets, 6);
            this.EmbryoStats.StdTotalDisplacement = NaN(NumSets, 6);
            this.EmbryoStats.StdDisplacementPerSecond = NaN(NumSets, 6);
            this.EmbryoStats.NCDivisionInfo = NaN(NumSets, 6);
            this.EmbryoStats.DivisionStdInfo = NaN(NumSets, 6);
            
            
            
            
            
            for SetIndex = 1:NumSets
                if ~ismember(SetIndex, this.ProcessedExperiments)
                    continue
                end
                
                liveExperiment = this.Experiments{SetIndex};
                HealthSummaryPath = [liveExperiment.resultsFolder, 'HealthSummary.mat'];
                load(HealthSummaryPath)
                for nc_idx = 1:length(this.IncludedNCs)
                    NC = this.IncludedNCs(nc_idx);
                    this.EmbryoStats.SchnitzCount(SetIndex, NC-8) = HealthSummary.SchnitzCount(nc_idx);
                    this.EmbryoStats.FractionSickNuclei(SetIndex, NC-8)  =HealthSummary.FractionSickNuclei(nc_idx);
                    this.EmbryoStats.FractionRejectedNuclei(SetIndex, NC-8)  = HealthSummary.FractionRejectedNuclei(nc_idx);
                    this.EmbryoStats.FractionCompleteNuclei(SetIndex, NC-8)  = HealthSummary.FractionCompleteNuclei(nc_idx);
                    this.EmbryoStats.FractionFirstLastFrameNuclei(SetIndex, NC-8)  = HealthSummary.FractionFirstLastFrameNuclei(nc_idx);
                    this.EmbryoStats.MeanTotalXDistanceTraveled(SetIndex, NC-8)   = HealthSummary.MeanTotalXDistanceTraveled(nc_idx);
                    this.EmbryoStats.MeanTotalYDistanceTraveled(SetIndex, NC-8)   = HealthSummary.MeanTotalYDistanceTraveled(nc_idx);
                    this.EmbryoStats.MeanTotalDistanceTraveled(SetIndex, NC-8)   = HealthSummary.MeanTotalDistanceTraveled(nc_idx);
                    this.EmbryoStats.MeanDistanceTraveledPerSecond(SetIndex, NC-8)   = HealthSummary.MeanDistanceTraveledPerSecond(nc_idx);
                    this.EmbryoStats.MeanTotalXDisplacement(SetIndex, NC-8)   = HealthSummary.MeanTotalXDisplacement(nc_idx);
                    this.EmbryoStats.MeanTotalYDisplacement(SetIndex, NC-8)   = HealthSummary.MeanTotalYDisplacement(nc_idx);
                    this.EmbryoStats.MeanTotalDisplacement(SetIndex, NC-8) = HealthSummary.MeanTotalDisplacement(nc_idx);
                    this.EmbryoStats.MeanDisplacementPerSecond(SetIndex, NC-8) = HealthSummary.MeanDisplacementPerSecond(nc_idx);
                    this.EmbryoStats.StdTotalXDistanceTraveled(SetIndex, NC-8) = HealthSummary.StdTotalXDistanceTraveled(nc_idx);
                    this.EmbryoStats.StdTotalYDistanceTraveled(SetIndex, NC-8) =HealthSummary.StdTotalYDistanceTraveled(nc_idx);
                    this.EmbryoStats.StdTotalDistanceTraveled(SetIndex, NC-8) = HealthSummary.StdTotalDistanceTraveled(nc_idx);
                    this.EmbryoStats.StdDistanceTraveledPerSecond(SetIndex, NC-8) =HealthSummary.StdDistanceTraveledPerSecond(nc_idx);
                    this.EmbryoStats.StdTotalXDisplacement(SetIndex, NC-8) = HealthSummary.StdTotalXDisplacement(nc_idx);
                    this.EmbryoStats.StdTotalYDisplacement(SetIndex, NC-8) = HealthSummary.StdTotalYDisplacement(nc_idx);
                    this.EmbryoStats.StdTotalDisplacement(SetIndex, NC-8) = HealthSummary.StdTotalDisplacement(nc_idx);
                    this.EmbryoStats.StdDisplacementPerSecond(SetIndex, NC-8) = HealthSummary.StdDisplacementPerSecond(nc_idx);
                    if NC < 14
                        this.EmbryoStats.NCDivisionInfo(SetIndex, NC-8) = HealthSummary.NCDivisionInfo(nc_idx);
                        this.EmbryoStats.DivisionStdInfo(SetIndex, NC-8) = HealthSummary.DivisionStdInfo(nc_idx);
                    end
                end
            end
  
        end
        
        
    end
end