classdef LiveProject
    %liveProject Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        Project = '';
        includedExperimentNames = [];
        includedExperiments = {};
        
        unhealthyNames = [];
        
        ignoredExperimentNames = [];
        
        dataStatus = {};
        
        hasSpots = [];
        hasParticles = [];
        hasSchnitzcells = [];
        hasCompiledParticles = [];
        
        anaphaseFramesAnnotated = [];
        
        
    end
    
    properties (Hidden = true)
        
        ignoredExperiments = {};
        
    end
    
    methods
        %% Constructors
        
        
        function this = liveProject(Project)
            %liveProject Construct an instance of this class
            %   Detailed explanation goes here
            this.Project = Project;
            
            [~, experiment, ~, ignoredExperiment, this.dataStatus] =...
                LoadMS2Sets(Project, 'justPrefixes', 'noCompiledNuclei');
            
            this.dataStatus = string(this.dataStatus);
            this.includedExperimentNames = string(experiment);
            this.ignoredExperimentNames = string(ignoredExperiment);
            
            for i = 1:length(experiment)
                this.includedExperiments{i} = LiveExperiment(experiment{i});
                isUnhealthy = this.includedExperiments{i}.isUnhealthy;
                if ~isnan(isUnhealthy) && isUnhealthy
                    this.unhealthyNames = [ this.unhealthyNames,...
                        string(experiment{i})];
                end
            end
            
            for j = 1:length(ignoredExperiment)
                this.ignoredExperiments{j} = LiveExperiment(ignoredExperiment{j});
            end
            
            this.hasSpots = haveSpots(this);
            this.hasParticles = haveParticles(this);
            this.hasSchnitzcells =  haveSchnitzcells(this);
            this.hasCompiledParticles = haveCompiledParticles(this);
            this.anaphaseFramesAnnotated = haveAnaphaseFrames(this);
            
            
        end
        
        
        
        
        %% Methods
        
        function hasSchnitzcells = haveSchnitzcells(this)
            
            numValidProjects = length(this.includedExperimentNames);
            
            hasSchnitzcells = zeros(1, numValidProjects);
            
            for k = 1:numValidProjects
                hasSchnitzcells(k) = ...
                    logical(this.includedExperiments{k}.hasSchnitzcellsFile);
            end
            
        end
        
        function hasParticles = haveParticles(this)
            
            numValidProjects = length(this.includedExperimentNames);
            
            hasParticles = zeros(1, numValidProjects);
            
            for k = 1:numValidProjects
                hasParticles(k) = ...
                    logical(this.includedExperiments{k}.hasParticlesFile);
            end
            
        end
        
        function hasSpots = haveSpots(this)
            
            numValidProjects = length(this.includedExperimentNames);
            
            hasSpots = zeros(1, numValidProjects);
            
            for k = 1:numValidProjects
                hasSpots(k) = ...
                    logical(this.includedExperiments{k}.hasSpotsFile);
            end
            
        end
        
        function hasCompiledParticles = haveCompiledParticles(this)
            
            numValidProjects = length(this.includedExperimentNames);
            
            hasCompiledParticles = zeros(1, numValidProjects);
            
            for k = 1:numValidProjects
                hasCompiledParticles(k) = ...
                    logical(this.includedExperiments{k}.hasCompiledParticlesFile);
            end
            
        end
        
        function anaphaseFramesAnnotated = haveAnaphaseFrames(this)
            
            numValidProjects = length(this.includedExperimentNames);
            
            anaphaseFramesAnnotated = zeros(1, numValidProjects);
            
            for k = 1:numValidProjects
                anaphaseFramesAnnotated(k) = ...
                    nansum(this.includedExperiments{k}.anaphaseFrames) > 1;
            end
            
        end
        
        
        
    end
end

