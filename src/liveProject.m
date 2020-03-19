classdef liveProject
    %liveProject Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        Project = '';
        includedExperimentNames = {};
        includedExperiments = {};

        ignoredExperimentNames = {};
        
        dataStatus = {};
        
        
    end
    
    properties (Hidden = true)
        
        ignoredExperiments = {};

    end
    
    methods
        %% Constructors
        
        
        function obj = liveProject(Project)
            %liveProject Construct an instance of this class
            %   Detailed explanation goes here
            obj.Project = Project;
            
            [~, experiment, ~, ignoredExperiment, obj.dataStatus] =...
                LoadMS2Sets(Project, 'justPrefixes', 'noCompiledNuclei');
            
            obj.includedExperimentNames = experiment;
            obj.ignoredExperimentNames = ignoredExperiment;

            
            for i = 1:length(experiment)
                obj.includedExperiments{i} = liveExperiment(experiment{i});
            end
            
            for j = 1:length(ignoredExperiment)
                obj.ignoredExperiments{j} = liveExperiment(ignoredExperiment{j});
            end
            
            
            
            
        end
        
        
        
        
        %% Methods
                
        
    end
end

