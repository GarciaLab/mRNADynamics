classdef liveProject
    %liveProject Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        Project = '';
        includedExperimentNames = {};
        includedExperiments = {};
        
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
            
            [~, experiment, ~] = LoadMS2Sets(Project, 'justPrefixes', 'noCompiledNuclei');
            
            obj.includedExperimentNames = experiment;
            
            for i = 1:length(experiment)
                obj.includedExperiments{i} = liveExperiment(experiment{i});
            end
            
            obj.ignoredExperiments = {};
            
        end
        
        
        
        
        %% Methods
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

