classdef PlotTraceSettings < handle
    properties
        AmpIntegral
        AmpIntegral3
        
        ErrorIntegral
        ErrorIntegral3
        
        backGround3
        
        GaussIntegral
        AmpIntegralGauss3D
        ErrorIntegralGauss3D
        
        Spots3D
        
        UseCompiledParticles
        
        UseTwinTraces
        
        PlotInputChannel
    end
    
    methods
        function this = PlotTraceSettings(UseTwinTraces, PlotInputChannel, UseCompiledParticles)
            if exist('UseTwinTraces', 'var')
                this.UseTwinTraces = UseTwinTraces;
            else
                this.UseTwinTraces = false;
            end
            
            if exist('PlotInputChannel', 'var')
                this.PlotInputChannel = PlotInputChannel;
            else
                this.PlotInputChannel = false;
            end
            
            if exist('UseCompiledParticles', 'var')
                this.UseCompiledParticles = UseCompiledParticles;
            else
                this.UseCompiledParticles = false;
            end
        end
    end

end
