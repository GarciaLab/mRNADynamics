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
    end
    
    methods
        function this = PlotTraceSettings(UseCompiledParticles)
            if exist('UseCompiledParticles', 'var')
                this.UseCompiledParticles = UseCompiledParticles;
            else
                this.UseCompiledParticles = false;
            end
        end
    end

end
