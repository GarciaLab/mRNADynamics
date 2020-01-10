 function inputHandler = FitApproveEventHandler(cptState)

    function fit_approve(~, ~)
        % Define the fitApproved as true
        cptState.fitApproved = 1; %JP: is this used by anyone?
        % save the fitted values (Coefficients and lineFitHandle) in Particles.mat
        % For now, I will save lineFitHandle (which is a plot for the initial
        % slope), but it might take up too much space, then I need a better
        % way to regenerate the plot, using the Coefficients and fittedFrames
        % Quality control : Check whether the slope is positive
        if cptState.Coefficients(1, 1) > 0
            cptState.Particles{cptState.CurrentChannel}(cptState.CurrentParticle).fitApproved = 1;
            cptState.Particles{cptState.CurrentChannel}(cptState.CurrentParticle).Coefficients = cptState.Coefficients;
            cptState.Particles{cptState.CurrentChannel}(cptState.CurrentParticle).fittedFrames = cptState.FrameIndicesToFit; % use the index of particle trace for convenience
        else
            cptState.Particles{cptState.CurrentChannel}(cptState.CurrentParticle).fitApproved = 0;
            cptState.Particles{cptState.CurrentChannel}(cptState.CurrentParticle).Coefficients = [];
            cptState.Particles{cptState.CurrentChannel}(cptState.CurrentParticle).fittedFrames = [];
        end
        
    end


    keyInputHandler = @fit_approve;
end


