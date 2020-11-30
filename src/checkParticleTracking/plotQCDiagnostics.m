function qcAxes = plotQCDiagnostics(qcAxes, cptState)

    if ~isempty(cptState)
        particle = cptState.Particles{cptState.CurrentChannelIndex}(cptState.CurrentParticle);

        inputScores = particle.qcScoreArray(particle.Frame==cptState.CurrentFrame,:);
        if isempty(inputScores)
          inputScores = zeros(1,size(particle.qcScoreArray,2));
        end

        nInputs = length(inputScores);
        inputIndex = 1:nInputs;

        labelCell = cptState.trackingOptions.qcFieldNames;

        inputScores(end+1) = inputScores(1);
        logLthresh = particle.qcThreshVec;
        logLthresh(end+1) = logLthresh(1);
        inputScores = inputScores./logLthresh;

        % calculate angles
        inputAngleIncrement = 2*pi/nInputs;
        inputAngles = (inputIndex-1)*inputAngleIncrement;
        inputAngles(end+1) = inputAngles(1);

        set(0,'CurrentFigure', qcAxes)
        title('Parameter QC breakdown')
        polarplot(inputAngles, ones(size(logLthresh)))
        hold on
        polarplot(inputAngles, inputScores)

        pax = gca;
        pax.ThetaAxisUnits = 'radians';
        thetaticks(inputAngles(1:end-1))
        thetaticklabels(labelCell)

        hold off
    else
        nInputs = 5;
        inputIndex = 1:nInputs;
        inputAngleIncrement = 2*pi/nInputs;
        inputAngles = (inputIndex-1)*inputAngleIncrement;
        inputAngles(end+1) = inputAngles(1);
        
        set(0,'CurrentFigure', qcAxes)
        title('Parameter QC breakdown')
        polarplot(inputAngles, ones(size(inputAngles)))
    end