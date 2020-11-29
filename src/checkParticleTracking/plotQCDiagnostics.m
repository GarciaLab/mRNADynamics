function qcAxes = plotQCDiagnostics(qcAxes)

    inputScores = rand(1,6).*10.*rand(1,6);
    nInputs = length(inputScores);
    inputIndex = 1:nInputs;

    labelCell = {'parameter1','parameter2','parameter3','parameter4','parameter5','parameter6'};

    inputScores(end+1) = inputScores(1);
    logLthresh = ones(size(inputScores));

    % calculate angles
    inputAngleIncrement = 2*pi/nInputs;
    inputAngles = (inputIndex-1)*inputAngleIncrement;
    inputAngles(end+1) = inputAngles(1);

    set(0,'CurrentFigure', qcAxes)

    polarplot(inputAngles, logLthresh)
    hold on
    polarplot(inputAngles, inputScores)

    pax = gca;
    pax.ThetaAxisUnits = 'radians';
    thetaticks(inputAngles(1:end-1))
    thetaticklabels(labelCell)