function x = displayParticlesFigure(app, particlesAxes, ParticlesFig, Spots, Channel, CurrentFrame, CurrentFrameFilter, PreProcPath, Prefix, SpotsChannel, FrameInfo, displayFigures)

% Get the positions of the spots in this frame
[x, y, ~] = SpotsXYZ(Spots{Channel}(CurrentFrame));

if displayFigures
    % Z plane to be displayed. We use the median of all particles found
    % in this frame
    if ~isempty(Spots{Channel}(CurrentFrame).Fits)
        CurrentZ = round(median([Spots{Channel}(CurrentFrame).Fits.brightestZ]));
    else
        CurrentZ = round(FrameInfo(1).NumberSlices / 2);
    end
    
    % Load the corresponding mRNA image. Check whether we have multiple
    % channels saved or not.
    FileNamePrefix = [PreProcPath, filesep, Prefix, filesep, Prefix, '_', iIndex(CurrentFrame, 3), '_z', iIndex(CurrentZ, 2)];
    particleImage = imread([FileNamePrefix, '_ch', iIndex(SpotsChannel(Channel), 2), '.tif']);
    
    FigureName = ['Ch', num2str(Channel), '  Frame: ', num2str(CurrentFrame), '/', num2str(length(Spots{Channel}))];
    
    if ~isempty(app)
        ax1 = app{1};
        title(ax1, FigureName)
    else
        ax1 = particlesAxes;
        set(ParticlesFig, 'Name', FigureName);
        title(ax1, CurrentFrame)
    end
    
    imshow(particleImage, [], 'Parent', ax1, 'InitialMagnification', 'fit')
    hold(ax1, 'on')
    plot(ax1, x(CurrentFrameFilter), y(CurrentFrameFilter), 'or', 'MarkerSize', 10)
    plot(ax1, x(~CurrentFrameFilter), y(~CurrentFrameFilter), 'ow', 'MarkerSize', 10)
    hold(ax1, 'off')
end

end