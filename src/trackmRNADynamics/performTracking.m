function [Particles, SpotFilter] = performTracking(Particles, scurrentChannelitzcells, NCh, Spots, app, SpotFilter, PreProcPath, Prefix, UseHistone, ParticlesFig, SpotsChannel, NDigits, NucleiFig, particlesAxes, nucAxes, Ellipses, PixelSize, SearchRadius, ExperimentType, FrameInfo, retrack, displayFigures)
  % Iterate over all channels
  for Channel = 1:NCh

    % Iterate over all frames
    for CurrentFrame = 1:length(Spots{Channel})

      if isempty(app) && displayFigures
        figure(ParticlesFig)
        set(ParticlesFig, 'units', 'normalized', 'position', [0.01, .55, .33, .33]);          
      end

      % Get the filter for this frame
      CurrentFrameFilter = logical(SpotFilter{Channel}(CurrentFrame, ~isnan(SpotFilter{Channel}(CurrentFrame, :))));

      xPos = displayParticlesFigure(app, particlesAxes, ParticlesFig, Spots, Channel, CurrentFrame, ...
        CurrentFrameFilter, PreProcPath, Prefix, SpotsChannel, FrameInfo, displayFigures);

      if UseHistone
        [Particles, SpotFilter] = trackParticlesBasedOnNuclei(PreProcPath, Prefix, CurrentFrame, NDigits, app, nucAxes, Ellipses, ...
          ExperimentType, Channel, scurrentChannelitzcells, Particles, Spots, SpotFilter, PixelSize, SearchRadius, retrack, displayFigures);
      else
        [Particles] = trackParticlesBasedOnProximity(Particles, Spots, xPos, SpotFilter, Channel, CurrentFrame, PixelSize, SearchRadius, retrack, displayFigures);
      end

    end

  end

  if isempty(app) && displayFigures
    close(ParticlesFig)
    if UseHistone
      close(NucleiFig)
    end
  end

  
  for currentChannel = 1:NCh
    
    if ~isfield(Particles{currentChannel}, 'FrameApproved')
      
      for i = 1:length(Particles{currentChannel})
        Particles{currentChannel}(i).FrameApproved = true(size(Particles{currentChannel}(i).Frame));
      end
      
    else
      
      for i = 1:length(Particles{currentChannel})
        
        if isempty(Particles{currentChannel}(i).FrameApproved)
          Particles{currentChannel}(i).FrameApproved = true(size(Particles{currentChannel}(i).Frame));
        end
        
      end
      
    end
    
    %AR: add x, y and z positions to the Particles structure. This was
    %originally only added by addParticlePosition, but it it's more useful
    %early on in the pipeline. 
    
        for i=1:length(Particles{currentChannel})
            for j=1:length(Particles{currentChannel}(i).Frame)
                [x,y,z]=SpotsXYZ(Spots{currentChannel}(Particles{currentChannel}(i).Frame(j)));
                if ~isempty(x) 
                    Particles{currentChannel}(i).xPos(j)=x(Particles{currentChannel}(i).Index(j));
                    Particles{currentChannel}(i).yPos(j)=y(Particles{currentChannel}(i).Index(j));
                    Particles{currentChannel}(i).zPos(j)=z(Particles{currentChannel}(i).Index(j));
                end
            end
        end
    
  end
  
  % If we only have one channel, then convert SpotFilter and Particles to a standard structure.
  if NCh == 1
    SpotFilter = SpotFilter{1};
    Particles = Particles{1};
  end

end

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

function [Particles, SpotFilter] = trackParticlesBasedOnNuclei(PreProcPath, Prefix, CurrentFrame, NDigits, app, nucAxes, Ellipses, ...
    ExperimentType, Channel, scurrentChannelitzcells, Particles, Spots, SpotFilter, PixelSize, SearchRadius, retrack, displayFigures)

if displayFigures
  hisImage = openHistoneImage(Prefix, PreProcPath, CurrentFrame, NDigits);

  if ~isempty(app)
    ax2 = app{2};
  else
    ax2 = nucAxes;
  end

  imshow(hisImage, [], 'Border', 'Tight', 'Parent', ax2, 'InitialMagnification', 'fit')
  hold(ax2, 'on')
  PlotHandle = [];
  [NEllipses, ~] = size(Ellipses{CurrentFrame});

  for EllipsesIndex = 1:NEllipses
    PlotHandle = [PlotHandle, ellipse(...
      Ellipses{CurrentFrame}(EllipsesIndex, 3), ...
      Ellipses{CurrentFrame}(EllipsesIndex, 4), ...
      Ellipses{CurrentFrame}(EllipsesIndex, 5), ...
      Ellipses{CurrentFrame}(EllipsesIndex, 1) + 1, ...
      Ellipses{CurrentFrame}(EllipsesIndex, 2) + 1, ...
      [], [], ax2)];

    text(ax2, Ellipses{CurrentFrame}(EllipsesIndex, 1) + 1, Ellipses{CurrentFrame}(EllipsesIndex, 2) + 1, ...
      num2str(EllipsesIndex), 'BackgroundColor', [.7 .9 .7]);
  end

  set(PlotHandle, 'Color', 'r')
  hold(ax2, 'off')
  title(ax2, CurrentFrame)
  drawnow
end

  if strcmp(ExperimentType, '2spot')
    SpotsPerNucleus = 2;
  else
     SpotsPerNucleus = 1;
  end

  [Particles{Channel}, SpotFilter{Channel}] = AssignParticle2Nucleus(scurrentChannelitzcells, Ellipses, ...
    Particles{Channel}, Spots{Channel}, SpotFilter{Channel}, CurrentFrame, PixelSize, SpotsPerNucleus, retrack);

end

function hisImage = openHistoneImage(Prefix, PreProcPath, CurrentFrame, NDigits)
  
  HistoneImageFileNamePrefix = [PreProcPath, filesep, Prefix, filesep, Prefix];
  HistoneImageFileNameSuffix = [iIndex(CurrentFrame, NDigits), '.tif'];

  try
    hisImage = imread([HistoneImageFileNamePrefix, '-His_', HistoneImageFileNameSuffix]);
  catch

    try
      hisImage = imread([HistoneImageFileNamePrefix, '_His_', HistoneImageFileNameSuffix]);
    catch
      hisImage = 0;
    end

  end

end
