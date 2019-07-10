function [Overlay, overlayAxes, snippetFigAxes, rawDataAxes, gaussianAxes, traceFigAxes, zProfileFigAxes,...
  zTraceAxes, HisOverlayFig,HisOverlayFigAxes] = checkParticleTracking_drawGUI(UseHistoneOverlay)
  Overlay = figure;
  
  if UseHistoneOverlay
    HisOverlayFig = figure;
    HisOverlayFigAxes = axes(HisOverlayFig);
  else
    HisOverlayFig = [];
    HisOverlayFigAxes = [];
  end


  overlayAxes = subplot(1, 2, 1, 'Parent', Overlay);
  tb = axtoolbar(overlayAxes);
  traceFigAxes = subplot(1, 2, 2, 'Parent', Overlay);

  zFig = figure;
  zProfileFigAxes = subplot(1, 2, 1, 'Parent', zFig);
  zTraceAxes = subplot(1, 2, 2, 'Parent', zFig);

  snipFig = figure();
  snippetFigAxes = subplot(1, 3, 1, 'Parent', snipFig);
  rawDataAxes = subplot(1, 3, 2, 'Parent', snipFig);
  gaussianAxes = subplot(1, 3, 3, 'Parent', snipFig);
  
  set(Overlay, 'units', 'normalized', 'position', [0.01, .45, .82, .33]);

  if UseHistoneOverlay
    set(HisOverlayFig, 'units', 'normalized', 'position', [0.01, 0.1, .33, .33]);
  end

  set(overlayAxes, 'units', 'normalized', 'position', [-.25 .06 .9 .9])
  set(traceFigAxes, 'units', 'normalized', 'position', [.48 .17 .48 .63])
  set(snipFig, 'units', 'normalized', 'position', [0.355, 0.15, 3 * (.2 / 2), .33 / 2]);
  set(zFig, 'units', 'normalized', 'position', [0.67, 0.15, .2, .33 / 2]);

end