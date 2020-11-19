function imStack = extractFrameStack(frame,ch,movieMatCh)

    if~isempty(mocieMatCh)
      imStack = movieMatCh(:, :, :,frame); 
    else
      imStack = getMovieFrame(liveExperiment, frame, ch);
    end