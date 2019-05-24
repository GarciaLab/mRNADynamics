function plotWindowTimings(movie)
 
    channel = 1; %no support for 2 channel
    %movie is the data set we want to look into. 
    compiledParticles = movie.CompiledParticles;
    if iscell(compiledParticles)
        compiledParticles = compiledParticles{channel};
    end
    nParticles = length(compiledParticles);
    onTimes = [];
    offTimes = [];
    durations = [];
    
    for i = 1:nParticles
        frames = compiledParticles(i).Frame; 
        onTimes(i) = frames(1);
        offTimes(i) = frames(end);
        duration(i) = frames(end) - frames(1);
    end
    
    figure('Name', 'timings')
    subplot(1, 3, 1)
    h = histogram(onTimes);
    title('on times')
    xlabel('on time (frames)')
    ylabel('counts')
    standardizeFigure(gca, [], 'red')
    subplot(1, 3, 2)
    h = histogram(offTimes);
    title('off times')
    xlabel('off time (frames)')
    ylabel('counts')
    standardizeFigure(gca, [], 'red')
    subplot(1, 3, 3)
    h = histogram(duration);
    title('duration')
    xlabel('duration of transcription (frames)')
    ylabel('counts')
    standardizeFigure(gca, [], 'red')


end