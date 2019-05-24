function analyzeContiguity(movie)

    %movie is the data set we want to look into. 
    
    compiledParticles = movie.CompiledParticles;
    nParticles = length(compiledParticles);
    contiguityLong = [];
    contiguity2Long = [];
    contiguity = zeros(1, nParticles);
    contiguity2 = zeros(1, nParticles);
    channel = 1; %no support for 2 channels at the moment
    
    for i = 1:nParticles
        if iscell(compiledParticles)
            compiledParticles = compiledParticles{channel};
        end
        frames = compiledParticles(i).Frame; len = length(frames); contiguity(i) = len; contiguity2(i) = len; missed = 0; missed2 = 0;
       
        if len > 1
           for j = 2:len
               frameInterval = frames(j) - frames(j-1);
               missed = missed + frameInterval - (len + 1);
               if frameInterval > 1 
                    missed2 = missed2 + 1; %counting gaps. 
               end
           end
       end
       
       contiguity(i) = (contiguity(i) - missed)/len; %every gap weighted by its duration
       contiguity2(i) = (contiguity2(i) - missed2)/len; %this counts the gaps and normalizes by full trace length

       if len > 1
           contiguityLong = [contiguityLong,contiguity(i)]; %every gap weighted by its duration but traces one frame long excluded
           contiguity2Long = [contiguity2Long,contiguity2(i)]; %analagous to above
       end
    end

    figure('Name', 'contiguity')
    subplot(2, 2, 1)
    h = histogram(contiguity);
    title({'contiguity of traces relative to';' trace length weighted by'; 'length of gaps'});
    xlabel('contiguity metric')
    ylabel('counts')
    standardizeFigure(gca, [], 'red')
    subplot(2, 2, 2)
    h = histogram(contiguityLong);
    title({'contiguity of traces > 1 frame';'relative to trace length';'weighted by lengths of gaps'});
    xlabel('contiguity metric')
    ylabel('counts')
    standardizeFigure(gca, [], 'red')
    subplot(2, 2, 3)
    h = histogram(contiguity2);
    title({'contiguity of traces';'relative to trace length'});
    xlabel('contiguity metric')
    ylabel('counts')
    standardizeFigure(gca, [], 'red')
    subplot(2, 2, 4)
    h = histogram(contiguity2Long);
    title({'contiguity of traces > 1 frame';' relative to trace length'});
    xlabel('contiguity metric')
    ylabel('counts')
    standardizeFigure(gca, [], 'red')
    
end