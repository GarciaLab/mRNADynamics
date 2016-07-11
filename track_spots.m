function Particles = track_spots(Particles, neighb)

% time tracking

    fields = fieldnames(Particles);
    changes = 1;
    while changes ~= 0
        changes = 0;
        for n = 1:length(Particles)-1 %particle of interest
            for j = n+1:length(Particles) %particle to compare to
                if Particles(n).frame(end) == (Particles(j).frame(end) -  1)
                    dist = sqrt( (Particles(n).xFit(end) - Particles(j).xFit(end))^2 + (Particles(n).yFit(end) - Particles(j).yFit(end))^2); 
                    if dist < neighb
                        for m = 1:numel(fields)-1 %do not include fields 'r'
                            Particles(n).(fields{m}) = [Particles(n).(fields{m}), Particles(j).(fields{m})];
                        end
                        Particles(j).r = 1;
                        changes = changes + 1;
                    end
                end
            end
        end
    Particles = Particles([Particles.r]~=1);
end

    for i = 1:length(Particles)
        if length(Particles(i).frame) > 70
            plot(Particles(i).frame, Particles(i).GaussianIntensity)
            i
            hold on
        end
    end
    MeanVectorAll = NaN(1, num_frames);
    for i = 1:length(Spots)
        for j = 1:length(Spots(i).Fits)
            MeanVectorAll(i) = 0;
            MeanVectorAll(i) = MeanVectorAll(i) + Spots(i).Fits(j).GaussianIntensity(Spots(i).Fits(j).brightestZ);
        end
        MeanVectorAll(i) = MeanVectorAll(i) / length(Spots(i).Fits);
    end
end