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