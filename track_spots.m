function Particles = track_spots(Particles, neighb)

% time tracking

fields = fieldnames(Particles);
changes = 1;
while changes ~= 0
    changes = 0;
    for n = 1:length(Particles)-1 %particle of interest
        for j = n+1:length(Particles) %particle to compare to
            if Particles(n).t(end) == (Particles(j).t(end) -  1)
                dist = sqrt( (Particles(n).x(end) - Particles(j).x(end))^2 + (Particles(n).y(end) - Particles(j).y(end))^2); 
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