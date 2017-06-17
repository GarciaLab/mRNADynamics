function Particles = track_spots(Particles, Spots, neighb, num_frames)
% track_spots(Particles, Spots, neighb, num_frames)
%
% DESCRIPTION
% Rough and ready tracking function that can be run directly after
% segmentSpots or a derivative. The output can be used for rudimentary
% analysis of the sort done during/after CompileParticles.
%
% ARGUMENTS
% Particles: Transcriptional loci that have been z-tracked.
% Spots: Spots is a restructured Particles made by segmentSpots. Mostly
% here for plotting purposes.
% neighb: A distance threshold for deciding how far a spot should move in
% one frame.
% num_frames: Here for plotting purposes. The number of frames in a movie. 
%
% OUTPUT
% Particles:  A structure array with a list of detected transcriptional
% loci. This will be restrucured by segmentSpots (or a derivative) into a
% list of particles per frame (the Spots structure).
%
% Author (contact): Armando Reimer (areimer@berkeley.edu)
% Created: 01/01/2016
% Last Updated: 12/31/2016
%
% Documented by: Armando Reimer (areimer@berkeley.edu)           

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
            hold on
        end
    end
    
    %Some rudimentary analysis and plotting abilities. Commented out to
    %to make this function more singular purpose.
    
%     MeanVectorAll = NaN(1, num_frames);
%     for i = 1:length(Spots)
%         for j = 1:length(Spots(i).Fits)
%             MeanVectorAll(i) = 0;
%             MeanVectorAll(i) = MeanVectorAll(i) + Spots(i).Fits(j).GaussianIntensity(Spots(i).Fits(j).brightestZ);
%         end
%         MeanVectorAll(i) = MeanVectorAll(i) / length(Spots(i).Fits);
%     end
end