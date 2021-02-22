function [Particles,falsePositives, Spots2] = findBrightestZ(Particles, num_shadows,...
    use_integral_center, force_z, Spots, varargin)
% Particles = findBrightestZ(Particles)
%
% DESCRIPTION
% Sub-function for segmentation that tracks transcription loci along the
% z-axis. Takes in a Spots structure and outputs a Spots structure.
%
% ARGUMENTS
% 'Particles':  A structure array with a list of detected transcriptional loci
% in each frame and their properties.
%
% OPTIONS
%
% OUTPUT
% 'Particles':  A structure array with a list of detected transcriptional loci
% in each frame and their properties.
%
% Author (contact): Nick Lammers (nlammers@berkeley.edu) & Armando Reimer (areimer@berkeley.edu)
% Created: 01/01/2016
% Last Updated: 7/19/2018
%
% Documented by: Armando Reimer (areimer@berkeley.edu)



numFrames = length(Spots);
Spots2 = repmat(struct('Fits', []), 1, numFrames);

falsePositives = 0;
for i = 1:length(Particles)
    
    z_vec = [Particles(i).z]; %convenience vector
    %pull intensity value from particle snippets
    RawIntensityVec = [Particles(i).FixedAreaIntensity];
    CentralIntensityVec = [Particles(i).CentralIntensity];
    %find slice with brightest pixel
    [~, MaxIndexCentral] = max(CentralIntensityVec);
    [~, MaxIndexRaw] = max(RawIntensityVec);
    % calculate convenience vectors
    z_grid = min(z_vec):max(z_vec);
    z_raw_values = zeros(size(z_grid));
    z_raw_values(ismember(z_grid,z_vec)) = RawIntensityVec;
    if ~use_integral_center
        CentralZ = z_vec(MaxIndexCentral);
        ZStackIndex = MaxIndexCentral;
    else
        if length(z_vec) < 3 %treat thinner spots separately
            CentralZ = z_vec(MaxIndexRaw);
            ZStackIndex = MaxIndexRaw;
        else
            % Convolve with gaussian filter to find "best" center
            g = [-1 0 1];
            gaussFilter = exp(-g .^ 2 / (2));
            RawRefVec = conv(gaussFilter,z_raw_values);
            RawRefVec = RawRefVec(2:end-1);
            RawRefVec(1) = NaN;
            RawRefVec(end) = NaN;
            RawRefVec = RawRefVec(ismember(z_grid,z_vec));
            [~, MaxIndexIntegral] = nanmax(RawRefVec);
            CentralZ = z_vec(MaxIndexIntegral);
            ZStackIndex = MaxIndexIntegral;
        end
    end
    
    %allow the function call to choose the "brightest" z plane rather
    %than automatically determining it
    if ~force_z
        Particles(i).brightestZ = CentralZ;
    else
        Particles(i).brightestZ = force_z;
    end
    
    
    Particles(i).FixedAreaIntensity3 = sum(z_raw_values(ismember(z_grid,Particles(i).brightestZ-1:Particles(i).brightestZ+1)));
    
    
    %use convolution kernel to look for shadows
    z_raw_binary = ~isnan(z_raw_values);
    z_shadow_vec = conv(z_raw_binary,[1 1 1],'same');
    z_shadow_vec = z_shadow_vec(ismember(z_grid,z_vec));
    n_shadows = z_shadow_vec(ZStackIndex)-1;
    
    if n_shadows < num_shadows
        Particles(i).discardThis = 1;
        falsePositives = falsePositives + 1;
    end
end



if isstruct(Particles)
    Particles = rmfield(Particles, 'r');
    Particles = rmfield(Particles, 'discardThis');
end

if ~isempty(Particles)
    Particles.snippet_size = Particles.snippet_size(1);
    Particles.intArea = Particles.intArea(1);
end

if ~isempty(Spots)
    falsePositives = 0;
    for frame = 1:length(Spots)
        nSpots = length(Spots(frame).Fits);
        for spotIndex = 1:nSpots
            z_vec = [Spots(frame).Fits(spotIndex).z]; %convenience vector
            %pull intensity value from particle snippets
            RawIntensityVec = [Spots(frame).Fits(spotIndex).FixedAreaIntensity];
            CentralIntensityVec = [Spots(frame).Fits(spotIndex).CentralIntensity];
            %find slice with brightest pixel
            [~, MaxIndexCentral] = max(CentralIntensityVec);
            [~, MaxIndexRaw] = max(RawIntensityVec);
            % calculate convenience vectors
            z_grid = min(z_vec):max(z_vec);
            z_raw_values = zeros(size(z_grid));
            z_raw_values(ismember(z_grid,z_vec)) = RawIntensityVec;
            z_raw_values(z_raw_values==0) = NaN;
            if ~use_integral_center
                CentralZ = z_vec(MaxIndexCentral);
                ZStackIndex = MaxIndexCentral;
            else
                if length(z_vec) < 3 %treat thinner spots separately
                    CentralZ = z_vec(MaxIndexRaw);
                    ZStackIndex = MaxIndexRaw;
                else
                    % Convolve with gaussian filter to find "best" center
                    g = [-1 0 1];
                    gaussFilter = exp(-g .^ 2 / (2));
                    RawRefVec = conv(gaussFilter,z_raw_values);
                    RawRefVec = RawRefVec(2:end-1);
                    RawRefVec(1) = NaN;
                    RawRefVec(end) = NaN;
                    RawRefVec = RawRefVec(ismember(z_grid,z_vec));
                    [~, MaxIndexIntegral] = nanmax(RawRefVec);
                    CentralZ = z_vec(MaxIndexIntegral);
                    ZStackIndex = MaxIndexIntegral;
                end
            end
            
            %allow the function call to choose the "brightest" z plane rather
            %than automatically determining it
            if ~force_z
                Spots(frame).Fits(spotIndex).brightestZ = uint8(CentralZ);
            else
                Spots(frame).Fits(spotIndex).brightestZ = uint8(force_z);
            end
            Spots(frame).Fits(spotIndex).FixedAreaIntensity3 = single(sum(z_raw_values(ismember(z_grid,Spots(frame).Fits(spotIndex).brightestZ-1:Spots(frame).Fits(spotIndex).brightestZ+1))));
            
            
            %use convolution kernel to look for shadows
            z_raw_binary = ~isnan(z_raw_values);
            z_shadow_vec = conv(z_raw_binary,[1 1 1],'same');
            z_shadow_vec = z_shadow_vec(ismember(z_grid,z_vec));
            n_shadows = z_shadow_vec(ZStackIndex)-1;
            
            if n_shadows < num_shadows
                Spots(frame).Fits(spotIndex).discardThis = true;
                falsePositives = falsePositives + 1;
            else
                Spots2(frame).Fits = [Spots2(frame).Fits, Spots(frame).Fits(spotIndex)];
            end
        end
    end
    
    for i = 1:length(Spots2)
        if isstruct(Spots2(i).Fits)
            Spots2(i).Fits = rmfield(Spots2(i).Fits, 'r');
            Spots2(i).Fits = rmfield(Spots2(i).Fits, 'discardThis');
        end
    end
    
end

end