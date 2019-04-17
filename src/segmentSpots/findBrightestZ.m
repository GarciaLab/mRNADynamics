 function [Particles,falsePositives, Spots2] = findBrightestZ(Particles, num_shadows, use_integral_center, force_z, Spots)
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
        
         %AR 7/19/2018- I'm appropriating these variables in order to
        %integrate within an ellipsoid volume. Not sure of the original
        %purpose but it seems deprecated. 
        
         % if there are insufficient slices, these metrics will register as NaNs
       % RawIntegral3 = mean(z_raw_values(ismember(z_grid,CentralZ-1:CentralZ+1)));
       % RawIntegral5 = mean(z_raw_values(ismember(z_grid,CentralZ-2:CentralZ+2)));
        %Particles(i).FixedAreaIntensity3 = RawIntegral3;
        %Particles(i).FixedAreaIntensity5 = RawIntegral5;
%             Particles(i).FixedAreaIntensity3 = Particles(i).FixedAreaIntensity(Particles(i).brightestZ - 1) + Particles(i).FixedAreaIntensity(Particles(i).brightestZ) + Particles(i).FixedAreaIntensity(Particles(i).brightestZ + 1);
%             Particles(i).FixedAreaIntensity5 = Particles(i).FixedAreaIntensity(Particles(i).brightestZ - 2) + Particles(i).FixedAreaIntensity(Particles(i).brightestZ - 1) + Particles(i).FixedAreaIntensity(Particles(i).brightestZ) + Particles(i).FixedAreaIntensity(Particles(i).brightestZ + 1) + Particles(i).FixedAreaIntensity(Particles(i).brightestZ + 2);
        Particles(i).FixedAreaIntensity3 = sum(z_raw_values(ismember(z_grid,Particles(i).brightestZ-1:Particles(i).brightestZ+1)));
        Particles(i).FixedAreaIntensity5 = sum(z_raw_values(ismember(z_grid,Particles(i).brightestZ-2:Particles(i).brightestZ+2)));
        try
            Particles(i).cylIntensity = Particles(i).cylIntensity(ZStackIndex);
        catch                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
        end
    
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
            Spots(frame).Fits(spotIndex).brightestZ = CentralZ;
        else
            Spots(frame).Fits(spotIndex).brightestZ = force_z;
        end
        
        %AR 7/19/2018- I'm appropriating these variables in order to
        %integrate within an ellipsoid volume. Not sure of the original
        %purpose but it seems deprecated.
        
        % if there are insufficient slices, these metrics will register as NaNs
        % RawIntegral3 = mean(z_raw_values(ismember(z_grid,CentralZ-1:CentralZ+1)));
        % RawIntegral5 = mean(z_raw_values(ismember(z_grid,CentralZ-2:CentralZ+2)));
        %Particles(i).FixedAreaIntensity3 = RawIntegral3;
        %Particles(i).FixedAreaIntensity5 = RawIntegral5;
        %             Particles(i).FixedAreaIntensity3 = Particles(i).FixedAreaIntensity(Particles(i).brightestZ - 1) + Particles(i).FixedAreaIntensity(Particles(i).brightestZ) + Particles(i).FixedAreaIntensity(Particles(i).brightestZ + 1);
        %             Particles(i).FixedAreaIntensity5 = Particles(i).FixedAreaIntensity(Particles(i).brightestZ - 2) + Particles(i).FixedAreaIntensity(Particles(i).brightestZ - 1) + Particles(i).FixedAreaIntensity(Particles(i).brightestZ) + Particles(i).FixedAreaIntensity(Particles(i).brightestZ + 1) + Particles(i).FixedAreaIntensity(Particles(i).brightestZ + 2);
        Spots(frame).Fits(spotIndex).FixedAreaIntensity3 = sum(z_raw_values(ismember(z_grid,Spots(frame).Fits(spotIndex).brightestZ-1:Spots(frame).Fits(spotIndex).brightestZ+1)));
        Spots(frame).Fits(spotIndex).FixedAreaIntensity5 = sum(z_raw_values(ismember(z_grid,Spots(frame).Fits(spotIndex).brightestZ-2:Spots(frame).Fits(spotIndex).brightestZ+2)));
        try
            Spots(frame).Fits(spotIndex).cylIntensity = Spots(frame).Fits(spotIndex).cylIntensity(ZStackIndex);
        catch
        end
        
        %use convolution kernel to look for shadows
        z_raw_binary = ~isnan(z_raw_values);
        z_shadow_vec = conv(z_raw_binary,[1 1 1],'same');
        z_shadow_vec = z_shadow_vec(ismember(z_grid,z_vec));
        n_shadows = z_shadow_vec(ZStackIndex)-1;
        
        if n_shadows < num_shadows
            Spots(frame).Fits(spotIndex).discardThis = 1;
            falsePositives = falsePositives + 1;
        else
            Spots2(frame).Fits = [Spots2(frame).Fits, Spots(frame).Fits(spotIndex)];            
        end
    end
end

end