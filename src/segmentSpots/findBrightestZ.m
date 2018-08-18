function [Particles,falsePositives] = findBrightestZ(Particles, num_shadows, use_integral_center, force_z)
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

    falsePositives = 0;
    for i = 1:length(Particles)
        
        z_vec = [Particles(i).z]; %convenience vector
        %pull intensity value from particle snippets
        RawIntensityVec = [Particles(i).FixedAreaIntensity];            
        CentralIntensityVec = [Particles(i).CentralIntensity]; 
        %find slice with brightest pixel
        [~, MaxIndexCentral] = max(CentralIntensityVec);            
        % calculate convenience vectors
        z_grid = min(z_vec):max(z_vec);
        z_raw_values = NaN(size(z_grid));            
        z_raw_values(ismember(z_grid,z_vec)) = RawIntensityVec;
        if ~use_integral_center                
            CentralZ = z_vec(MaxIndexCentral); 
            ZStackIndex = MaxIndexCentral;
        else
            % Convolve with gaussian filter to find "best" center
            g = [-1 0 1];
            gaussFilter = exp(-g .^ 2 / (2 ));
            RawRefVec = conv(gaussFilter,z_raw_values);
            RawRefVec = RawRefVec(2:end-1);
            RawRefVec(1) = NaN;
            RawRefVec(end) = NaN;
            RawRefVec = RawRefVec(ismember(z_grid,z_vec));
            [~, MaxIndexIntegral] = max(RawRefVec);
            CentralZ = z_vec(MaxIndexIntegral);               
            ZStackIndex = MaxIndexIntegral;
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
        Particles(i).FixedAreaIntensity3 = nansum(z_raw_values(ismember(z_grid,Particles(i).brightestZ-1:Particles(i).brightestZ+1)));
        Particles(i).FixedAreaIntensity5 = nansum(z_raw_values(ismember(z_grid,Particles(i).brightestZ-2:Particles(i).brightestZ+2)));
        
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


end