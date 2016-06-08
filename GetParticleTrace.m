function [Frame,AmpIntegral,AmpGaussian,Offset,...
    ErrorIntegral,ErrorGauss,optFit,FitType]=...
    GetParticleTrace(CurrentParticle,Particles,Spots)

%This function uses the total intensity mask to calculate the particles
%intensity and subtracts the background obtained from the fit.

%First, get the different intensity values corresponding to this particle.
for i=1:length(Particles(CurrentParticle).Frame)
        Frame(i)=Particles(CurrentParticle).Frame(i);
        
        %Determine the brightest Z plane of this particle
        zIndex=find(Spots(Particles(CurrentParticle).Frame(i)).Fits(Particles(CurrentParticle).Index(i)).brightestZ==...
            Spots(Particles(CurrentParticle).Frame(i)).Fits(Particles(CurrentParticle).Index(i)).z);
        
        %Offset obtained using Gaussian fitting.
        Offset(i)=...
            Spots(Particles(CurrentParticle).Frame(i)).Fits(Particles(CurrentParticle).Index(i)).Offset(zIndex);
        %Intensity obtained by integrating over the area. This has
        %magnitude already has the offset (obtained from the Gaussian fit)
        %subtracted.
        AmpIntegral(i)=...
            Spots(Particles(CurrentParticle).Frame(i)).Fits(Particles(CurrentParticle).Index(i)).FixedAreaIntensity(zIndex);
        %Intensity obtained by integrating over the Gaussian fit. This
        %already has subtracted the offset.
        AmpGaussian(i)=...
            Spots(Particles(CurrentParticle).Frame(i)).Fits(Particles(CurrentParticle).Index(i)).GaussianIntensity(zIndex);
        %Amplitude of the Gaussian fit. This is an intensity per pixel.
        IntensityMaxGauss(i)=...
            Spots(Particles(CurrentParticle).Frame(i)).Fits(Particles(CurrentParticle).Index(i)).CentralIntensity(zIndex);
end


%Do a spline fit to the offset and use it to estimate the error

%If we have more than five data points, we fit a spline.
if length(Frame)>5
    FitType='spline';
    
    nBreaks=max([floor(length(Frame)/3),3]);
    if nBreaks>5
        nBreaks=5;
    end

    try
        optFit = adaptiveSplineFit(double(Frame),double(Offset),nBreaks);
        OffsetFit=ppval(optFit,double(Frame));

        %Calculate the error in the offset
        OffsetError=std(Offset-OffsetFit);

        %Now, estimate the error in the signal.
        %For the Gaussian fit, we use the average area to get an overall
        %error for the whole trace. We could also just define an error for
        %each data point.
        ErrorGauss=OffsetError*sqrt(2)*...
            cellfun(@mean,Spots(Particles(CurrentParticle).Frame(i)).Fits(Particles(CurrentParticle).Index(i)).Area);
        %For the Integral, we just use the area of the snippet, which is a
        %constant for all time points.
        ErrorIntegral=OffsetError*sqrt(2)*...
            sum(sum(Spots(Particles(CurrentParticle).Frame(i)).Fits(Particles(CurrentParticle).Index(i)).snippet_mask{1}));
    catch
        ErrorGauss=[];
        ErrorIntegral=[];
        optFit=[];
    end
    
%If we have between 3 and five data points, we fit a line.
elseif length(Frame)>2
      
    FitType='line';
    optFit = polyfit(double(Frame),double(Offset),1);

    %Calculate the error in the offset
    OffsetError=std(Offset-polyval(optFit,double(Frame)));
    
%If we have less than 2 data points, we just take the mean.    
else
    error('HG: Adapt this')
    FitType='mean';
  
    OffFit=mean(double(Offset));
    optFit=OffFit;
    
    %Calculate the error in the offset
    OffsetError=std(OffFit-optFit);
end

if exist('OffsetError')
    %Now, estimate the error in the signal.
    %For the Gaussian fit, we use the average area to get an overall
    %error for the whole trace. We could also just define an error for
    %each data point.
    ErrorGauss=OffsetError*sqrt(2)*...
        cellfun(@mean,Spots(Particles(CurrentParticle).Frame(i)).Fits(Particles(CurrentParticle).Index(i)).Area);
    %For the Integral, we just use the area of the snippet, which is a
    %constant for all time points.
    ErrorIntegral=OffsetError*sqrt(2)*...
        sum(sum(Spots(Particles(CurrentParticle).Frame(i)).Fits(Particles(CurrentParticle).Index(i)).snippet_mask{1}));
else
    ErrorGauss=[];
    ErrorIntegral=[];
    optFit=[];
end
