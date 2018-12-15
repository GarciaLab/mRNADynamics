function [Frame,AmpIntegral,AmpIntegral3,AmpIntegral5,AmpGaussian,Offset,...
    ErrorIntegral,ErrorGauss,optFit,FitType,ErrorIntegral3, ErrorIntegral5,backGround3, AmpIntegralGauss3D, ErrorIntegralGauss3D]=...
    GetParticleTrace(CurrentParticle,Particles,Spots)

%function [Frame,AmpIntegral,AmpIntegral3,AmpIntegral5,AmpGaussian,Offset,...
%    ErrorIntegral,ErrorGauss,optFit,FitType,ErrorIntegral3, ErrorIntegral5]=...
%    GetParticleTrace(CurrentParticle,Particles,Spots)
%
%This function uses the total intensity mask to calculate the particles
%intensity and subtracts the background obtained from the fit.
%
% AR 9/3/2018- Currently, AmpIntegral3 is being calculated over a
% cylindrical volume and not an accordion volume as ampintegral5 is. 



%
%First, get the different intensity values corresponding to this particle.

ErrorIntegral = NaN;
ErrorIntegral3 = NaN;
ErrorIntegral5 = NaN;
backGround3 = NaN;

defaultArea = 109; %109 pixels is the default area when the pixels are assumed to be 212nm x 212 nm AR 9/3/18

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
        %Check to see it multi-slice integration was performed for this set
        fields = fieldnames(Spots(Particles(CurrentParticle).Frame(i)).Fits(Particles(CurrentParticle).Index(i)));
        try
            AmpIntegral3(i)=...
            Spots(Particles(CurrentParticle).Frame(i)).Fits(Particles(CurrentParticle).Index(i)).FixedAreaIntensity3;
% %             AmpIntegral3(i)=...
%             Spots(Particles(CurrentParticle).Frame(i)).Fits(Particles(CurrentParticle).Index(i)).cylIntensity;
%               AmpIntegral3(i)=...
%               Spots(Particles(CurrentParticle).Frame(i)).Fits(Particles(CurrentParticle).Index(i)).gauss3DIntensity;    
        catch
            AmpIntegral3(i)= NaN;
        end
        try 
            AmpIntegralGauss3D(i)=...
              Spots(Particles(CurrentParticle).Frame(i)).Fits(Particles(CurrentParticle).Index(i)).gauss3DIntensity;    
        catch
            AmpIntegralGauss3D(i)= NaN;
        end
        try
            AmpIntegral5(i)=...
            Spots(Particles(CurrentParticle).Frame(i)).Fits(Particles(CurrentParticle).Index(i)).FixedAreaIntensity5;
        catch
            AmpIntegral5(i)= NaN;
        end
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
    catch
        ErrorGauss=[];
        ErrorIntegral=[];
        ErrorIntegral3=[];
        ErrorIntegral5 = [];
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
    FitType='mean';
  
    OffFit=mean(double(Offset));
    optFit=OffFit;
    
    %Calculate the error in the offset
    OffsetError=std(OffFit);
end

if exist('OffsetError')
    %Now, estimate the error in the signal.
    %For the Gaussian fit, we use the average area to get an overall
    %error for the whole trace. We could also just define an error for
    %each data point.
    try
        ErrorGauss=OffsetError*sqrt(2)*...
         mean(cell2mat(Spots(Particles(CurrentParticle).Frame(i)).Fits(Particles(CurrentParticle).Index(i)).Area));
    catch
        try
             ErrorGauss=OffsetError*sqrt(2)*...
             mean(Spots(Particles(CurrentParticle).Frame(i)).Fits(Particles(CurrentParticle).Index(i)).Area);
        catch
            ErrorGauss=OffsetError*sqrt(2)*defaultArea;
        end
    end
    try
        ErrorIntegralGauss3D=OffsetError*sqrt(2)*...
         mean(cell2mat(Spots(Particles(CurrentParticle).Frame(i)).Fits(Particles(CurrentParticle).Index(i)).Area));
    catch
        try
             ErrorIntegralGauss3D=OffsetError*sqrt(2)*...
             mean(Spots(Particles(CurrentParticle).Frame(i)).Fits(Particles(CurrentParticle).Index(i)).Area);
        catch
            ErrorIntegralGauss3D=OffsetError*sqrt(2)*defaultArea;
        end
    end
    %For the Integral, we just use the area of the snippet, which is a
    %constant for all time points.
    if isfield(Spots(Particles(CurrentParticle).Frame(i)).Fits(Particles(CurrentParticle).Index(i)), 'intArea') && ~isempty(Spots(Particles(CurrentParticle).Frame(i)).Fits(Particles(CurrentParticle).Index(i)).intArea)
        intArea = Spots(Particles(CurrentParticle).Frame(i)).Fits(Particles(CurrentParticle).Index(i)).intArea(1);
        ErrorIntegral=OffsetError*sqrt(2)*intArea;
        if ~isnan(AmpIntegral3(i))
           backGround3 = 3*Offset*intArea;
           ErrorIntegral3=OffsetError*sqrt(2)*3*intArea;%since this integration is actually done as a mean over the available slices, multiplying by 5 is definitely wrong. this error should be taken with 
           %a grain of salt. AR 9/3/18
        end
        if ~isnan(AmpIntegral5(i))
           ErrorIntegral5=OffsetError*sqrt(2)*5*intArea; %since this integration is actually done as a mean over the available slices, multiplying by 5 is definitely wrong. this error should be taken with 
           %a grain of salt. AR 9/3/18
        end
    else
        ErrorIntegral=OffsetError*sqrt(2)*defaultArea; 
    end
    
else
    ErrorGauss=[];
    ErrorIntegral=[];
    ErrorIntegral3=[];
    ErrorIntegral5 = [];
    ErrorIntegralGauss3D = [];
    optFit=[];
end
