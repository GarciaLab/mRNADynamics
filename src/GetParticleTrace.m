function [Frame,AmpIntegral,AmpIntegral3,AmpIntegral5,AmpGaussian,Offset,...
    ErrorIntegral,ErrorGauss,optFit,FitType,ErrorIntegral3, ErrorIntegral5,backGround3,...
    AmpIntegralGauss3D, ErrorIntegralGauss3D, AmpDog, AmpDogMax, ampdog3, ampdog3Max] =...
    ...
    GetParticleTrace(...
    ...
    CurrentParticle,Particles,Spots)

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
Spots = castStructNumbersToDoubles(Spots);

ErrorIntegral = NaN;
ErrorIntegral3 = NaN;
ErrorIntegral5 = NaN;
backGround3 = NaN;
AmpIntegralGauss3D = NaN;
ErrorIntegralGauss3D = NaN;
ampdog3 = NaN;
ampdog3Max = NaN;

defaultArea = 109; %109 pixels is the default area when the pixels are assumed to be 212nm x 212 nm AR 9/3/18

warning('off', 'MATLAB:rankDeficientMatrix'); %suppress the spline fitting warnings

for i=1:length(Particles(CurrentParticle).Frame)
    
    Frame(i)=Particles(CurrentParticle).Frame(i);
    spot = Spots(Frame(i)).Fits(Particles(CurrentParticle).Index(i));
    
    %Determine the brightest Z plane of this particle
    zIndex=find(spot.brightestZ == spot.z);
    
    %Offset obtained using Gaussian fitting.
    Offset(i) = spot.Offset(zIndex);
    %Intensity obtained by integrating over the area. This has
    %magnitude already has the offset (obtained from the Gaussian fit)
    %subtracted.
    AmpIntegral(i) = spot.FixedAreaIntensity(zIndex);
    %Intensity obtained by integrating over the Gaussian fit. This
    %already has subtracted the offset.
    AmpGaussian(i)= spot.GaussianIntensity(zIndex);
    %Amplitude of the Gaussian fit. This is an intensity per pixel.
    IntensityMaxGauss(i)= spot.CentralIntensity(zIndex);
    %Check to see it multi-slice integration was performed for this set
    fields = fieldnames(spot);
    try
        AmpIntegral3(i) = spot.FixedAreaIntensity3;
    catch
        AmpIntegral3(i)= NaN;
    end
    try
        AmpIntegralGauss3D(i)=...
            spot.gauss3DIntensity;
    catch
        AmpIntegralGauss3D(i)= NaN;
    end
    if isfield(spot, 'gauss3DIntensityCI95')
        try
            ErrorIntegralGauss3D(i) = spot.gauss3DIntensityCI95;
        catch
            ErrorIntegralGauss3D(i) = NaN;
        end
    else
        ErrorIntegralGauss3D(i) = NaN;
%         warning('gauss3d intensities calculated but not their errors. Re-run fit3dgaussianstoallspots if this is desired.');
    end
    try
        AmpIntegral5(i)=...
            spot.FixedAreaIntensity5;
    catch
        AmpIntegral5(i)= NaN;
    end
    try
        AmpDog(i) =  spot.dogFixedAreaIntensity(zIndex);
    catch
        AmpDog(i) = NaN;
    end
    try
        AmpDogMax(i) =  spot.DOGIntensity(zIndex);
    catch
        AmpDogMax(i) = NaN;
    end
     try
        ampdog3(i) =  spot.ampdog3;
    catch
        ampdog3(i) = NaN;
     end
     try
        ampdog3Max(i) =  spot.ampdog3Max;
    catch
        ampdog3Max(i) = NaN;
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
            mean(spot.Area);
    catch
        try
            ErrorGauss=OffsetError*sqrt(2)*...
                mean(spot.Area);
        catch
            ErrorGauss=OffsetError*sqrt(2)*defaultArea;
        end
    end
    %      ErrorIntegralGauss3D=OffsetError*sqrt(2)*...
    %          mean(spot.Area);
    
    
    %For the Integral, we just use the area of the snippet, which is a
    %constant for all time points.
    if isfield(spot, 'intArea') && ~isempty(spot.intArea)
        intArea = spot.intArea(1);
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
