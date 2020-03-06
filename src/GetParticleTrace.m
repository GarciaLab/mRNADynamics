function [Frame, AmpGaussian,Offset, ErrorGauss,optFit,FitType,...
    AmpDog, AmpDogMax, ampdog3, ampdog3Max] =...
    GetParticleTrace(CurrentParticle, Particles, Spots, plotTraceSettings, noSpline)
% This function uses the total intensity mask to calculate the particles
% intensity and subtracts the background obtained from the fit.

% First, get the different intensity values corresponding to this particle.

plotTraceSettings.ErrorIntegral = NaN;
plotTraceSettings.ErrorIntegral3 = NaN;
plotTraceSettings.backGround3 = NaN;
plotTraceSettings.AmpIntegralGauss3D = NaN;
plotTraceSettings.ErrorIntegralGauss3D = NaN;

doSpline = ~noSpline;

ampdog3 = NaN;
ampdog3Max = NaN;
optFit = [];
FitType = [];
ErrorGauss = [];

defaultArea = 109; %109 pixels is the default area when the pixels are assumed to be 212nm x 212 nm AR 9/3/18

warning('off', 'MATLAB:rankDeficientMatrix'); %suppress the spline fitting warnings

for i=1:length(Particles(CurrentParticle).Frame)
    
    Frame(i)=Particles(CurrentParticle).Frame(i);
    spot = Spots(Frame(i)).Fits(Particles(CurrentParticle).Index(i));
    
    %Determine the brightest Z plane of this particle
    zIndex=find(spot.brightestZ == spot.z);
    
    %Offset obtained using Gaussian fitting.
    Offset(i) = double(spot.Offset(zIndex));
    %Intensity obtained by integrating over the area. This has
    %magnitude already has the offset (obtained from the Gaussian fit)
    %subtracted.
    plotTraceSettings.AmpIntegral(i) = double(spot.FixedAreaIntensity(zIndex));
    %Intensity obtained by integrating over the Gaussian fit. This
    %already has subtracted the offset.
    AmpGaussian(i)= double(spot.GaussianIntensity(zIndex));
    %Amplitude of the Gaussian fit. This is an intensity per pixel.
    IntensityMaxGauss(i)= double(spot.CentralIntensity(zIndex));
    %Check to see it multi-slice integration was performed for this set
    fields = fieldnames(spot);
    try
        plotTraceSettings.AmpIntegral3(i) = double(spot.FixedAreaIntensity3);
    catch
        plotTraceSettings.AmpIntegral3(i)= NaN;
    end
    try
        plotTraceSettings.AmpIntegralGauss3D(i)=...
            double(spot.gauss3DIntensityRaw);
    catch
        plotTraceSettings.AmpIntegralGauss3D(i)= NaN;
    end
    if isfield(spot, 'gauss3DIntensitySE')
        try
            plotTraceSettings.ErrorIntegralGauss3D(i) = double(spot.gauss3DIntensitySE);
        catch
            plotTraceSettings.ErrorIntegralGauss3D(i) = NaN;
        end
    else
        plotTraceSettings.ErrorIntegralGauss3D(i) = NaN;
    end
  
    try
        AmpDog(i) =  double(spot.dogFixedAreaIntensity(zIndex));
    catch
        AmpDog(i) = NaN;
    end
    try
        AmpDogMax(i) =  double(spot.DOGIntensity(zIndex));
    catch
        AmpDogMax(i) = NaN;
    end
     try
        ampdog3(i) =  double(spot.ampdog3);
    catch
        ampdog3(i) = NaN;
     end
     try
        ampdog3Max(i) =  double(spot.ampdog3Max);
    catch
        ampdog3Max(i) = NaN;
    end
end

if doSpline
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
            ErrorGauss = [];
            plotTraceSettings.ErrorIntegral = [];
            plotTraceSettings.ErrorIntegral3 = [];
            optFit = [];
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


if exist('OffsetError', 'var')
    %Now, estimate the error in the signal.
    %For the Gaussian fit, we use the average area to get an overall
    %error for the whole trace. We could also just define an error for
    %each data point.
    try
        ErrorGauss=OffsetError*sqrt(2)*...
            mean(double(spot.Area));
    catch
        try
            ErrorGauss=OffsetError*sqrt(2)*...
                mean(double(spot.Area));
        catch
            ErrorGauss=OffsetError*sqrt(2)*defaultArea;
        end
    end
    
    
    %For the Integral, we just use the area of the snippet, which is a
    %constant for all time points.
    if isfield(spot, 'intArea') && ~isempty(spot.intArea)
        intArea = double(spot.intArea(1));
        plotTraceSettings.ErrorIntegral = OffsetError*sqrt(2)*intArea;
        if ~isnan(plotTraceSettings.AmpIntegral3(i))
            plotTraceSettings.backGround3 = 3*Offset*intArea;
            plotTraceSettings.ErrorIntegral3 = OffsetError*sqrt(2)*3*intArea; % since this integration is actually done as a mean over the available slices, multiplying by 5 is definitely wrong. this error should be taken with
            %a grain of salt. AR 9/3/18
        end
    else
        plotTraceSettings.ErrorIntegral = OffsetError*sqrt(2)*defaultArea;
    end
    
end

end
