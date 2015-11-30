function [Frame,Amp,Off,Off2,Amp2,AmpOld,AmpRaw,Error,optFit1,FitType]=GetParticleTrace(CurrentParticle,Particles,fad)

%This function uses the total intensity mask to calculate the particles
%intensity and subtracts the background obtained from the fit.

for i=1:length(Particles(CurrentParticle).Frame)
        Frame(i)=Particles(CurrentParticle).Frame(i);

        Off(i)=fad.channels(Particles(CurrentParticle).Frame(i)).fits.off(Particles(CurrentParticle).Index(i));

        Off2(i)=fad.channels(Particles(CurrentParticle).Frame(i)).fits.off2(Particles(CurrentParticle).Index(i));


        AmpRaw(i)=sum(sum(double(immultiply(fad.channels(Particles(CurrentParticle).Frame(i)).fits.snippets(:,:,Particles(CurrentParticle).Index(i)),...
            fad.channels(Particles(CurrentParticle).Frame(i)).fits.maskUsedForTotalInt))));


        AmpOld(i)=sum(sum(double(immultiply(fad.channels(Particles(CurrentParticle).Frame(i)).fits.snippets(:,:,Particles(CurrentParticle).Index(i)),...
            fad.channels(Particles(CurrentParticle).Frame(i)).fits.maskUsedForTotalInt))))-...
            fad.channels(Particles(CurrentParticle).Frame(i)).fits.off(Particles(CurrentParticle).Index(i))*...
            sum(sum(fad.channels(Particles(CurrentParticle).Frame(i)).fits.maskUsedForTotalInt));
end




%Do a spline fit to the offset and use it to calculate the actual amplitude

%If we have less the five data points just fit a line.
if length(Frame)>5
    FitType='spline';
    
    nBreaks=max([floor(length(Frame)/3),3]);
    if nBreaks>5
        nBreaks=5;
    end

    try
        optFit1 = adaptiveSplineFit(double(Frame),double(Off),nBreaks);
        Off1Fit=ppval(optFit1,double(Frame));

        %Calculate the error in the offset
        OffsetError=std(Off-Off1Fit);

        Error=OffsetError*sqrt(2)*sum(sum(fad.channels(Particles(CurrentParticle).Frame(i)).fits.maskUsedForTotalInt));

        for i=1:length(Frame)
            Amp(i)=AmpRaw(i)-Off1Fit(i)*sum(sum(fad.channels(Particles(CurrentParticle).Frame(i)).fits.maskUsedForTotalInt));
        end
    catch
        Amp=double(AmpOld);
        Error=std(Off)*sqrt(2)*sum(sum(fad.channels(Particles(CurrentParticle).Frame(i)).fits.maskUsedForTotalInt));
        optFit1=[];
    end

    try

        optFit2 = adaptiveSplineFit(double(Frame),double(Off2),nBreaks);
        Off2Fit=ppval(optFit2,double(Frame));


        for i=1:length(Frame)
            Amp2(i)=AmpRaw(i)-Off2Fit(i)*sum(sum(fad.channels(Particles(CurrentParticle).Frame(i)).fits.maskUsedForTotalInt));
        end
    catch
        Amp2=[];
    end
elseif length(Frame)>2
    FitType='line';
    optFit1=polyfit(double(Frame),double(Off),1);
    

    %Calculate the error in the offset
    OffsetError=std(Off-polyval(optFit1,double(Frame)));
    
    Error=OffsetError*sqrt(2)*sum(sum(fad.channels(Particles(CurrentParticle).Frame(i)).fits.maskUsedForTotalInt));

    
    optFit2=polyfit(double(Frame),double(Off2),1);
    
    
    
    i=1;
    Amp=AmpRaw-polyval(optFit1,double(Frame))*sum(sum(fad.channels(Particles(CurrentParticle).Frame(i)).fits.maskUsedForTotalInt));
    
    
    Amp2=AmpRaw-polyval(optFit2,double(Frame))*sum(sum(fad.channels(Particles(CurrentParticle).Frame(i)).fits.maskUsedForTotalInt));
    
    
    
else
    FitType='mean';
    
    
    Off1Fit=mean(double(Off));
    optFit1=Off1Fit;
    
    
    %Calculate the error in the offset
    OffsetError=std(Off-Off1Fit);

    Error=OffsetError*sqrt(2)*sum(sum(fad.channels(Particles(CurrentParticle).Frame(i)).fits.maskUsedForTotalInt));

    i=1;
    Amp=AmpRaw-Off1Fit*sum(sum(fad.channels(Particles(CurrentParticle).Frame(i)).fits.maskUsedForTotalInt));
    
    
    Off2Fit=nanmean(double(Off2));
    Amp2=AmpRaw-Off2Fit*sum(sum(fad.channels(Particles(CurrentParticle).Frame(i)).fits.maskUsedForTotalInt));
end
