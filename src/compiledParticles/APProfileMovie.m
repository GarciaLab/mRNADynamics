function APProfileMovie(MeanVectorAP, NParticlesAP, MinParticles, numFrames, ...
    APbinID, SDVectorAP, FrameInfo, ElapsedTime, DropboxFolder, Prefix)
%APPROFILEMOVIE Summary of this function goes here
%   Detailed explanation goes here

for ChN=1:NChannels

    APMovieFig = figure();

    MaxValue=max(max(MeanVectorAP{ChN}));
    NParticlesAPFilter=NParticlesAP{ChN}>=MinParticles;

    for i=1:numFrames

        PlotHandle=errorbar(APbinID(NParticlesAPFilter(i,:)),...
            MeanVectorAP{ChN}(i,NParticlesAPFilter(i,:)),SDVectorAP{ChN}(i,NParticlesAPFilter(i,:)),'.-k');
        hold on
        PlotHandle=[PlotHandle,errorbar(APbinID(NParticlesAPFilter(i,:)),...
            MeanVectorAP{ChN}(i,NParticlesAPFilter(i,:)),...
            SDVectorAP{ChN}(i,NParticlesAPFilter(i,:))./sqrt(NParticlesAP{ChN}(i,NParticlesAPFilter(i,:))),'-k')];
        hold off

        if isfield(FrameInfo, 'nc')
            currentNC = num2str(FrameInfo(i).nc);
            iStr = num2str(i);
            elapsedStr = num2str(round(ElapsedTime(i)*10)/10);
            if eval(['nc',currentNC]) > 0
                title(['nc', currentNC,'. Time into nc: ',num2str(round((ElapsedTime(i)-...
                    round(ElapsedTime(eval(['nc',currentNC]))*10)/10))),' min. Total time: ',...
                    elapsedStr,' min (Frame ',iStr,').'])
            else
                title(['nc',currentNC,'. Total time: ',...
                    elapsedStr,' min (Frame ',iStr,').'])
            end
        else
            title(['nc',currentNC,'. Total time: ',...
                elapsedStr,' min (Frame ',iStr,').'])
        end
        xlim([0.1,0.8]);
        ylim([0,MaxValue]);
        xlabel('AP position (x/L)');
        ylabel('Mean fluorescence');

        StandardFigure(PlotHandle,gca)
        saveas(gcf,[DropboxFolder,filesep,Prefix,filesep,'APMovie',filesep,iIndex(i,3),'_ch',iIndex(ChN,2),'.tif']);
    end
end
end

