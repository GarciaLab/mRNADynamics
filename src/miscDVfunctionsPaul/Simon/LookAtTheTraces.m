
Mean = MeanVectorAll{1};
MovieTime = 1:length(FrameInfo);%[FrameInfo.Time]./60;
Error = SDVectorAll{1}./NParticlesAll{1};
%shadedErrorBar(Time,Mean,Error,'lineProps',{'Color','k','LineWidth',2})
hold on
maxParticleFluo = nanmax([CompiledParticles{1}.Fluo]);

%%
for p = 1:length(CompiledParticles{1})
    
    particleFrames = CompiledParticles{1}(p).Frame;
    if length(particleFrames) >3 && particleFrames(1)>nc13 && particleFrames(end)<nc14
        firstFrame = particleFrames(1);
        lastFrame = particleFrames(end);
        particleFluo = CompiledParticles{1}(p).Fluo;
        plot(particleFrames,particleFluo,'r-','LineWidth',1.5)
        hold on
        plot([nc13,nc13],[0,maxParticleFluo]);
        plot([nc14,nc14],[0,maxParticleFluo]);
        shadedErrorBar(MovieTime(nc13:nc14),Mean(nc13:nc14),Error(nc13:nc14))
        hold off
        ylim([0 maxParticleFluo])
        waitforbuttonpress
    end
end
hold off

%% nc14
for p = 1:length(CompiledParticles{1})
    
    particleFrames = CompiledParticles{1}(p).Frame;
    if length(particleFrames)>3 && particleFrames(end)>nc14 
        firstFrame = particleFrames(1);
        lastFrame = particleFrames(end);
        particleFluo = CompiledParticles{1}(p).Fluo;
        
        figure(1)
        plot(particleFrames,particleFluo,'ro-','LineWidth',1.5)
        hold on
        plot([nc14,nc14],[0,maxParticleFluo]);
        shadedErrorBar(MovieTime(nc14:end),Mean(nc14:end),Error(nc14:end))
        hold off
        ylim([0 maxParticleFluo])
        saveas(gcf,num2str(p))
        
        figure(2)
        hold on
        plot([nc13,nc13],[0,maxParticleFluo]);
        plot([nc14,nc14],[0,maxParticleFluo]);
        plot(particleFrames,particleFluo,'o-','LineWidth',1.5)
        hold off       
        
        waitforbuttonpress
    end
end
hold off

%% NC13
for p = 1:length(CompiledParticles{1})
    
    particleFrames = CompiledParticles{1}(p).Frame;
    if length(particleFrames) >3 && particleFrames(1)>nc13 && particleFrames(end)<nc14 
        firstFrame = particleFrames(1);
        lastFrame = particleFrames(end);
        particleFluo = CompiledParticles{1}(p).Fluo;
        
        figure(1)
        plot(particleFrames,particleFluo,'ro-','LineWidth',1.5)
        hold on
        plot([nc13,nc13],[0,maxParticleFluo]);
        plot([nc14,nc14],[0,maxParticleFluo]);
        shadedErrorBar(MovieTime(nc13:nc14),Mean(nc13:nc14),Error(nc13:nc14))
        hold off
        ylim([0 maxParticleFluo])
        saveas(gcf,num2str(p))
        
        figure(2)
        hold on
        plot([nc13,nc13],[0,maxParticleFluo],'k');
        plot([nc14,nc14],[0,maxParticleFluo],'k');
        plot(particleFrames,particleFluo,'o-','LineWidth',1.5)
        hold off       
        
        waitforbuttonpress
    end
end
hold off
    