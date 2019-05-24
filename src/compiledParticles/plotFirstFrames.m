function plotFirstFrames(NChannels, HistoneChannel, nc13, nc14, ...
    CompiledParticles, DropboxFolder, Prefix, ElapsedTime, ExperimentAxis)
%plotFirstFramesSummary of this function goes here
%   Detailed explanation goes here
for ChN=1:NChannels
    
    %How does the first frame from schnitzcells compare to the general one set
    %by just looking at the movie? Do this only for nuclei where we have
    %identified the parent nucleus.
    
    if HistoneChannel
        if ~SkipAll
            figure(11)
            clf
            xRange=linspace(13.5,14.5);
            plot(xRange,ones(size(xRange))*nc14,'-k')
            hold on
            for i=1:length(CompiledParticles{ChN})
                if (CompiledParticles{ChN}(i).nc==14)&(CompiledParticles{ChN}(i).PParticle>0) %#ok<*AND2>
                    plot(14*(1+(rand(1)-0.5)/100),CompiledParticles{ChN}(i).NucStart,'.k')
                end
            end
        

            xRange=linspace(12.5,13.5);
            plot(xRange,ones(size(xRange))*nc13,'-k')
            hold on
            for i=1:length(CompiledParticles{ChN})
                if (CompiledParticles{ChN}(i).nc==13)&(CompiledParticles{ChN}(i).PParticle>0)
                    plot(13*(1+(rand(1)-0.5)/100),CompiledParticles{ChN}(i).NucStart,'.k')
                end
            end
            hold off
            try
                ylim([nc13-5,nc14+5])
            end
            set(gca,'XTick',[13,14])
            xlabel('nc')
            ylabel('Frame')
            title('Division time set by eye vs. actual division times of nuclei')
            saveas(gcf,[DropboxFolder,filesep,Prefix,filesep,'Various',filesep,'DivisionTimes_ch',...
                iIndex(ChN,2),'.tif'])



            %Is there any correlation between the first frame and the time of division?


            figure(12)
            subplot(1,2,1)
            hold on
            for i=1:length(CompiledParticles{ChN})
                if (CompiledParticles{ChN}(i).nc==13)&(CompiledParticles{ChN}(i).PParticle>0)
                    plot(CompiledParticles{ChN}(i).NucStart*(1+(rand(1)-0.5)/100),CompiledParticles{ChN}(i).FirstFrame,'.k')
                end
            end
            hold off
            axis square
            xlabel('Nuclear birth (frame)')
            ylabel('First particle frame')
            title('nc13')


            subplot(1,2,2)
            hold on
            for i=1:length(CompiledParticles{ChN})
                if (CompiledParticles{ChN}(i).nc==14)&(CompiledParticles{ChN}(i).PParticle>0)
                    plot(CompiledParticles{ChN}(i).NucStart*(1+(rand(1)-0.5)/100),CompiledParticles{ChN}(i).FirstFrame,'.k')
                end
            end
            hold off
            axis square
            title('nc14')
            xlabel('Nuclear birth (frame)')
            ylabel('First particle frame')
            saveas(gcf,[DropboxFolder,filesep,Prefix,filesep,'Various',filesep,'DivisionVsFirstFrame_ch',...
                iIndex(ChN,2),'.tif'])

            %First frame and AP position
            if strcmpi(ExperimentAxis,'AP')
                figure(13)
                clf
                hold on
                for i=1:length(CompiledParticles{ChN})
                    plot(CompiledParticles{ChN}(i).MeanAP,....
                        ElapsedTime(CompiledParticles{ChN}(i).FirstFrame)-...
                        ElapsedTime(nc14),'.k')
                end
                hold off
                box on
                xlabel('AP position (x/L)')
                ylabel('Particle first frame (min)')
                %         if length(ElapsedTime) > nc14+20
                %             ylim([0,ElapsedTime(nc14+20)-ElapsedTime(nc14)])
                %         elseif (ElapsedTime(end) - ElapsedTime(nc14))>0
                %             ylim([0, ElapsedTime(end) - ElapsedTime(nc14)])
                %         end
                % ES 2014-01-05 Testing early-nc14-only movies
                saveas(gcf,[DropboxFolder,filesep,Prefix,filesep,'Various',filesep,'FirstFrameVsAP_ch',...
                    iIndex(ChN,2),'.tif'])
            end
        end
    end
end
end

