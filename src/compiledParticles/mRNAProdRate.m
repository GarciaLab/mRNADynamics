function CompiledParticles = mRNAProdRate(NChannels, CompiledParticles, ...
    ncFilter, ElapsedTime, SkipFits, DropboxFolder, Prefix)
%MRNAPRODRATE Summary of this function goes here
%   Detailed explanation goes here
for ChN=1:NChannels
    
    %Plot the results from fitting the individual traces. Notice that I'm
    %redoing the fits here using the range given by FrameRange. This is because
    %the fluorescence values were calculated slightly differently such that the
    %offset could vary.
    
    if isfield(CompiledParticles{ChN},'Fit')
        
        %First, find the maximum in intensity over traces for each cycle
        nc13Max=max([CompiledParticles{ChN}(ncFilter(:,end-1)).Fluo]);
        nc14Max=max([CompiledParticles{ChN}(ncFilter(:,end)).Fluo]);
        
        figure(10)
        clf
        
        for i=1:length(CompiledParticles{ChN})
            if ~isempty(CompiledParticles{ChN}(i).Fit)
                
                %Redo the fit and obtain the parameters
                
                FrameRange=CompiledParticles{ChN}(i).Fit.FrameRange;
                
                FramesRangeFilter=ismember(CompiledParticles{ChN}(i).Frame,[FrameRange(1):FrameRange(2)]);
                
                [a, b, sigma_a, sigma_b] = york_fit(ElapsedTime(CompiledParticles{ChN}(i).Frame(FramesRangeFilter)),...
                    CompiledParticles{ChN}(i).Fluo(FramesRangeFilter),...
                    ones(1,sum(FramesRangeFilter))*mean(diff(ElapsedTime))/2,...
                    ones(1,sum(FramesRangeFilter))*CompiledParticles{ChN}(i).FluoError);
                
                
                CompiledParticles{ChN}(i).Fit.Intercept=a;
                CompiledParticles{ChN}(i).Fit.SDIntercept=sigma_a;
                
                CompiledParticles{ChN}(i).Fit.Slope=b;
                CompiledParticles{ChN}(i).Fit.SDSlope=sigma_b;
                
                
                
                if ncFilter(i,end-1)
                    StartFrame=nc13;
                    EndFrame=nc14;
                    MaxFluo=nc13Max;
                elseif ncFilter(i,end)
                    StartFrame=nc14;
                    EndFrame=length(ElapsedTime);
                    MaxFluo=nc14Max;
                end
                
                if ~SkipFits
 
                    xRange=linspace(ElapsedTime(FrameRange(1)),...
                        ElapsedTime(FrameRange(end)));
                    
                    
                    plot(xRange,b*xRange+a,'-k','LineWidth',3)
                    
                    hold on
                    errorbar(ElapsedTime(CompiledParticles{ChN}(i).Frame),CompiledParticles{ChN}(i).Fluo,...
                        ones(1,length(CompiledParticles{ChN}(i).Frame))*CompiledParticles{ChN}(i).FluoError,'.-r')
                    hold off
                    xlabel('Time (min)')
                    ylabel('Fluorescence (au)')
                    title(['Compiled particle ',num2str(i),', nc',num2str(CompiledParticles{ChN}(i).nc)])
                    
                    xlim([ElapsedTime(StartFrame),ElapsedTime(EndFrame)])
                    ylim([0,MaxFluo])
                    
                    saveas(gcf,[DropboxFolder,filesep,Prefix,filesep,'Fits',filesep,'Fit',iIndex(i,3),'-nc',...
                        num2str(CompiledParticles{1}(i).nc),'_ch',iIndex(ChN,2),'.tif'])
                end
            end
        end
        
        close(10)
    end
end
end

