function [MeanVector,SDVector,NParticles]=AverageTracesNuclei(FrameInfo,Nuclei,NChannels)

%Average and SD over each time point. In order to do this we'll generate a
%cell array with all the values for a given time point

if ~exist('NChannels', 'var')
    NChannels=1;
end

if ~isempty(Nuclei)

    for ChN=1:NChannels

        TraceCell=cell(length(FrameInfo),1);

        for nuc=1:length(Nuclei)
            for frame=1:length(Nuclei(nuc).Frames)
                TraceCell{Nuclei(nuc).Frames(frame)}=[TraceCell{Nuclei(nuc).Frames(frame)},...
                    Nuclei(nuc).FluoMax(frame,ChN)];
            end
        end

        %Average inside each element of TraceCell. I need to be careful with the
        %NaNs
        for frame=1:length(TraceCell)
            NanFilter=~isnan(TraceCell{frame});
            MeanVector{ChN}(frame)=mean(TraceCell{frame}(NanFilter));
            SDVector{ChN}(frame)=std(TraceCell{frame}(NanFilter));
            NParticles{ChN}(frame)=sum(NanFilter);
        end
    end

%If there were no nuclei
else
    for ch=1:NChannels
        MeanVector{ch}=nan(1,length(FrameInfo));
        SDVector{ch}=nan(1,length(FrameInfo));
        NParticles{ch}=nan(1,length(FrameInfo));
    end
end

%If we only had one channel, then turn it into a vector
if NChannels==1
    MeanVector=MeanVector{1};
    SDVector=SDVector{1};
    NParticles=NParticles{1};
end
