function [MeanVector,SDVector,NParticles]=AverageTracesNuclei(FrameInfo,Nuclei,NChannels)

%Average and SD over each time point. In order to do this we'll generate a
%cell array with all the values for a given time point

% TraceCell=cell(length(FrameInfo),1);
% 
% for i=1:length(Nuclei)
%     for j=1:length(Nuclei(i).Frames)
%         TraceCell{Nuclei(i).Frames(j)}=[TraceCell{Nuclei(i).Frames(j)},...
%             Nuclei(i).FluoMax(j)];
%     end
% end

if ~exist('NChannels')
    NChannels=1;
end

%Determine the number of protein channels by looking at the dimensionality
%of FluoMax.
if ~isempty(Nuclei)

    for ChN=1:NChannels

        TraceCell=cell(length(FrameInfo),1);

        for i=1:length(Nuclei)
            for j=1:length(Nuclei(i).Frames)
                TraceCell{Nuclei(i).Frames(j)}=[TraceCell{Nuclei(i).Frames(j)},...
                    Nuclei(i).FluoMax(j,ChN)];
            end
        end

        %Average inside each element of TraceCell. I need to be careful with the
        %NaNs
        for i=1:length(TraceCell)
            NanFilter=~isnan(TraceCell{i});
            MeanVector{ChN}(i)=mean(TraceCell{i}(NanFilter));
            SDVector{ChN}(i)=std(TraceCell{i}(NanFilter));
            NParticles{ChN}(i)=sum(NanFilter);
        end
    end

%If there were no nuclei
else
    for i=1:NChannels
        MeanVector{i}=nan(1,length(FrameInfo));
        SDVector{i}=nan(1,length(FrameInfo));
        NParticles{i}=nan(1,length(FrameInfo));
    end
end

%If we only had one channel, then turn it into a vector
if NChannels==1
    MeanVector=MeanVector{1};
    SDVector=SDVector{1};
    NParticles=NParticles{1};
end
