function [ncFilterID, ncFilter, APFilter, APFilter_ROI, APFilter_nonROI, ...
    DVFilter, DVFilter_ROI, DVFilter_nonROI]...
    ...
    = binParticles(nc9, nc10, nc11, nc12,...
    ...
    nc13, nc14, NChannels, CompiledParticles, ExperimentAxis, ROI,...
    APbinID, DVbinID, CompiledParticles_ROI, CompiledParticles_nonROI)

%CREATEFILTERS Summary of this function goes here
%   Detailed explanation goes here

ncFilter = [];
APFilter = cell(1, NChannels);
APFilter_ROI = cell(1, NChannels); APFilter_nonROI = cell(1, NChannels);
DVFilter = cell(1, NChannels); DVFilter_ROI = cell(1, NChannels);  DVFilter_nonROI = cell(1, NChannels);


if ~isnan(nc9)|~isnan(nc10)|~isnan(nc11)|~isnan(nc12)|~isnan(nc13)|~isnan(nc14)
    %ncFilterID just tells you the identity of the different
    %filters stored in the cell ncFilter
    ncFilterID=[];
    if nc9~=0
        ncFilterID=9;
    end
    if nc10~=0
        ncFilterID=[ncFilterID,10];
    end
    if nc11~=0
        ncFilterID=[ncFilterID,11];
    end
    if nc12~=0
        ncFilterID=[ncFilterID,12];
    end
    if nc13~=0
        ncFilterID=[ncFilterID,13];
    end
    if nc14~=0
        ncFilterID=[ncFilterID,14];
    end
    %Add the first nc
    ncFilterID=[min(ncFilterID)-1,ncFilterID];
    
    
    %Create the filters for the CompiledParticles
    for ChN=1:NChannels
        if isempty(CompiledParticles{ChN})
            warning(['No compiled particles found in channel ',num2str(ChN),'. Did you mean to run the code with ApproveAll?',...
                'Also try running CheckParticleTracking again and saving the results. If that doesn''t work and you think there should be particles, talk to HG.'])
        else
            ncFilter=false(length(CompiledParticles{ChN})...
                ,length(ncFilterID)); %AR 6/16/17: I think multi-channel data might require this to be a cell? Something for the future.
            for i=1:length(CompiledParticles{ChN})
                %Sometimes CompiledParticles{1}(i).nc is empty. This is because of some
                %problem with FrameInfo! In that case we'll pull the information out of
                %the XLS file.
                if ~isfield(CompiledParticles{ChN}(i), 'nc')
                    CompiledParticles{ChN}(i).nc = [];
                end
                if ~isempty(CompiledParticles{ChN}(i).nc)
                    ncFilter(i,CompiledParticles{ChN}(i).nc==ncFilterID)=true;
                else
                    ncsFound=find(CompiledParticles{ChN}(i).Frame(1)>=[nc9,nc10,nc11,nc12,nc13,nc14]);
                    if ncsFound(end)==1
                        CompiledParticles{ChN}(i).nc=9;
                        ncFilter(i,ncFilterID==9)=true;
                    elseif ncsFound(end)==2
                        CompiledParticles{ChN}(i).nc=10;
                        ncFilter(i,ncFilterID==10)=true;
                    elseif ncsFound(end)==3
                        CompiledParticles{ChN}(i).nc=11;
                        ncFilter(i,ncFilterID==11)=true;
                    elseif ncsFound(end)==4
                        CompiledParticles{ChN}(i).nc=12;
                        ncFilter(i,ncFilterID==12)=true;
                    elseif ncsFound(end)==5
                        CompiledParticles{ChN}(i).nc=13;
                        ncFilter(i,ncFilterID==13)=true;
                    elseif ncsFound(end)==6
                        CompiledParticles{ChN}(i).nc=14;
                        ncFilter(i,ncFilterID==14)=true;
                    end
                    
                end
            end
            
            %AP filters:
            if strcmpi(ExperimentAxis,'AP') || strcmpi(ExperimentAxis,'DV')
                %Divide the AP axis into boxes of a certain AP size. We'll see which
                %particle falls where.
                
                
                if ROI
                    %Define two APFilters for ROI and non-ROI respectively
                    APFilter_ROI{ChN}=false(length(CompiledParticles_ROI{ChN}),length(APbinID));
                    APFilter_nonROI{ChN}=false(length(CompiledParticles_nonROI{ChN}),length(APbinID));
                    APFilter{ChN}=false(length(CompiledParticles{ChN}),length(APbinID));
                    
                    for i=1:length(CompiledParticles{ChN})
                        APFilter{ChN}(i,find(APbinID<=CompiledParticles{ChN}(i).MeanAP, 1, 'last' ))=1;
                    end
                    
                    for i=1:length(CompiledParticles_ROI{ChN})
                        APFilter_ROI{ChN}(i,find(APbinID<=CompiledParticles_ROI{ChN}(i).MeanAP, 1, 'last' ))=1;
                    end
                    
                    for i=1:length(CompiledParticles_nonROI{ChN})
                        APFilter_nonROI{ChN}(i,find(APbinID<=CompiledParticles_nonROI{ChN}(i).MeanAP, 1, 'last' ))=1;
                    end
                    
                end
                
                
                APFilter{ChN}=false(length(CompiledParticles{ChN}),length(APbinID));
                for i=1:length(CompiledParticles{ChN})
                    APFilter{ChN}(i,find(APbinID<=CompiledParticles{ChN}(i).MeanAP, 1, 'last' ))=1;
                end
                
            end
            
            %DV filters:
            if strcmpi(ExperimentAxis,'DV')
                %Divide the DV axis into boxes of a certain DV size. We'll see which
                %particle falls where.
                
                if ROI
                    %Define two DVFilters for ROI and non-ROI respectively
                    DVFilter_ROI{ChN}=false(length(CompiledParticles_ROI{ChN}),length(DVbinID));
                    DVFilter_nonROI{ChN}=false(length(CompiledParticles_nonROI{ChN}),length(DVbinID));
                    DVFilter{ChN}=false(length(CompiledParticles{ChN}),length(DVbinID));
                    
                    for i=1:length(CompiledParticles{ChN})
                        DVFilter{ChN}(i,max(find(DVbinID<=abs(CompiledParticles{ChN}(i).MeanDV))))=1; %JAKE: Added abs for DV
                    end
                    
                    for i=1:length(CompiledParticles_ROI{ChN})
                        DVFilter_ROI{ChN}(i,max(find(DVbinID<=abs(CompiledParticles_ROI{ChN}(i).MeanDV))))=1; %JAKE: Added abs for DV
                    end
                    
                    for i=1:length(CompiledParticles_nonROI{ChN})
                        DVFilter_nonROI{ChN}(i,max(find(DVbinID<=abs(CompiledParticles_nonROI{ChN}(i).MeanDV))))=1; %JAKE: Added abs for DV
                    end
                    
                end
                
                
                DVFilter{ChN}=false(length(CompiledParticles{ChN}),length(DVbinID));
                
                for i=1:length(CompiledParticles{ChN})
                    DVFilter{ChN}(i,max(find(DVbinID<=abs(CompiledParticles{ChN}(i).MeanDV)))) = 1;
                end
                
                
            end
        end
    end
end

end
