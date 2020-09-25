function CorrectDivisionTimes(Prefix, varargin)

%The idea is to determine division times with a higher spatial resolution
%by doing it per AP bin.

%m n: Move between nuclear cycles
%, .: Move between frames
%Click: Division of clicked AP bin in current frame
%r  : Reset the information for the current nuclear cycle
%s  : Save the information
%x  : Save and quit


%Find out which computer this is. That will determine the folder structure.
%Information about about folders
%% 

close all
resetAnaphaseFrames = false;
if ~isempty(varargin)
    if lower(varargin{1}) == 'resetanaphaseframes'
        
        resetAnaphaseFrames = true;
    else
        error('Unrecognized flag passed to function')
    end
end
liveExperiment = LiveExperiment(Prefix);
FrameInfo = getFrameInfo(liveExperiment); 
schnitzcells = getSchnitzcells(liveExperiment);
schnitzAnaphaseFrames = zeros(1, length(schnitzcells));
if ~isfield(schnitzcells, 'anaphaseFrame') | resetAnaphaseFrames
    for i=1:length(schnitzcells)
        schnitzcells(i).anaphaseFrame = [];
    end
else
    for i=1:length(schnitzcells)
        if ~isempty(schnitzcells(i).anaphaseFrame)
            schnitzAnaphaseFrames(i) = schnitzcells(i).anaphaseFrame;
        end
    end
end

if ~isfield(schnitzcells, 'inferredAnaphaseFrame') | resetAnaphaseFrames
    for i=1:length(schnitzcells)
        schnitzcells(i).inferredAnaphaseFrame = false;
    end
end
Ellipses = getEllipses(liveExperiment);
broken = checkSchnitzCellErrors(schnitzcells, Ellipses);
if ~isempty(find(broken ~= 0, 1))
    [schnitzcells, Ellipses] = correctSchnitzCellErrors(schnitzcells, Ellipses);
end

[SourcePath, FISHPath, DefaultDropboxFolder, DropboxFolder, MS2CodePath, PreProcPath,...
    configValues, movieDatabasePath] = DetermineAllLocalFolders(Prefix);



try 
    hisMat = getHisMat(liveExperiment); 
    ZoomImage = hisMat(:,:,end);
    nFrames = size(hisMat, 3);
catch 
    D=dir([PreProcPath,filesep,Prefix,filesep,Prefix,'-His*.tif']);
    ZoomImage=imread([PreProcPath,filesep,Prefix,filesep,D(end).name]);
    nFrames = length(D);
end






%Get the information about the AP axis as well as the image shifts
%used for the stitching of the two halves of the embryo
load([DropboxFolder,filesep,Prefix,filesep,'APDetection.mat'])


if ~exist('coordPZoom', 'var')
    warning('AddParticlePosition should have been run first. Running it now.')
    AddParticlePosition(Prefix, 'ManualAlignment')
    load([DropboxFolder,filesep,Prefix,filesep,'APDetection.mat'])
end


%Angle between the x-axis and the AP-axis
APAngle=atan2((coordPZoom(2)-coordAZoom(2)),(coordPZoom(1)-coordAZoom(1)));
DVAngle = APAngle + pi/2;
APLength=sqrt((coordPZoom(2)-coordAZoom(2))^2+(coordPZoom(1)-coordAZoom(1))^2);
% Slope and intercept
m = (coordPZoom(2)-coordAZoom(2))/(coordPZoom(1)-coordAZoom(1));
b = coordA(2)-coordA(1)*m;


%% 

APPosImage=zeros(size(ZoomImage));
DVPosImage=zeros(size(ZoomImage));
[Rows,Columns]=size(ZoomImage);

for i=1:Rows
    for j=1:Columns
        Angle = atan2((i-coordAZoom(2)),(j-coordAZoom(1)));
        APDistance=sqrt((coordAZoom(2)-i).^2+(coordAZoom(1)-j).^2);
        APPosition=APDistance.*cos(Angle-APAngle);
        DVPosition = (-m*j + i - b )/sqrt(m^2 + 1);
        APPosImage(i,j)=APPosition/APLength;
        DVPosImage(i,j)=DVPosition/APLength;
    end
end
%% 

[DateFromDateColumn, ExperimentType, ExperimentAxis, CoatProtein, StemLoop, APResolution,...
    Channel1, Channel2, Objective, Power, DataFolder, DropboxFolderName, Comments,...
    nc9, nc10, nc11, nc12, nc13, nc14, CF] = getExperimentDataFromMovieDatabase(Prefix, DefaultDropboxFolder);

ncs=[nc9,nc10,nc11,nc12,nc13,nc14];

APbinID=0:APResolution:1;
DVBinMinIdx = round(min(min(DVPosImage))/APResolution)-1;
DVBinMin = DVBinMinIdx*APResolution;
DVBinMaxIdx = round(max(max(DVPosImage))/APResolution)+1;
DVBinMax = DVBinMaxIdx*APResolution;
DVbinID = DVBinMin:APResolution:DVBinMax;

APPosBinImage=zeros(size(APPosImage));
DVPosBinImage = zeros(size(DVPosImage));
for i=1:(length(APbinID)-1)
    FilteredMask=(APbinID(i)<=APPosImage)&(APbinID(i+1)>APPosImage);
    
    APPosBinImage=APPosBinImage+FilteredMask*i;
end
for j=1:(length(DVbinID)-1)
    FilteredMask=(DVbinID(j)<=DVPosImage)&(DVbinID(j+1)>DVPosImage);
    DVPosBinImage=DVPosBinImage+FilteredMask*j;
end


%Load the division information if it's already there
if ~resetAnaphaseFrames
    if exist([DropboxFolder,filesep,Prefix,filesep,'GridDivision.mat'], 'file')
        load([DropboxFolder,filesep,Prefix,filesep,'GridDivision.mat'], 'GridDivision')
        %Check if we changed the number of AP bins
        if (size(GridDivision,2)~=length(APbinID)) | ((size(GridDivision, 3) ~= length(DVbinID)))
            GridDivision=zeros(14,length(APbinID), length(DVbindID));
        end
    else
        %Matrix where we'll store the information about the divisions
        GridDivision=zeros(14,length(APbinID), length(DVbinID));
    end
else
    GridDivision=zeros(14,length(APbinID), length(DVbinID));
end



%Show the frames on top of the AP bins
figureOverlay=figure;
axOverlay = axes(figureOverlay);


try
    CurrentFrame = min(ncs(ncs~=0));
catch
    CurrentFrame=1;
end

CurrentNC=find(ncs,1)+8;

cc=1;

%% 

APmin = min(min(APPosBinImage));
APmax = max(max(APPosBinImage));
DVmin = min(min(DVPosBinImage));
DVmax = max(max(DVPosBinImage));
APboxes = [];
DVboxes = [];

for apb=1:length(APbinID)-1
    
    DVPosBinSubset = DVPosBinImage(APPosBinImage == apb);
    dvb_subset = unique(unique(DVPosBinSubset));
    for j=1:length(dvb_subset)
        dvb = dvb_subset(j);
        counter = length(APboxes) + 1;
        APboxes(counter) = apb;
        DVboxes(counter) = dvb;
    end
end



%% 



SelectionMode = 1;
SelectionLabels = {'AP', 'DV', 'AP/DV', 'Ellipses', 'Unselect Ellipses'};
debug_mode = false;

%% 

while (cc~='x')
    
    if length(axOverlay.Children) > 1
        delete(axOverlay.Children(1:end-1))
    end
    
    
    
    
    %Generate the image with the information about divisions
    %Green: This AP bin divided in the current frame
    %Blue: This AP bin has been determined, but divided in another frame
    %than the current one
    %Red: The division of this AP bin has not been determined.
    
    BlueImage=zeros(size(APPosBinImage));
    GreenImage=zeros(size(APPosBinImage));
    RedImage=zeros(size(APPosBinImage));
    
    for i=1:length(APbinID)
        for j = 1:length(DVbinID)
            if GridDivision(CurrentNC,i, j)
                if GridDivision(CurrentNC,i,j)==CurrentFrame
                    GreenImage((APPosBinImage==i) & (DVPosBinImage == j))=i+j;
                else
                    BlueImage((APPosBinImage==i) & (DVPosBinImage == j))=i+j;
                end
            else
                RedImage((APPosBinImage==i) & (DVPosBinImage == j))=i+j;
            end
        end
    end
    
    unselectedSchnitzCells = [];
    currentFrameSchnitzCells = [];
    otherFrameSchnitzCells = [];

    for i=1:length(schnitzcells)
        frame_idx = find(schnitzcells(i).frames == CurrentFrame, 1);
        
        if ~isempty(frame_idx)
            if schnitzAnaphaseFrames(i) == CurrentFrame
                currentFrameSchnitzCells = [currentFrameSchnitzCells, schnitzcells(i).cellno(frame_idx)];
            elseif schnitzAnaphaseFrames(i) ~= 0
                otherFrameSchnitzCells = [otherFrameSchnitzCells, schnitzcells(i).cellno(frame_idx)];
            else
                unselectedSchnitzCells= [unselectedSchnitzCells,schnitzcells(i).cellno(frame_idx)];
            end
        end
        
    end
    if ~isnan(CurrentFrame)
        try
            HisImage = hisMat(:, :, CurrentFrame);
        catch
            HisImage=imread([PreProcPath,filesep,Prefix,filesep,D(CurrentFrame).name]);
        end
    end
    
    
    
    %Combine the images
    HisOverlay=cat(3,mat2gray(HisImage)+mat2gray(RedImage),...
        mat2gray(HisImage)+mat2gray(GreenImage)/2,...
        mat2gray(HisImage)+mat2gray(BlueImage));
    if debug_mode
%         for i=1:length(schnitzcells)
%             frame_idx = find(schnitzcells(i).frames == CurrentFrame, 1);
%             if ~isempty(frame_idx)
%                 if schnitzAnaphaseFrames(i) > 0
%                     HisOverlay = insertText(HisOverlay,...
%                         [schnitzcells(i).cenx(frame_idx),schnitzcells(i).ceny(frame_idx)],...
%                         [num2str(schnitzAnaphaseFrames(i))],'FontSize',14 ,'TextColor','white', 'BoxOpacity', 0);
%                 end
%             end
% 
%         end

        for i=1:length(APbinID)-1
            for j = 1:length(DVbinID)-1
                [ypos, xpos] = find((APPosBinImage == i) & (DVPosBinImage == j));
                if ~isempty(xpos)
                    middlex = min(xpos)+(max(xpos)-min(xpos))/4;
                    middley = min(ypos)+(max(ypos)-min(ypos))/4;
                    HisOverlay = insertText(HisOverlay,...
                        [middlex,middley],...
                        [num2str(i), ', ', num2str(j), ', ', num2str(GridDivision(CurrentNC,i,j))],...
                        'FontSize',14 ,'TextColor','black', 'BoxOpacity', 0);

                end
            end
        end
    end
    if isempty(axOverlay.Children)
        hisImageHandle = imshow(HisOverlay, 'Parent', axOverlay);
    else
        hisImageHandle.CData = HisOverlay;
    end
    
    hold(axOverlay, 'on')
    
    
   
    
    

    EllipseHandle = plotEllipses(Ellipses, CurrentFrame, unselectedSchnitzCells, 'w', 20, axOverlay);
    EllipseGreenHandle = plotEllipses(Ellipses, CurrentFrame, currentFrameSchnitzCells, 'y', 20, axOverlay);
    EllipseBlueHandle = plotEllipses(Ellipses, CurrentFrame, otherFrameSchnitzCells, 'b', 20, axOverlay);

    

    set(figureOverlay,'Name',(['Frame: ',num2str(CurrentFrame),'/',num2str(nFrames),...
        '. Current nc:',num2str(CurrentNC)]));
    title(['Selection Mode: ', SelectionLabels{SelectionMode}])
    hold(axOverlay, 'off')
    %     figure(Overlay)
    ct=waitforbuttonpress;
    cc=get(figureOverlay,'currentcharacter');
    cm=get(axOverlay,'CurrentPoint');
    
    %Move frames
    if (ct~=0)&(cc=='.')&(CurrentFrame<nFrames)
        CurrentFrame=CurrentFrame+1;
    elseif (ct~=0)&(cc==',')&(CurrentFrame>1)
        CurrentFrame=CurrentFrame-1;
    elseif (ct~=0)&(cc=='>')& (CurrentFrame+10)< nFrames
        CurrentFrame=CurrentFrame+10;
    elseif (ct~=0)&(cc=='<')&( CurrentFrame-10) > 1
        CurrentFrame=CurrentFrame-10;
        %Move nc
    elseif (ct~=0)&(cc=='m')&(CurrentNC<14) & ~isnan(ncs(CurrentNC+1-8))
        NCTest = [];
        for k=1:length(APboxes)
            NCTest(k) = GridDivision(CurrentNC, APboxes(k), DVboxes(k));
        end
        if min(NCTest) > 0
            CurrentNC=CurrentNC+1;
            eval(['CurrentFrame=nc',num2str(CurrentNC)]);
        else
            incomplete_bins = find(NCTest == 0);
            disp('Finish selecting anaphase frames for the current nuclear cycle.')
            for k =1:length(incomplete_bins)
                disp(['AP bin: ', num2str(APboxes(incomplete_bins(k))), ', DV bin: ', num2str(DVboxes(incomplete_bins(k)))]);
            end
            prompt = ['Do you want to autofill remaining regions (y/n)? '];
            ID = input(prompt,'s');
            if ID == 'y' 
               assigned_rows = [];
               assigned_columns = [];
               associated_anaphase = [];
               for k=1:length(APboxes)
                    if GridDivision(CurrentNC, APboxes(k), DVboxes(k)) > 0
                      [binrows, bincolumns] = find((APPosBinImage ==  APboxes(k)) &  ...
                          (DVPosBinImage ==  DVboxes(k)));
                      assigned_rows = [assigned_rows; binrows];
                      assigned_columns = [assigned_columns; bincolumns];
                      associated_anaphase = [associated_anaphase; GridDivision(CurrentNC, APboxes(k), DVboxes(k))*ones(size(binrows))];
                    end
               end
               for l=1:length(incomplete_bins)
                   idx = incomplete_bins(l);
                   closest_anaphase_frames = [];
                   [current_binrows, current_bincolumns] = find((APPosBinImage ==  APboxes(idx)) &  ...
                          (DVPosBinImage ==  DVboxes(idx)));
                    for m=1:length(current_binrows)
                        Distances=sqrt((assigned_rows-current_binrows(m)).^2+...
                            (assigned_columns-current_bincolumns(m)).^2);
                        closest_anaphase_frames(m) = associated_anaphase(find(Distances == min(Distances), 1));
                    end
                    GridDivision(CurrentNC, APboxes(idx), DVboxes(idx)) = round(median(closest_anaphase_frames));
               end
               CurrentNC=CurrentNC+1;
               eval(['CurrentFrame=nc',num2str(CurrentNC)]);
            end
        end
    elseif (ct~=0)&(cc=='n')&(CurrentNC>8)&eval(['nc',num2str(CurrentNC-1),'~=0'])
        CurrentNC=CurrentNC-1;
        eval(['CurrentFrame=nc',num2str(CurrentNC)]);
    elseif (ct~=0)&(cc=='d')
        if debug_mode
            debug_mode = false;
        else
            debug_mode = true;
        end
    elseif (ct~=0)&(cc=='p')
        SelectionMode = SelectionMode + 1;
        if SelectionMode > 5
            SelectionMode = 1;
        end
    elseif (ct~=0)&(cc=='o')
        SelectionMode = SelectionMode - 1;
        if SelectionMode < 1
            SelectionMode = 5;
        end
        
        %Reset the information
    elseif (ct~=0)&(cc=='r')
        GridDivision(CurrentNC,:,:)=0;
        
        %Save
    elseif (ct~=0)&(cc=='s')
        save([DropboxFolder,filesep,Prefix,filesep,'GridDivision.mat'],'GridDivision')
        disp('GridDivision.mat saved.');
        %Select a time for division
    elseif (ct==0)&(strcmp(get(figureOverlay,'SelectionType'),'normal'))
        cc=1;
        if (cm(1,2)>0)&(cm(1,1)>0)&(cm(1,2)<=Rows)&(cm(1,1)<=Columns)
            if SelectionMode == 3
                GridDivision(CurrentNC,APPosBinImage(round(cm(1,2)),round(cm(1,1))),...
                    DVPosBinImage(round(cm(1,2)),round(cm(1,1))))=CurrentFrame;
            elseif SelectionMode == 1
                dv_set_range = find(GridDivision(CurrentNC,APPosBinImage(round(cm(1,2)),round(cm(1,1))),:) == 0);
                GridDivision(CurrentNC,APPosBinImage(round(cm(1,2)),round(cm(1,1))),dv_set_range)=CurrentFrame;
            elseif SelectionMode == 2
                ap_set_range = find(GridDivision(CurrentNC,:,DVPosBinImage(round(cm(1,2)),round(cm(1,1)))) == 0);
                GridDivision(CurrentNC,ap_set_range,DVPosBinImage(round(cm(1,2)),round(cm(1,1))))=CurrentFrame;
            elseif SelectionMode == 4
                schnitz_idx = SelectNucleus(schnitzcells, Ellipses, CurrentFrame, [cm(1,1), cm(1,2)]);
                schnitzAnaphaseFrames(schnitz_idx) = CurrentFrame;
            elseif SelectionMode == 5
                schnitz_idx = SelectNucleus(schnitzcells, Ellipses, CurrentFrame, [cm(1,1), cm(1,2)]);
                schnitzAnaphaseFrames(schnitz_idx) = 0;
            end
            
            
        end
        
        %Debug mode
    elseif (ct~=0)&(cc=='9')
        keyboard;
    end
end

NCTest = [];
for k=1:length(APboxes)
    NCTest(k) = GridDivision(CurrentNC, APboxes(k), DVboxes(k));
end
if min(NCTest) == 0
    incomplete_bins = find(NCTest == 0);
    disp('Finish selecting anaphase frames for the current nuclear cycle.')
    for k =1:length(incomplete_bins)
    disp(['AP bin: ', num2str(APboxes(incomplete_bins(k))), ', DV bin: ', num2str(DVboxes(incomplete_bins(k)))]);
    end
    prompt = ['Do you want to autofill remaining regions (y/n)? '];
    ID = input(prompt,'s');
    if ID == 'y' 
        assigned_rows = [];
        assigned_columns = [];
        associated_anaphase = [];
        for k=1:length(APboxes)
            if GridDivision(CurrentNC, APboxes(k), DVboxes(k)) > 0
              [binrows, bincolumns] = find((APPosBinImage ==  APboxes(k)) &  ...
                  (DVPosBinImage ==  DVboxes(k)));
              assigned_rows = [assigned_rows; binrows];
              assigned_columns = [assigned_columns; bincolumns];
              associated_anaphase = [associated_anaphase; GridDivision(CurrentNC, APboxes(k), DVboxes(k))*ones(size(binrows))];
            end
        end
        for l=1:length(incomplete_bins)
           idx = incomplete_bins(l);
           closest_anaphase_frames = [];
           [current_binrows, current_bincolumns] = find((APPosBinImage ==  APboxes(idx)) &  ...
                  (DVPosBinImage ==  DVboxes(idx)));
            for m=1:length(current_binrows)
                Distances=sqrt((assigned_rows-current_binrows(m)).^2+...
                    (assigned_columns-current_bincolumns(m)).^2);
                closest_anaphase_frames(m) = associated_anaphase(find(Distances == min(Distances), 1));
            end
            GridDivision(CurrentNC, APboxes(idx), DVboxes(idx)) = round(median(closest_anaphase_frames));
        end
    end
end
close(figureOverlay);
%% 
 

schnitz_ncs = [schnitzcells(:).cycle];
for nc=9:14
    if isempty(find([schnitzcells(:).cycle] == nc, 1))
        continue
    end
    for ap_idx=1:(length(APbinID)-1)
        APbinMin = APbinID(ap_idx);
        APbinMax = APbinID(ap_idx+1);
        for dv_idx=1:(length(DVbinID)-1)
            DVbinMin = DVbinID(dv_idx);
            DVbinMax = DVbinID(dv_idx+1);
            anaFrame = GridDivision(nc, ap_idx, dv_idx);
            if anaFrame > 0
                
                %[matchrows, matchcolumns] = find((APPosBinImage == ap_idx) & (DVPosBinImage == dv_idx));
                %boundary1 = boundary(matchcolumns, matchrows);
                
                ellipses_x = Ellipses{anaFrame}(:,1);
                ellipses_y = Ellipses{anaFrame}(:,2);
                
                %schnitzes_idx = Ellipses{anaFrame}(schnitz_match, 9);
                for k=1:size(Ellipses{anaFrame}, 1)
                    sc = Ellipses{anaFrame}(k, 9);
                    if schnitzAnaphaseFrames(sc) == 0
                        xpos = min([max([round(ellipses_x(k)), 1]), size(APPosImage, 2)]);
                        ypos = min([max([round(ellipses_y(k)), 1]), size(APPosImage, 1)]);
                        schnitzAP = APPosImage(ypos, xpos);
                        schnitzDV = DVPosImage(ypos, xpos);
                        %boundary2 = boundary([matchcolumns; round( Ellipses{anaFrame}(k, 1))],...
                        %    [matchrows; round( Ellipses{anaFrame}(k, 2))]);
                        %if isempty(find(length(matchrows)+1 == boundary2, 1))
                        if (schnitzAP >= APbinMin) & (schnitzAP <= APbinMax) & ...
                                (schnitzDV >= DVbinMin) & (schnitzDV <= DVbinMax)
                            %disp(['yes', num2str(k)])
                            schnitzAnaphaseFrames(sc) = anaFrame;
                            
                        end

                    end
                end
            end
            
        end
    end 
end

%% 
for i = 1:length(schnitzcells)
    if schnitzAnaphaseFrames(i) > 0
        schnitzcells(i).anaphaseFrame = schnitzAnaphaseFrames(i);
    end
end
%% 
% Some checks that things are working here
IDed_sc = [];
for i = 1:length(schnitzcells)
    if ~isempty(schnitzcells(i).anaphaseFrame) 
        IDed_sc(length(IDed_sc)+1) = i;
    end  
end
goodAF = zeros(1, length(IDed_sc));
goodFF = zeros(1, length(IDed_sc));
for j=1:length(IDed_sc)
    idx = IDed_sc(j);
    goodAF(j) = schnitzcells(idx).anaphaseFrame;
    goodFF(j) = schnitzcells(idx).frames(1);
end
    


%% 

adjust_idx = 0;
for i =1:length(schnitzcells)
    sc_idx = i + adjust_idx;
    if ~isempty(schnitzcells(sc_idx).anaphaseFrame)
        if schnitzcells(sc_idx).frames(1) ~= schnitzcells(sc_idx).anaphaseFrame
            anaFrame = schnitzcells(sc_idx).anaphaseFrame;
            approved_status = schnitzcells(sc_idx).Approved;
            schnitzcells=SeparateNuclearTraces(sc_idx,anaFrame,schnitzcells, FrameInfo, ncs);
            % First part of split 
            schnitzcells(sc_idx).cycle = schnitzcells(sc_idx).cycle -1;
            schnitzcells(sc_idx).anaphaseFrame = [];
            schnitzcells(sc_idx).inferredAnaphaseFrame = false;
            schnitzcells(sc_idx).Approved = 0;
            % Second part of split 
            schnitzcells(sc_idx+1).anaphaseFrame = anaFrame;
            schnitzcells(sc_idx+1).inferredAnaphaseFrame = false;
            schnitzcells(sc_idx+1).Approved = approved_status;
            anaTime = FrameInfo(anaFrame).Time/60;
            for j=1:length(schnitzcells(sc_idx+1).frames)
                schnitzcells(sc_idx+1).timeSinceAnaphase(j) = ...
                    FrameInfo(schnitzcells(sc_idx+1).frames(j)).Time/60-anaTime;
            end
            adjust_idx = adjust_idx + 1;
        end
    end
end
%% 
[schnitzcells, Ellipses] = correctSchnitzCellErrors(schnitzcells, Ellipses);

%% Infer anaphaseFrames for missing anaphases
schnitz_ncs = [schnitzcells(:).cycle];
hasAnaphaseFrame = zeros(1, length(schnitzcells));
for i = 1:length(schnitzcells)
    if ~isempty(schnitzcells(i).anaphaseFrame)
        hasAnaphaseFrame(i) = 1;
    end
end
for i=1:length(ncs)
    nc = i+8;
    if ncs(i) == 0
        continue
    end
    matchAnaFrameIdx = find((schnitz_ncs == nc) & hasAnaphaseFrame);
    matchNoAnaFrameIdx = find((schnitz_ncs == nc) & ~hasAnaphaseFrame);
    for j=1:length(matchNoAnaFrameIdx)
        numframes = length(schnitzcells(matchNoAnaFrameIdx(j)).frames);
        frameidx = min([20, numframes]);
        inferenceframe = schnitzcells(matchNoAnaFrameIdx(j)).frames(frameidx);
        framecell= schnitzcells(matchNoAnaFrameIdx(j)).cellno(frameidx);
        xpos = schnitzcells(matchNoAnaFrameIdx(j)).cenx(frameidx);
        ypos = schnitzcells(matchNoAnaFrameIdx(j)).ceny(frameidx);
        allschnitz_idx = Ellipses{inferenceframe}(:,9);
        goodschnitz_idx = [];
        for k=1:length(allschnitz_idx)
            if ismember(allschnitz_idx(k), matchAnaFrameIdx)
                goodschnitz_idx(length(goodschnitz_idx)+1) = k;
            end
        end
        if isempty(goodschnitz_idx)
            continue
        end
        good_xpos = Ellipses{inferenceframe}(goodschnitz_idx,1);
        good_ypos = Ellipses{inferenceframe}(goodschnitz_idx,2);
        Distances=sqrt((good_xpos-xpos).^2+(good_ypos-ypos).^2);
        inferenceidx = find(Distances == min(Distances));
        inferenceschnitz = Ellipses{inferenceframe}(goodschnitz_idx(inferenceidx), 9);
        anaphaseFrame = schnitzcells(inferenceschnitz).anaphaseFrame;
        schnitzcells(matchNoAnaFrameIdx(j)).anaphaseFrame = anaphaseFrame;
        schnitzcells(matchNoAnaFrameIdx(j)).inferredAnaphaseFrame = true;
        anaTime = FrameInfo(anaphaseFrame).Time/60;
        for l=1:length(schnitzcells(matchNoAnaFrameIdx(j)).frames)
            schnitzcells(matchNoAnaFrameIdx(j)).timeSinceAnaphase(l) = ...
                FrameInfo(schnitzcells(matchNoAnaFrameIdx(j)).frames(l)).Time/60-anaTime;
        end
    end
end


%% 
[schnitzcells, Ellipses] = correctSchnitzCellErrors(schnitzcells, Ellipses);

%% 
save([DropboxFolder,filesep,Prefix,filesep,'GridDivision.mat'],'GridDivision');
disp('GridDivision.mat saved.');
FilePrefix = [Prefix, '_'];
if whos(var2str(schnitzcells)).bytes < 2E9

    save([DropboxFolder, filesep, FilePrefix(1:end-1), filesep, FilePrefix(1:end-1), '_lin.mat'], 'schnitzcells', '-v6')
else
    save([DropboxFolder, filesep, FilePrefix(1:end-1), filesep, FilePrefix(1:end-1), '_lin.mat'], 'schnitzcells', '-v7.3', '-nocompression')
end
save2([DropboxFolder,filesep,Prefix,filesep,'Ellipses.mat'], Ellipses); 

