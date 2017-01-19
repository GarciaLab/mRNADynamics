function ImportTr2d(Prefix,ObjectFile,TrackFile)

ObjectFile='tr2d_objects.csv';
TrackFile='tr2d_tracks.csv';

%This function grabs the data exported by tr2d and saves it into the
%corresponding Prefix in the format of Ellipses and lineages.

%To do: How does tr2d know where to save the files?

%Get folder and movie length information for Prefix
[SourcePath,FISHPath,DropboxFolder,MS2CodePath,PreProcPath]=...
    DetermineLocalFolders(Prefix);
load([DropboxFolder,filesep,Prefix,filesep,'FrameInfo.mat'])

%Load the tr2d files
Objects=csvread(ObjectFile,2,0);    %Read starting at the second row
Tracks=csvread(TrackFile,2,0);

%Create the Ellipses structure. The format is:
%(x, y, a, b, theta, maxcontourvalue, time, particle_id)

%The tr2d object format is (note that time and IDs start from 0):
%# t	 id	 area	 com_x	 com_y	 angle	 r1	 r2

%Note that we're mapping Area in tr2d to maxcontourvalue in Ellipses. We
%were not using the latter at all in the code.

Ellipses=cell(length(FrameInfo));
for i=1:length(FrameInfo)
    %Find all objects in this frame
    FrameFilter=(Objects(:,1)==i-1);
    %Copy the information
    Ellipses{i}=Objects(FrameFilter,[4,5,7,8,6,3]);
    Ellipses{i}(:,7)=i; %Frame
    Ellipses{i}(:,8)=Objects(FrameFilter,2)+1;
end

%Create the schnitzcells structure. The format is:
% P
% E
% D
% frames
% cenx
% ceny
% len
% cellno

%The tr2d format is
%# tracklet_id	 parent_tracklet_id	 child_tracklat_id1	 child_tracklat_id2	 (time	 object_id)...
%Note that time and ID are given as a pair of numbers. Also, the csv is
%padded with zeroes at the end of the list of frames and IDs.

for i=1:size(Tracks,1)
    %Find where the list of frames and IDs ends. I need to be careful with
    %to potentially problematic cases
    %1) There is no padding because this is the longest track
    %2) The track exists for only the first frame (denoted by 0 in tr2d)
    
    %Check that this track is not one of the longest ones. To do this,
    %we'll look for the first pair of [0,0] in the frame/index part of
    %Tracks
    
    if ~isempty(findstr(Tracks(i,7:end),[0,0]))
        MaxColumn=min(findstr(Tracks(i,7:end),[0,0]))+5;
        %If MaxColumn is an odd number, it means there was a 0 for the id
        %before the last pair of zeros that mark the end of the track
        if mod(MaxColumn,2)
            MaxColumn=MaxColumn+1;
        end
    else
        MaxColumn=size(Tracks,2);
    end
    
    %Actual track information:
    schnitzcells(i).frames=Tracks(i,5:2:MaxColumn)+1;
    schnitzcells(i).cellno=Tracks(i,6:2:MaxColumn)+1;
    schnitzcells(i).P=Tracks(i,2)+1;
    schnitzcells(i).D=Tracks(i,3)+1;
    schnitzcells(i).E=Tracks(i,4)+1;
    
    %Information coming from Ellipses. Our code should be modified to not
    %need to use duplicated information
    for j=1:length(schnitzcells(i).frames)
        EllipseToPull=...
            find(Ellipses{schnitzcells(i).frames(j)}(:,8)==...
            schnitzcells(i).cellno(j));
        
        schnitzcells(i).cenx(j)=...
            Ellipses{schnitzcells(i).frames(j)}(EllipseToPull,1);
        schnitzcells(i).ceny(j)=...
            Ellipses{schnitzcells(i).frames(j)}(EllipseToPull,2);
        %Calculate the length as the major axis length
        schnitzcells(i).len(j)=...
            max([Ellipses{schnitzcells(i).frames(j)}(EllipseToPull,3),...
                Ellipses{schnitzcells(i).frames(j)}(EllipseToPull,4)]);
    end
end


%Save the information
save([DropboxFolder,filesep,Prefix,filesep,'Ellipses.mat'],'Ellipses')
save([DropboxFolder,filesep,Prefix,filesep,Prefix,'_lin.mat'],...
    'schnitzcells')

