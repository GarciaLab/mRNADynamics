function GenerateFakeNuclei(Prefix)

%When dealing with red proteins we can get the outlines of the nuclei.
%We'll use this information to generate a fake nuclear image that we can
%use to do tracking.

%Get the relevant
[SourcePath,FISHPath,DropboxFolder,MS2CodePath,SchnitzcellsFolder]=...
    DetermineLocalFolders(Prefix);

%Load the information about the experiment
load([DropboxFolder,filesep,Prefix,filesep,'FrameInfo.mat'])


%This assumes that channel 2 is the one with the red protein showing the
%nuclear outlines.

%Let's start by using the center half of the frames only. We'll do a
%maximum projection and then we'll invert the image

StackCenter=round((FrameInfo(1).NumberSlices-1)/2);
StackRange=[StackCenter-round(StackCenter/2):StackCenter+round(StackCenter/2)];

h=waitbar(0,'Generating fake nuclear images');
for i=1:length(FrameInfo)
    waitbar(i/length(FrameInfo),h)
    Images=[];
    StackCounter=1;
    for j=StackRange
        Images(:,:,StackCounter)=imread([FISHPath,filesep,'Data',filesep,...
            Prefix,filesep,Prefix,'_',iIndex(i,3),'_z',iIndex(j,2),'_ch02.tif']);
        StackCounter=StackCounter+1;
    end
    MaxImage=max(Images,[],3);
    InverseImage=1-mat2gray(MaxImage);
    
%     DiskFilter=strel('disk',2);
%     FilterImage=imclose(InverseImage,DiskFilter);
%     
%     imshow(imadjust(FilterImage))
    imwrite(imadjust(InverseImage),[FISHPath,filesep,'Data',filesep,...
            Prefix,filesep,Prefix,'_',iIndex(i,3),'_FakeHis.tif'])
end
close(h)
