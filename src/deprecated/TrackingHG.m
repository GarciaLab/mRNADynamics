function TrackingHG

%% Information about about folders

[Dummy,XLS]=xlsread('ComputerFolders.xlsx');

%Find out which computer this is. That will determine the folder structure.
[ret, name] = system('hostname');  
if ret ~= 0,  
   if ispc  
      name = getenv('COMPUTERNAME');  
   else  
      name = getenv('HOSTNAME');  
   end  
end  
name = lower(name); 


%Find which computer we are dealing with:
ComputerColumn=find(strcmp(XLS(1,:),name(1:end-1)));

%Now load the corresponding folders
SourceRow=find(strcmp(XLS(:,1),'SourcePath'));
FISHRow=find(strcmp(XLS(:,1),'FISHPath'));
DropboxRow=find(strcmp(XLS(:,1),'DropboxFolder'));
SchnitzRow=find(strcmp(XLS(:,1),'SchnitzcellsFolder'));



%Assign the folders
SourcePath=XLS{SourceRow,ComputerColumn};
FISHPath=XLS{FISHRow,ComputerColumn};
DropboxFolder=XLS{DropboxRow,ComputerColumn};
SchnitzcellsFolder=XLS{SchnitzRow,ComputerColumn};


%% Initial tracking

%START segmentation and tracking :
Prefix = '2012-06-16-MCP(10)TM3-X1';

%Create the cell array with the names.
D=dir([FISHPath,'\Data\',Prefix,'\*His*.tif']);
for i=1:length(D)
    names{i}=[FISHPath,'\Data\',Prefix,filesep,D(i).name];
end

%Pull the mitosis information from the APDivision MAT file. I'll have to
%do this differently later.
load([DropboxFolder,filesep,Prefix,'\APDivision.mat'])
APDivision(APDivision==0)=nan;

MinAPDivision=min(APDivision')';
MinAPDivision=MinAPDivision(~isnan(MinAPDivision));
MaxAPDivision=max(APDivision')';
MaxAPDivision=MaxAPDivision(~isnan(MaxAPDivision));
clear APDivision

%Note that I'm adding a range of two frames frames before and after the
%determines mitosis
indMit=[MinAPDivision-2,MaxAPDivision+2];

%Embryo mask
ImageTemp=imread(names{1});
embryo_mask=true(size(ImageTemp));
clear ImageTemp



[nuclei, centers] = mainTracking(names,'indMitosis',indMit,'embryoMask', embryo_mask);
% names is a cell array containing the names of all frames in the movie in order.
% indMitosis is an nx2 array containing the first and last frame of mitosis in every row.
% embryoMask is the possible mask of the embryo. If no embryo edge is visible,
% true(size(an_image_from_the_movie)) can be given as input.


%% Convert the results to compatible structures and save them

%PUT circles on the nuclei

[Ellipses] = putCirclesOnNuclei(Prefix,centers,names,indMit);

%Convert nuclei structure into schnitzcell structure
[schnitzcells] = convertNucleiToSchnitzcells(nuclei); 

%Save the information
save([DropboxFolder,filesep,Prefix,'\Ellipses.mat'],'Ellipses')
save([DropboxFolder,filesep,Prefix,filesep,Prefix,'_lin.mat'],'schnitzcells')



%whos('-file',[DropboxFolder,filesep,Prefix,'\OldCode\Ellipses.mat'])
%whos('-file',[DropboxFolder,filesep,Prefix,'\OldCode\2012-06-16-MCP(10)TM3-X1_lin.mat'])


%% Extra

%Approve all centers to being with
for i=1:length(centers)
    approvedCenters{i}=ones(size(centers{i},1),1);
end
[schnitzcells] = convertNucleiToSchnitzcells(nuclei,approvedCenters); 
