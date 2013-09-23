filename='D:\MS2LiveImaging\LiveImaging700\SnaEThsPMovie';

name = 'SnaEThsPMovie';

Prefix=name;

filesep='\';

nametosave='SnaEThsPMovie'

foldertosavein='D:\MS2LiveImaging\FISHPath\Data\RawData-2013-09-12'

OutputFolder=foldertosavein;

%%%%%%%%%%%%%%%% Creation of FrameInfo %%%%%%%%%%%%%%%%%%%

load(filename)

N=Datas.LSM_info.DimensionTime; % Number of time intervals

for i=1:N
    
FrameInfo(i).LinesPerFrame=Datas.LSM_info.DimensionX;
FrameInfo(i).PixelsPerLine=Datas.LSM_info.DimensionY;

FrameInfo(i).ZoomFactor=1;

FrameInfo(i).Rotation=0;

FrameInfo(i).ScanAmplitudeX=Datas.LSM_info.DisplayAspectX;
FrameInfo(i).ScanAmplitudeY=Datas.LSM_info.DisplayAspectY;

%%%%% for time string

n = datenum('30-December-1899');

TimeDays = n+Datas.LSM_info.AllData.EntryNumber35.Data +(Datas.LSM_info.TimeInterval/2 + Datas.LSM_info.Timeinfo.(['Posi', num2str(i)])-Datas.LSM_info.Timeinfo.Posi1)/60/60/24;

FrameInfo(i).TimeString=datestr(TimeDays,'dd/mm/yyyy HH:MM:SS');

%%%%%%%%%%%%%%%%

FrameInfo(i).XPosition=Datas.LSM_info.OriginX*10^6;
FrameInfo(i).YPosition=Datas.LSM_info.OriginY*10^6;
FrameInfo(i).ZPosition=Datas.LSM_info.OriginZ*10^6;

FrameInfo(i).NumberSlices=Datas.LSM_info.DimensionZ;

FrameInfo(i).ZStep=Datas.LSM_info.VoxelSizeZ*10^6;

FrameInfo(i).nc=13; % Dummy for the moment.

end

save([foldertosavein,'\FrameInfo', name, '.mat'],...
    'FrameInfo')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Saving of Tifs

for i=1:N

Im = LoadLsmToMat('D:\MS2LiveImaging\LiveImaging700\SnaEThsPMovie',i)

for j=1:Datas.LSM_info.DimensionZ
    
    II=Im.(['Slice', num2str(j)]);
    
    nametif = ([foldertosavein, '\', name,'_', num2str(floor(i/100)),num2str(floor((i-floor(i/100)*100)/10)),num2str(i-floor(i/100)*100-floor((i -floor(i/100)*100)/10)*10),'_z',num2str(floor(j/10)),num2str(floor(j)-10*floor(j/10)),'.tif'])
    
    imwrite(II(:,:,1),nametif,'tiff')
    
end
end



%%%%%%%%%%  Generating and saving tag file %%%%%%%%%%%%%


Output{1}=['id ',Prefix,'_'];
Output{2}='';
Output{3}='1';
Output{4}=['frames ',num2str(length(FrameInfo)),':1:',num2str(Datas.LSM_info.DimensionZ+2)];
Output{5}=['suffix ???_z??'];
Output{6}=['flat FF'];

%Create the TAG file for the FISH analysis code
fid = fopen([OutputFolder,filesep,Prefix,'.tag'], 'wt');

for i=1:length(Output)
    fprintf(fid, '%s \n', Output{i});
end

fclose(fid);

    

