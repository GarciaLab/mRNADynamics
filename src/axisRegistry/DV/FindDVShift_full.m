function  [DV_shift] = FindDVShift_full(Prefix, varargin)

%% Part 1: Read image data

[SourcePath, ~, DefaultDropboxFolder, DropboxFolder, ~, ~,...
~, ~] = DetermineAllLocalFolders(Prefix);

%Figure out what type of experiment we have. Note: we name the var "DateFromDateColumn" to avoid shadowing previously defined "Date" var.
[DateFromDateColumn, ExperimentType, ExperimentAxis, CoatProtein, StemLoop, APResolution,...
Channel1, Channel2, Objective, Power, DataFolder, DropboxFolderName, Comments,...
nc9, nc10, nc11, nc12, nc13, nc14, CF] = getExperimentDataFromMovieDatabase(Prefix, DefaultDropboxFolder)

%Find out the date it was taken
Dashes=findstr(Prefix,'-');
Date=Prefix(1:Dashes(3)-1);
EmbryoName=Prefix(Dashes(3)+1:end);

%Find Zoom Factor
    FileMode='LIFExport';
    %Find the zoomed movie pixel size
    D=dir([SourcePath,filesep,Date,filesep,EmbryoName,filesep,'*.',FileMode(1:3)]);

    %Load only the metadata from the zoomed images
    MetaReader=bfGetReader([SourcePath,filesep,Date,filesep,EmbryoName,filesep,D(end).name]);
    MetaZoom=MetaReader.getMetadataStore();
    try
        PixelSizeZoom=str2double(MetaZoom.getPixelsPhysicalSizeX(0).value);
    catch
        PixelSizeZoom=str2double(MetaZoom.getPixelsPhysicalSizeX(0));
    end

    %Find the full embryo pixel size and load the image
    D=dir([SourcePath,filesep,Date,filesep,EmbryoName,filesep,'FullEmbryo',filesep,'*surf*.',FileMode(1:3)]);

    ImageTemp=bfopen([SourcePath,filesep,Date,filesep,EmbryoName,filesep,'FullEmbryo',filesep,D(end).name]);
    MetaFullEmbryo= ImageTemp{:, 4};
    PixelSizeFullEmbryo=str2double(MetaFullEmbryo.getPixelsPhysicalSizeX(0) );
    try
        PixelSizeFullEmbryo=str2double(MetaFullEmbryo.getPixelsPhysicalSizeX(0).value);
    catch
        PixelSizeFullEmbryo=str2double(MetaFullEmbryo.getPixelsPhysicalSizeX(0));
    end

    %Zoom factor
    MovieZoom=PixelSizeFullEmbryo(1)/PixelSizeZoom(1);
    SurfZoom=1;     %We'll call the zoom of the full embryo image 1


% Read full embryo image (surf, mid)
FullEmbryo=imread([DropboxFolder,filesep,Prefix,filesep,'APDetection',filesep,'FullEmbryo.tif']);

if size(FullEmbryo, 1) == 1024
    AreaThresh = 20;
    AreaMax = 100;
elseif size(FullEmbryo, 2) == 2048
    AreaThresh = 100;
    AreaMax = 450;
else
    disp('full embryo resolution found not supported. talk to jake.')
end


FullEmbryoSurf=imread([DropboxFolder,filesep,Prefix,filesep,'DV',filesep,'surf_max.tif']);



%% Part 2: Label Image (Classified with Weka)
%I = imread('./test2/binary.tif');
classifiedFullEmbryoImage = imread([DropboxFolder,filesep,Prefix,filesep,'DV',filesep,'Classified_image.tif']);
invertedClassifiedFullEmbryoImage = 1-classifiedFullEmbryoImage;    %Inverse image

figure(1)
imshow(invertedClassifiedFullEmbryoImage,[]);

%First, we need to give each regions of 1s a unique identity.
ImLabel=bwlabel(invertedClassifiedFullEmbryoImage,4);
%imshow(ImLabel,[])

%We need to calculate the area of each region found by bwlabel.
%To get the area, we'll use regionprops
ImProps=regionprops(ImLabel,'Area');

%Generate a list of all areas we have
Areas=[ImProps.Area];
%Plot histagram for area
%edges = [0:5:100];
edges = 0:5:AreaMax+20;
figure(2)
histogram(Areas,edges);
xlabel('area (pixels)')
ylabel('number of regions')


%% Part 3: Quality Control

%Note that most of our regions have one or two pixels.
%We want to get rid of them by filtering for areas
%that are larger than 20 and smaller than 100 pixels.


%We are going to extract the regions with an area larger
%than AreaThresh. To do this, we'll copy those selected
%regions onto a blank image.
NewImage=zeros(size(invertedClassifiedFullEmbryoImage));
                                              
%Loop through each regions in ImLabel
for i=1:length(Areas)
    if (Areas(i)>AreaThresh) && (Areas(i)<AreaMax) %copy the region to NewImage
        %Isolate the single cell
        SingleCellImage=ImLabel==i;
        %Add it to the NewImage while preserving
        %what I had before in NewImage
        NewImage=NewImage+SingleCellImage;
        %This is eye candy, you can comment it out
%         figure(1)
%         imshow(NewImage)
%         pause(0.5)
    end    
end

%Plot new thresholded and filtered image
figure(3);
imshow(NewImage);
ImLabel2=bwlabel(NewImage); %relabel image using bwlabel
%imshow(ImLabel2,[])


%% Part 4: Calculate average fluorescence and position for each cell
%Load surface image
ImFluo = FullEmbryoSurf;
%ImFluo = imread('./test2/surf.tif');

%Finally, we want to calculate the mean fluorescence per
%pixel per cell for all the cells on this image. Dividing
%by area will make our results a little bit more
%insensitive to cell clumps.

%Let's get a list of the new areas then
ImProps2=regionprops(ImLabel2,'Area');
Areas2=[ImProps2.Area];

[x_pos, y_pos] = meshgrid(1:size(classifiedFullEmbryoImage,2), 1:size(classifiedFullEmbryoImage,1));

%For-loop to calculate the total fluorescence and average position of each cell
parfor i=1:length(Areas2)
    %Generate the mask for the i-th cell
    ImMask=(ImLabel2==i);
    %Multiply the mask by the fluorescence image
    SingleCellImFluo=immultiply(ImMask,ImFluo);
    %Calculate and store the total fluorescence
    CellFluo(i)=sum(sum(SingleCellImFluo));
    Im_x = immultiply(ImMask,x_pos);
    Im_y = immultiply(ImMask,y_pos);
    x_ave(i) = sum(sum(Im_x))/Areas2(i);
    y_ave(i) = sum(sum(Im_y))/Areas2(i);
end

%Now divide the total cell fluorescence in CellFluo
%by the cell area
CellFluoPerArea=CellFluo./Areas2;

%Plot the histogram
figure(4);
histogram(CellFluoPerArea);
xlabel('fluorescence per pixel');
ylabel('number of cells');

%% Part 3: Calculate Ellipse Position

%Get the information about the AP axis as well as the image shifts
%used for the stitching of the two halves of the embryo
load([DropboxFolder,filesep,Prefix,filesep,'APDetection.mat'])
load([DropboxFolder,filesep,Prefix,filesep,'Ellipses.mat'], 'Ellipses')

APLength=sqrt((coordPZoom(2)-coordAZoom(2))^2+(coordPZoom(1)-coordAZoom(1))^2);
APLength_full=sqrt((coordP(2)-coordA(2))^2+(coordP(1)-coordA(1))^2);
AP_Ratio=APLength/APLength_full;
APAngle = atan2((coordP(2)-coordA(2)),(coordP(1)-coordA(1)));

AP_max = 0;
AP_min = 1;

% Find AP position range for Ellipses
for i=1:length(Ellipses)
    for j=1:size(Ellipses{i},1)

        %Angle between the x-axis and the particle using the A position as a
        %zero

        Angles=atan2((Ellipses{i}(j,2)-coordAZoom(2)),(Ellipses{i}(j,1)-coordAZoom(1)));

        %Distance between the points and the A point
        Distances=sqrt((coordAZoom(2)-Ellipses{i}(j,2)).^2+(coordAZoom(1)-Ellipses{i}(j,1)).^2);
        APPositions=Distances.*cos(Angles-APAngle);
        EllipsePos{i}(j)=APPositions/APLength;

        DVPositions=Distances.*sin(Angles-APAngle);
        EllipsePos_DV{i}(j)=DVPositions;
        
        if EllipsePos{i}(j) < AP_min
            AP_min = EllipsePos{i}(j);
        end
        if EllipsePos{i}(j) > AP_max
            AP_max = EllipsePos{i}(j);
        end
        
    end
end

AP_max = AP_max+0.1;
AP_min = AP_min-0.1;

%% Part 4: Calculate DV position

ImMid = FullEmbryo;

%{
AP=ginput(2); %Click once on the anterior most tip of the embryo and once on
%the posterior most tip of the embryo and then press enter. The two points
%will be stored as coordinates in a matrix called AP.

%Let's define two vectors for our clicks. The first will be our x-values.
APx=[AP(1,1),AP(2,1)];
%The second vector will be for our y-values
APy=[AP(1,2),AP(2,2)];
%}


APx = [coordA(1),coordP(1)];
APy = [coordA(2),coordP(2)];

imshow(ImMid,[0,2]);
hold on
plot(APx,APy,'o-r','LineWidth',4);

%Find the equations of the lines
mAP=(APy(1)-APy(2))/(APx(1)-APx(2));
bAP=APy(1)-mAP*APx(1);

%Invert mAP to get tangent line that is perpendicular
mAP_i = -1/mAP;

num = 0;

for i=1:length(Areas2)
    b = y_ave(i)-mAP_i*x_ave(i);
    x_int = (b-bAP)/(mAP-mAP_i);
    y_int = mAP*x_int+bAP;
    
    Angles = atan2((y_ave(i)-coordA(2)),...
                    (x_ave(i)-coordA(1)));
    %DVpos(i) = sqrt((x_ave(i)-x_int)^2+(y_ave(i)-y_int)^2);
    Distances = sqrt((coordA(2)-y_ave(i)).^2+(coordA(1)-x_ave(i)).^2);
    
    APpos_temp = Distances.*cos(Angles-APAngle)/APLength*AP_Ratio;
    APpos(i) = APpos_temp;
    if (APpos_temp<AP_max) && (APpos_temp>AP_min)
        num = num + 1;
        APpos(num+1) = APpos_temp;
        DVpos(num+1) = Distances.*sin(Angles-APAngle)*AP_Ratio;
        Cell_Fluo(num+1) = CellFluoPerArea(i);
    end

    %if (x_ave(i)-x_int)<0
    %    DVpos(i) = -DVpos(i);
    %end
    %plot([x_int,x_ave(i)],[y_int,y_ave(i)],'o-r','LineWidth',0.5);
    %pause(0.1);
end

figure(5);
plot3(x_ave,y_ave,CellFluoPerArea,'o');

%% Part 6: Fit Pattern to Gaussian

%Fit data using curve fitting toolbox
mygauss = fittype('a1*exp(-((x-b1)/c1)^2)',...
            'dependent',{'y'},'independent',{'x'},...
            'coefficients',{'a1','b1','c1'});
        
%If you change the number of parameters, make sure to update the
%lower and upper bounds to include ranges for added/removed
%parameters

options = fitoptions(mygauss);
%options.Lower = [0 -Inf 100];
%options.Upper = [Inf Inf 500];
[temp_gauss, temp_gof] = fit(DVpos',Cell_Fluo',mygauss,options);

%Plot fitted data
figure(6);
plot(temp_gauss,DVpos, Cell_Fluo);

%DV_shift = temp_gauss.b1*MovieZoom;
DV_shift = temp_gauss.b1;

end
