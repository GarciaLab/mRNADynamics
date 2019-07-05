function  [DV_shift] = FindDVShift_full(Prefix)

%% Part 1: Read image data

[SourcePath, ~, DefaultDropboxFolder, DropboxFolder, ~, ~,...
~, ~] = DetermineAllLocalFolders(Prefix);

%Figure out what type of experiment we have. Note: we name the var "DateFromDateColumn" to avoid shadowing previously defined "Date" var.
[DateFromDateColumn, ExperimentType, ExperimentAxis, CoatProtein, StemLoop, APResolution,...
Channel1, Channel2, Objective, Power, DataFolder, DropboxFolderName, Comments,...
nc9, nc10, nc11, nc12, nc13, nc14, CF] = getExperimentDataFromMovieDatabase(Prefix, DefaultDropboxFolder)

% Read full embryo image (surf, mid)
midImage=imread([DropboxFolder,filesep,Prefix,filesep,'APDetection',filesep,'FullEmbryo.tif']);

load([DropboxFolder,filesep,Prefix,filesep,'APDetection.mat'])
load([DropboxFolder,filesep,Prefix,filesep,'Ellipses.mat'], 'Ellipses')

if size(midImage, 1) == 1024
    AreaThresh = 50;
    AreaMax = 250;
    %AreaThresh = 10;
    %AreaMax = 100;
elseif size(midImage, 1) == 2048
    AreaThresh = 100;
    AreaMax = 450;
else
    disp('full embryo resolution found not supported. talk to jake.')
end


FullEmbryoSurf=imread([DropboxFolder,filesep,Prefix,filesep,'DV',filesep,'surf_max.tif']);



%% Part 2: Label Image (Classified with Weka)
figure(1)
classifiedFullEmbryoImage = imread([DropboxFolder,filesep,Prefix,filesep,'DV',filesep,'Classified_image.tif']);
invertedClassifiedFullEmbryoImage = 1-classifiedFullEmbryoImage;    %Inverse image
imshow(invertedClassifiedFullEmbryoImage,[]);
ImLabel=bwlabel(invertedClassifiedFullEmbryoImage,4); %First, we need to give each regions of 1s a unique identity.
%imshow(ImLabel,[])

figure(2)
ImProps=regionprops(ImLabel,'Area'); %We need to calculate the area of each region found by bwlabel.
Areas=[ImProps.Area]; %Generate a list of all areas we have
edges = 0:5:AreaMax+20; %Plot histagram for area
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
for i=1:length(Areas) %Loop through each regions in ImLabel
    if (Areas(i)>AreaThresh) && (Areas(i)<AreaMax) %copy the region to NewImage
        SingleCellImage=ImLabel==i;         %Isolate the single cell
        NewImage=NewImage+SingleCellImage;         %Add it to the NewImage while preserving what I had before in NewImage
%         figure         %This is eye candy, you can comment it out
%         imshow(NewImage)
%         pause(0.5)
    end    
end

figure(3); %Plot new thresholded and filtered image
imshow(NewImage);
ImLabel2=bwlabel(NewImage); %relabel image using bwlabel
%imshow(ImLabel2,[])


%% Part 4: Calculate average fluorescence and position for each cell
ImFluo = FullEmbryoSurf; %Load surface image
%Finally, we want to calculate the mean fluorescence per
%pixel per cell for all the cells on this image. Dividing
%by area will make our results a little bit more
%insensitive to cell clumps.

ImProps2=regionprops(ImLabel2,'Area'); %Let's get a list of the new areas then
Areas2=[ImProps2.Area];
[x_pos, y_pos] = meshgrid(1:size(classifiedFullEmbryoImage,2), 1:size(classifiedFullEmbryoImage,1));
for i=1:length(Areas2) %For-loop to calculate the total fluorescence and average position of each cell
    ImMask=(ImLabel2==i);     %Generate the mask for the i-th cell
    SingleCellImFluo=immultiply(ImMask,ImFluo);     %Multiply the mask by the fluorescence image
    CellFluo(i)=sum(sum(SingleCellImFluo));     %Calculate and store the total fluorescence
    Im_x = ImMask.*x_pos;
    Im_y = ImMask.*y_pos;
    x_ave(i) = sum(sum(Im_x))/Areas2(i);
    y_ave(i) = sum(sum(Im_y))/Areas2(i);
end
CellFluoPerArea=CellFluo./Areas2; %Now divide the total cell fluorescence in CellFluo by the cell area

figure(4); %Plot the histogram
histogram(CellFluoPerArea);
xlabel('fluorescence per pixel');
ylabel('number of cells');

%% Part 3: Calculate Ellipse Position

%Get the information about the AP axis as well as the image shifts
%used for the stitching of the two halves of the embryo
load([DropboxFolder,filesep,Prefix,filesep,'APDetection.mat'])
load([DropboxFolder,filesep,Prefix,filesep,'Ellipses.mat'], 'Ellipses')

APLength = norm(coordPZoom - coordAZoom);
APLength_full = norm(coordP-coordA);
AP_Ratio=APLength/APLength_full;
APAngle = atan2((coordP(2)-coordA(2)),(coordP(1)-coordA(1)));

AP_max = 0;
AP_min = 1;

% Find AP position range for Ellipses
for i=1:length(Ellipses)
    for j=1:size(Ellipses{i},1)
        angleToAPAxis=atan2((Ellipses{i}(j,2)-coordAZoom(2)),(Ellipses{i}(j,1)-coordAZoom(1))); %Angle between the x-axis and the particle using the A position as a zero
        distToAPAxis=sqrt((coordAZoom(2)-Ellipses{i}(j,2)).^2+(coordAZoom(1)-Ellipses{i}(j,1)).^2); %Distance between the points and the A point
        APPositions=distToAPAxis.*cos(angleToAPAxis-APAngle);
        EllipsePos{i}(j)=APPositions/APLength;

        DVPositions=distToAPAxis.*sin(angleToAPAxis-APAngle);
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
close all
redLineFig = figure();
redLineAxes = axes(redLineFig);

imshow(midImage,[0,2]);
hold on

APx = [coordA(1),coordP(1)];
APy = [coordA(2),coordP(2)];
plot(APx,APy,'o-r','LineWidth',2);

A = [APx(1), APy(1)]; 
P = [APx(2), APy(2)];
AP = P - A; 
APMid = (P+A)*.5;
plot(APMid(1), APMid(2), 'o-y')

DV = [AP(2), -AP(1)];
% APAngle = atan2((coordP(2)-coordA(2)),(coordP(1)-coordA(1)));
% DVAngle = atan2(DV(1), DV(2));
% dvmag = norm(AP);
% DV2 = [DV(1)-dvmag*cos(DVAngle), DV(2)-dvmag*sin(DVAngle)]
% plot([APMid(1),DV2(1)],[APMid(2),DV2(2)],'o-r','LineWidth',2);
% 


%Find the equations of the lines
mAP=(APy(1)-APy(2))/(APx(1)-APx(2)); %slope of the ap axis
bAP=APy(1)-mAP*APx(1); %intercept of the ap axis

%Invert mAP to get tangent line that is perpendicular
mDV = -1/mAP;
dvx = APMid(1);
dvy = dvx*mDV + APMid(2);
plot([APMid(1), dvx], [APMid(2), dvy], 'o-y')
goodNucleus = 0;

for i=1:length(Areas2)
    b = y_ave(i)-mDV*x_ave(i);
    x_int = (b - bAP) / (mAP - mDV);
    y_int = mAP*x_int + bAP; 
    
    angleToAPAxis = atan2((y_ave(i)-coordA(2)),...
                    (x_ave(i)-coordA(1)));
    %DVpos(i) = sqrt((x_ave(i)-x_int)^2+(y_ave(i)-y_int)^2);
    distToAPAxis = norm(coordA - [x_ave(i), y_ave(i)]);
    
    APpos_temp = distToAPAxis.*cos(angleToAPAxis-APAngle)/APLength*AP_Ratio;
    APpos(i) = APpos_temp;
    if (APpos_temp<AP_max) && (APpos_temp>AP_min)
        APpos(goodNucleus+1) = APpos_temp;
        DVpos(goodNucleus+1) = distToAPAxis.*sin(angleToAPAxis-APAngle)*AP_Ratio;
        Cell_Fluo(goodNucleus+1) = CellFluoPerArea(i);
        goodNucleus = goodNucleus + 1;
    end

    %if (x_ave(i)-x_int)<0
    %    DVpos(i) = -DVpos(i);
    %end
%     plot([x_int,x_ave(i)],[y_int,y_ave(i)],'o-r','LineWidth',0.5);
%     pause(0.1);
end

figure(7);
plot3(x_ave,y_ave,CellFluoPerArea,'o');
%% Part 6: Fit Pattern to Gaussian- %Fit data using curve fitting toolbox

mygauss = fittype('a1*exp(-((x-b1)/c1)^2)',...
            'dependent',{'y'},'independent',{'x'},...
            'coefficients',{'a1','b1','c1'});
        
%If you change the number of parameters, make sure to update the
%lower and upper bounds to include ranges for added/removed
%parameters

options = fitoptions(mygauss);
options.Lower = [0 -Inf 100];
options.Upper = [Inf Inf Inf];
[temp_gauss, temp_gof] = fit(DVpos',Cell_Fluo',mygauss,options);

figure(8); %Plot fitted data
plot(temp_gauss,DVpos, Cell_Fluo);
xlabel('distance to AP axis (pixels)');
ylabel('intensity (au)');
legend('dorsal-venus (au)', 'Gaussian fit');


DV_shift = temp_gauss.b1;


end