
function LabelNucsCore=SegmentNucleiLive(folder,DivisionTimes,FilterRadius,FilterSize)
  
 load([folder,'MaxNuclei.mat'],'MaxNuclei');
    
 TotalTime=length(fieldnames(MaxNuclei));  

if exist([folder,'LabelNucsCore.mat'])>0

     load([folder,'LabelNucsCore.mat']);
end
 
 
%%% Set up time dependent radius of filtering.

if exist('Radi')
else
    if isempty(DivisionTimes)
        
        SizeG=FilterSize;
        Radi=FilterRadius;
        
    else
        
        SizeG=FilterSize(end)*ones(1,TotalTime);
        Radi=FilterRadius(end)*ones(1,TotalTime);
        
        
        for i=1:length(DivisionTimes)
            
            SizeG(1:DivisionTimes(end-i+1))=FilterSize(end-i);
            Radi(1:DivisionTimes(end-i+1))=FilterRadius(end-i);
            
        end
    end
end

%%%%% Major filtering step

        for i=1:TotalTime % i is the global time index
            
            
            I = MaxNuclei.(['Time', num2str(i)]); % Image to use to filter
            
            NucCourseFiltT=segmentnucleiCoreLive(I,SizeG(i),Radi(i),0); % Find cores of nuclei, maybe need to adjust size with time
    
            NucCourseFiltT = bwmorph(NucCourseFiltT,'thicken',1);  % Add thickness                       
            
            LabelNucsCore.(['Time', num2str(i)]).Image=bwlabel(NucCourseFiltT); % Label the obejcts
            
            LabelNucsCore.(['Time', num2str(i)]).RegionProps=regionprops(LabelNucsCore.(['Time', num2str(i)]).Image,'PixelIdxList','Centroid'); % Linear indices
            
            
        end

%User input on splitting and joining nuclei  

f=figure;

imshow(label2rgbBackdropLive((LabelNucsCore.(['Time', num2str(1)]).Image>0)+1,[1,0,0;0,1,0],[1,1,1],2*imadjust(MaxNuclei.(['Time', num2str(1)]))), 'Border','tight','InitialMagnification',100)
        
   but = 1;
   
   textholder=[];
   
   Nlayer=1;
   
    while but ~= 3 % User input keep loop going until user presses right mouse
        
        figure(f)
        [xi,yi,but] = ginput(f);
        
        if but==30 % Up arrow
            
           Nlayer=Nlayer+1;
            
           if Nlayer==TotalTime+1; % Start back at layer 1 once go above 3
               Nlayer=1;
           else
           end
            
           clf(f)
           
           disp(['Image time', num2str(Nlayer)])
           
           imshow(label2rgbBackdropLive((LabelNucsCore.(['Time', num2str(Nlayer)]).Image>0)+1,[1,0,0;0,1,0],[1,1,1],2*imadjust(MaxNuclei.(['Time', num2str(Nlayer)]))), 'Border','tight','InitialMagnification',100), 
            

        elseif but==31 % Down arrow
            
           Nlayer=Nlayer-1;
            
           if Nlayer==0;
                Nlayer=TotalTime;
           else
           end
            
           clf(f)
           
           disp(['Image time', num2str(Nlayer)])

           imshow(label2rgbBackdropLive((LabelNucsCore.(['Time', num2str(Nlayer)]).Image>0)+1,[1,0,0;0,1,0],[1,1,1],2*imadjust(MaxNuclei.(['Time', num2str(Nlayer)]))), 'Border','tight','InitialMagnification',100), 

        elseif but==114 % r is pressed join nuclei that are selected
        
            textholder='r';
            
        elseif but == 98 % b is pressed split nuclei
            
            textholder='b';
            
         elseif but == 115 % s is pressed split nuclei
            
            textholder='s';
            
         elseif but == 97 % s is pressed split nuclei
            
            textholder='a';            
         
        elseif but == 105 % s is pressed split nuclei
            
            textholder='i';            
            
         elseif but==1 % Do something if left mouse button pressed
              
             if strcmp(textholder,'r')
                 
    Lab=LabelNucsCore.(['Time', num2str(Nlayer)]).Image;                 
    Bw=LabelNucsCore.(['Time', num2str(Nlayer)]).Image>0;
    BwDil=bwmorph(Bw,'dilate',1);
    BWErode=bwmorph(BwDil,'erode',1);
    LabBWErode=bwlabeln(BWErode);


    xy = []; % initially the list of points is empty
    NucleiIndex=[]; % initialize nuclei index

    NucleiIndex=sub2ind(size(Bw),round(yi),round(xi));            

    BWSelect = ismember(LabBWErode,LabBWErode(NucleiIndex));

    AllNuc=nonzeros(unique(Lab(:)));

    NucToGetRid=nonzeros(unique(Lab(find(BWSelect))));

    NucToKeep=setdiff(AllNuc,NucToGetRid);

    BwKeep=ismember(Lab,NucToKeep);

    BwFullKeep=BwKeep+BWSelect;

    imshow(BwFullKeep,'Border','tight','InitialMagnification',100), 

    LabelMatrix = bwlabeln(BwFullKeep);

    LabelNucsCore.(['Time', num2str(Nlayer)]).Image=LabelMatrix;

    RegionLabelMatrix=regionprops(LabelMatrix,'PixelIdxList','Centroid');

    LabelNucsCore.(['Time', num2str(Nlayer)]).RegionProps=RegionLabelMatrix;

    
    elseif strcmp(textholder,'b') % For removing nuclei


    Lab=LabelNucsCore.(['Time', num2str(Nlayer)]).Image;                 
    Bw=LabelNucsCore.(['Time', num2str(Nlayer)]).Image>0;

    xy = []; % initially the list of points is empty
    NucleiIndex=[]; % initialize nuclei index

    NucleiIndex=sub2ind(size(Bw),round(yi),round(xi));            

    BWSelect = ismember(Lab,Lab(NucleiIndex));

    AllNuc=nonzeros(unique(Lab(:)));

    NucToGetRid=Lab(NucleiIndex);

    NucToKeep=setdiff(AllNuc,NucToGetRid);

    BwKeep=ismember(Lab,NucToKeep);

    imshow(BwKeep,'Border','tight','InitialMagnification',100), 

    LabelMatrix = bwlabeln(BwKeep);

    LabelNucsCore.(['Time', num2str(Nlayer)]).Image=LabelMatrix;

    RegionLabelMatrix=regionprops(LabelMatrix,'PixelIdxList','Centroid');

    LabelNucsCore.(['Time', num2str(Nlayer)]).RegionProps=RegionLabelMatrix;                
                
   
           elseif strcmp(textholder,'a') %%%%% For changing radius of filter
 

              str=input(['Input new radius, old is ', num2str(Radi(Nlayer)), ' :'],'s');
                    
              NewRad=str2num(str);
              
              Radi(Nlayer)=NewRad;
               
             I = MaxNuclei.(['Time', num2str(Nlayer)]); % Image to use to filter
             
             NucCourseFiltT=segmentnucleiCoreLive(I,SizeG(Nlayer),Radi(Nlayer),0); % Find cores of nuclei, maybe need to adjust size with time
    
            NucCourseFiltT = bwmorph(NucCourseFiltT,'thicken',1);  % Add thickness                       
            
            LabelNucsCore.(['Time', num2str(Nlayer)]).Image=bwlabel(NucCourseFiltT); % Label the obejcts
            
            LabelNucsCore.(['Time', num2str(Nlayer)]).RegionProps=regionprops(LabelNucsCore.(['Time', num2str(Nlayer)]).Image,'PixelIdxList','Centroid'); % Linear indices
            
               
            imshow(label2rgbBackdropLive((LabelNucsCore.(['Time', num2str(Nlayer)]).Image>0)+1,[1,0,0;0,1,0],[1,1,1],2*imadjust(MaxNuclei.(['Time', num2str(Nlayer)]))), 'Border','tight','InitialMagnification',100),                
                    
           elseif strcmp(textholder,'i') %%%%% For changing size of filter
               
                 str=input(['Input new Size, old is ', num2str(SizeG(Nlayer)), ' :'],'s');
                    
              NewSize=str2num(str);
              
              SizeG(Nlayer)=NewSize;
               
             I = MaxNuclei.(['Time', num2str(Nlayer)]); % Image to use to filter
             
             NucCourseFiltT=segmentnucleiCoreLive(I,SizeG(Nlayer),Radi(Nlayer),0); % Find cores of nuclei, maybe need to adjust size with time
    
            NucCourseFiltT = bwmorph(NucCourseFiltT,'thicken',1);  % Add thickness                       
            
            LabelNucsCore.(['Time', num2str(Nlayer)]).Image=bwlabel(NucCourseFiltT); % Label the obejcts
            
            LabelNucsCore.(['Time', num2str(Nlayer)]).RegionProps=regionprops(LabelNucsCore.(['Time', num2str(Nlayer)]).Image,'PixelIdxList','Centroid'); % Linear indices
            
               
            imshow(label2rgbBackdropLive((LabelNucsCore.(['Time', num2str(Nlayer)]).Image>0)+1,[1,0,0;0,1,0],[1,1,1],2*imadjust(MaxNuclei.(['Time', num2str(Nlayer)]))), 'Border','tight','InitialMagnification',100),                

            end
        end
    end
    

%%%% Tracking the nuclei based on overlap alone

NumNucs = max(LabelNucsCore.Time1.Image(:));

LabelNucsCore.Time1.ImageZ=LabelNucsCore.Time1.Image;

for i=1:TotalTime-1
    
    NumNucsi = max(LabelNucsCore.(['Time', num2str(i+1)]).Image(:));
    [lenn,widd]=size(LabelNucsCore.(['Time', num2str(i+1)]).Image);
    
    LabelNucsCore.(['Time', num2str(i+1)]).ImageZ=zeros(lenn,widd);
    
    for j=1:NumNucsi

        IndexOfjthNuc = LabelNucsCore.(['Time', num2str(i+1)]).RegionProps(j).PixelIdxList;
        MostCommonIndex = mode(nonzeros(LabelNucsCore.(['Time', num2str(i)]).ImageZ(IndexOfjthNuc)));
      
        if MostCommonIndex>0
        
        LabelNucsCore.(['Time', num2str(i+1)]).ImageZ(IndexOfjthNuc) = MostCommonIndex;
        
        else
            
        LabelNucsCore.(['Time', num2str(i+1)]).ImageZ(IndexOfjthNuc) = NumNucs+1;
      
        NumNucs=NumNucs+1;
        
        end
    end
end


%%%% Check for doubling of indices
 NumNucs=[];

for i=1:TotalTime
    NumNucs=[unique(LabelNucsCore.(['Time', num2str(i)]).ImageZ(:)); NumNucs];
end

NumNucs=max(NumNucs)

for i=1:TotalTime-1
    
    NumNucsi = max(LabelNucsCore.(['Time', num2str(i+1)]).Image(:)); % Number of nuclei in the i+1 th time
    
    [lengthh,widthh]=size(LabelNucsCore.Time1.Image);
    
    LabelNucsCore.(['Time', num2str(i+1)]).ImageZ=zeros(lengthh,widthh); % Start with empty matrix to populate
    
    NewIndex=zeros(NumNucsi,1);
    
    for j=1:NumNucsi          % Run loop over total number of nuclei in i+1 th time

        
        IndexOfjthNuc = LabelNucsCore.(['Time', num2str(i+1)]).RegionProps(j).PixelIdxList; % Linear index of jth nuc
        MostCommonIndex = mode(nonzeros(LabelNucsCore.(['Time', num2str(i)]).ImageZ(IndexOfjthNuc))); % Most common index in jth nuclei in i th time
        
      
        if MostCommonIndex>0 % If there was a nucleus in the i'th layer, same one in the i+1'th layer
        
        LabelNucsCore.(['Time', num2str(i+1)]).ImageZ(IndexOfjthNuc) = MostCommonIndex;
        
        NewIndex(j)=MostCommonIndex;
        
        else                 % If no nucleus in the i'th layer make new label
            
        LabelNucsCore.(['Time', num2str(i+1)]).ImageZ(IndexOfjthNuc) = NumNucs+1;
        
        NewIndex(j)=NumNucs+1;
        
        NumNucs=NumNucs+1;
        
        end
    end
        
        
        %%%%% See if there is a doubling up of indexes
        
        SortNewIndex=sort(NewIndex);
        
        IndexesDup = find((SortNewIndex-circshift(SortNewIndex,1))==0);
        
        ActualIndex=SortNewIndex(IndexesDup);
        
        
        
        if isempty(IndexesDup)
        else
            for ii=1:length(ActualIndex)
                
                Indicess = find(NewIndex==ActualIndex(ii));
                
                for iii=1:length(Indicess)
                    
                    LabelNucsCore.(['Time', num2str(i+1)]).ImageZ(LabelNucsCore.(['Time', num2str(i+1)]).RegionProps(Indicess(iii)).PixelIdxList) = NumNucs+1;
                    
                    NumNucs=NumNucs+1;
                end
            end
            
        end
  
end

%%%%
    
 save([folder,'LabelNucsCore.mat'],'LabelNucsCore','Radi','SizeG');    
 
        
