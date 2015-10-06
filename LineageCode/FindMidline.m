% We find the midline of a frame, if it exists, by analysing all the nuclei in
% nuclear cycle 14. We do this by selecting the outliers of neighbourhood
% distances. The neighbourhood distance is defined as the vertical distance between a
% nucleus and its closest nucleus vertically

function typosns=FindMidline(PreProcPath, prefix, Ellipses, cf, dimx, lrep)
% This is the same as find midline but without graphics

%% Parameters - please redefine as needed
close all;
% We ignore the particles in the margins
Margins=0.05;
% We will be identifying the midline by looking at nuclei in a particular x
% bin. The following is the size of the AP bin
XBinSize=20;
ConfidenceInterval=2; % Confidence interval for identifying outliers (suggested is 2)
% filename='outliers4.avi'; % Filename to save the movie of the midline
LinesPerFrame=dimx;% Height of frames
NumberOfFrames=length(Ellipses);

%% Finding Neighbours
% Start a movie to watch where the neighbours occur
% mov=VideoWriter(filename);
% mov.FrameRate = 5;
% open(mov);

h=waitbar(0,'Finding Midline: Finding neighbours');

% We find the maximum dimension that the matrix of x-distances and
% y-distances can have to avoid resizing matrices every time
MaximumNumberOfNuclei=0;
for i=1:NumberOfFrames
    if length(Ellipses{i,1}(:,1))>MaximumNumberOfNuclei
        MaximumNumberOfNuclei=length(Ellipses{i,1}(:,1));
    end
end
xdist=cell(cf,MaximumNumberOfNuclei);
ydist=cell(cf,MaximumNumberOfNuclei);
dists1=zeros(cf,MaximumNumberOfNuclei);
dists2=zeros(cf,MaximumNumberOfNuclei);
for i =cf:NumberOfFrames
    waitbar((i-cf)/(NumberOfFrames-cf));
    x=Ellipses{i,1}(:,1); % x distance array for the frame
    y=Ellipses{i,1}(:,2); % y distance array for the frame
    
    % Choose a nucleus
    for j=1:length(y)
        %         I=imread(strcat('~/Documents/tgLab/LivemRNA/Data/PreProcessedData/2015-02-25-SogSha-hz/2014-12-30-SogSha-hz-His_',num2str(i),'.tif'));
        %         circles=[Ellipses{i,1}(j,1) Ellipses{i,1}(j,2) 15];
        %         I=insertShape(I, 'circle', circles,'Color','w');
        
        % Choose the nuclei within selected margins
        if(x(j)>Margins*LinesPerFrame && x(j)<(1-Margins)*LinesPerFrame)% & y(j)>margin*ht & y(j)<(1-margin)*ht)
            
            % Compute arrays of distances to other nuclei in frame
            xdist{i,j}=x-Ellipses{i,1}(j,1);
            ydist{i,j}=y-Ellipses{i,1}(j,2);
            ydist{i,j}(j)=inf; % Sets the distance to itself to infinity
            
            % slices is the array of indices of neighbours that are within a
            % selected x slice
            slices=find(abs(xdist{i,j})<XBinSize);
            
            % Find the closest on the y-axis within the x-distance margins
            xslicecheck=true;
            infcheck=true;
            while(xslicecheck || infcheck)
                [miny,indy]=min(abs(ydist{i,j}));
                xslicecheck=abs(xdist{i,j}(indy))>XBinSize; % Force loop again if not in xslice
                infcheck = sum(~isinf(ydist{i,j}))==0;
                if xslicecheck
                    ydist{i,j}(indy)=inf;
                end
                if infcheck
                    ydist{i,j}(j)=0;
                end
            end
            
            dists1(i,j) = xdist{i,j}(indy)^2+ydist{i,j}(indy)^2;
            %             circles=[Ellipses{i,1}(indy,1) Ellipses{i,1}(indy,2) 15];
            %             I=insertShape(I, 'circle', circles,'Color','r');
            
            % Find the 2nd closest
            prev=ydist{i,j}(indy);
            ydist{i,j}(indy)=inf;
            
            % Keeps checking if the 2nd closest is in the opposite direction as previous point
            % Find the closest on the y-axis within the x-distance margins
            xslicecheck=true;
            oppcheck=true;
            %infcheck=true;
            while(xslicecheck || oppcheck)
                [miny,indy]=min(abs(ydist{i,j}));
                xslicecheck=abs(xdist{i,j}(indy))>XBinSize; % Force loop again if not in xslice
                oppcheck=prev*ydist{i,j}(indy)>0;
                infcheck = sum(~isinf(ydist{i,j}))==0;
                if (xslicecheck || oppcheck)
                    ydist{i,j}(indy)=inf;
                end
                if infcheck
                    ydist{i,j}(j)=0;
                end
            end
            xslicecheck=true;
            oppcheck=true;
            infcheck=true;
            
            if(prev*ydist{i,j}(indy)<0)
                dists2(i,j) = xdist{i,j}(indy)^2+ydist{i,j}(indy)^2;
            else
                dists2(i,j)=min([y(j)^2,(512-y(j))^2]);
            end
            %             circles=[Ellipses{i,1}(indy,1) Ellipses{i,1}(indy,2) 15];
            %             I=insertShape(I, 'circle', circles,'Color','g');
        end
        %         imshow(I)
    end
    %     writeVideo(mov,getframe(gca));
end
close(h)
% close(mov)

%% Data Analysis
% Now, combine the data and do the data analysis
dists=[dists1 dists2];
valid_dists=dists(dists~=0 & ~isnan(dists) & ~isinf(dists));
stdev=std(valid_dists);
avg=mean(valid_dists);
outliers=zeros(size(dists));
density=zeros(size(dists,1));
out_nucleus=zeros(sum(dists(i,j)~=inf && dists(i,j)>(avg+ConfidenceInterval*stdev)),1);
out_frame=zeros(sum(dists(i,j)~=inf && dists(i,j)>(avg+ConfidenceInterval*stdev)),1);
cnt0=1;
% out_frame=[];
% out_nucleus=[];
for i=1:size(dists,1)
    for j=1:size(dists,2)
        density(i)=sum(dists(i,:)~=inf & dists(i,:)~=0 & dists(i,:)>(avg+ConfidenceInterval*stdev))/sum(dists(i,:)~=inf & dists(i,:)~=0);
        if (dists(i,j)~=inf && dists(i,j)>(avg+ConfidenceInterval*stdev))
            outliers(i,j)=dists(i,j);
            out_frame(cnt0)=i;
            out_nucleus(cnt0)=j;
            cnt0=cnt0+1;
        end
    end
end

%% Decide on the midline
out_posnx=zeros(length(out_frame));
out_posny=zeros(length(out_frame));
numnuc=length(dists)/2;
for i=1:length(out_frame)
    out_posnx(i)=Ellipses{out_frame(i),1}(mod(out_nucleus(i),numnuc),1);
    out_posny(i)=Ellipses{out_frame(i),1}(mod(out_nucleus(i),numnuc),2);
end

meanys=zeros(NumberOfFrames);
for i=cf:NumberOfFrames
    yposns=out_posny(out_frame==i);
    meanys(i)=mean(yposns);
    stdev2=std(yposns);
    %    yposns=yposns(yposns~=0);
    if sum(yposns==0)~=0
        display(i)
    end
    yposns=yposns(abs(yposns-meanys(i))<2*stdev2);
    if length(yposns)>6
        [ids,centroids]=kmeans(yposns',2);
        if sum(ids==1)>0.3*length(yposns)
            meanys(i)=mean(centroids);
        else
            meanys(i)=mean(yposns);
        end
    else
        meanys(i)=mean(yposns);
    end
end

% Smoothen midline by time-averaging
meanys(meanys==0)=nan;
typosns=zeros(NumberOfFrames,1);
tnum=3; % Number of frames by which to average
for i=cf:NumberOfFrames
    den=0;
    if i<cf+tnum
        for j=cf:i+tnum
            typosns(i)=typosns(i)+meanys(j);
            den=den+1;
        end
        typosns(i)=typosns(i)/den;
    else if i>NumberOfFrames-tnum
            for j=i-tnum:NumberOfFrames
                
                typosns(i)=typosns(i)+meanys(j);
                den=den+1;
            end
            typosns(i)=typosns(i)/den;
        else
            for j=i-tnum:i+tnum
                typosns(i)=typosns(i)+meanys(j);
                den=den+1;
            end
            typosns(i)=typosns(i)/den;
        end
    end
end

%% Create a movie
% Start a movie to watch where the outliers occur
typosns(isnan(typosns))=0;
mov=VideoWriter('Midline');
mov.FrameRate=2;
open(mov);
count=1;
%random_color=RGB::randn();
for i=cf:NumberOfFrames
    I=imread(strcat(PreProcPath,filesep,prefix,filesep,prefix,'-His_',num2str(i),'.tif'));
    for j=1:sum(out_frame==i)
        circles=[out_posnx(count) out_posny(count) 15];
        %         if mod(count,2)==1
        %             random_color=RGB::random();
        %         end
        I=insertShape(I, 'circle', circles);%, 'Color', RGB::random());
        if ~isnan(typosns(i))
            ln=[0 typosns(i) LinesPerFrame typosns(i)];
            I=insertShape(I, 'Line', ln, 'Color', 'r');
        end
        if (exist('lrep','var'))
            for ii=1:length(lrep)
                if ~isnan(lrep(ii))
                rep=[0 lrep(ii) LinesPerFrame lrep(ii)];
                I=insertShape(I, 'Line', rep, 'Color', 65536/2*[0 ii 0]);
                end
            end
        end
        count=count+1;
    end
    imshow(I)
    truesize([512,512]);
    % title('Movie of outliers');
    
    writeVideo(mov,getframe(gca));
end
close(mov)