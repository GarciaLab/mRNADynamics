function GenerateFigures(varargin)

%Setting the input variables
nc=varargin{1};
data_set=varargin{2};
data=LoadMS2Sets(data_set);
try APBinLowerrange=varargin{3};
end
try APBinUpperrange=varargin{4};
catch
    disp('Since no lower limit was provided, I''ll assume the range starts at 0 if that''s alright Sir')
end
%Parameters:
GeneLength5=5.296;      %Distance from the first MS2 site to the end of the
                        %TUB3'UTR in kb for the 5' constrcut.  
                        
%Get the FitResultsMultiple data        
[SourcePath,FISHPath,DropboxFolder,MS2CodePath,PreProcPath]=...
    DetermineLocalFolders;

%This is all to read the DataStatus file and pull the right folders. Errors
%here are likely due to discrepancies in file names versus the excel sheet
[xldata, xltext] = xlsread('DataStatus.xlsx',data_set);
dimensionsofxltext=size(xltext);
Prefixes=[];
for i=1:dimensionsofxltext(2)-1
    prefixi=xltext{6,i+1};
    prefixi=cellstr(prefixi);
    Prefixes=[Prefixes,prefixi];
end
i=1;
checkdataset=Prefixes{i}(22:21+length(data_set));
samedataset=isequal(data_set,checkdataset);
while samedataset==0
    i=i+1;
    checkdataset=Prefixes{i}(22:21+length(data_set));
    samedataset=isequal(data_set,checkdataset);
end
folderlocation=Prefixes{i}(11:end-1);
load([DropboxFolder,filesep,folderlocation,'\FitResultsMultiple.mat']);
 
%% Average fit results
%We'll assume that mRNA's are terminated in a non-random way with a certain
%time delay. The idea is to be able to calculate what part of the average
%signal change is due to transcription initiation by removing the part
%corresponding to termination.

close all
% 
% Constructs arrays of integrated mRNA production and associated error,
% respectively, using formula (pol loading rate) * (transcriptional time
% window = toff - ton).


%IGNORE THIS FOR NOW
% NumberOfDatasets=1+numberNCs(2);
% Alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
% letter=Alphabet(NumberOfDatasets);
% range=strcat('B6:',letter,'6');
% filelocation=xlsread([DropboxFolder,filesep,'DataStatus.xlsx'],data_set);

%This creates all the data needed for histograms of TimeStart and TimeEnd
numberNCs=size(data(1).MeanFits);
AllTimeStarts=[];
AllTimeEnds=[];
AllDurations=[];
lowerlimitexist=exist('APBinLowerrange','var'); %Checking existance of the inputs
upperlimitexist=exist('APBinUpperrange','var');

%This runs the code for two bounds provided
if upperlimitexist==1
    for dsn=1:length(FitResultsMultiple); %dsn is short for data set number
        for particlenumber=1:length(FitResultsMultiple{1,dsn})
            for cycle=1:nc-11
                if FitResultsMultiple{1,dsn}(particlenumber,cycle).Approved==1
                    TimeStart=FitResultsMultiple{1,dsn}(particlenumber,cycle).TimeStart;
                    TimeEnd=FitResultsMultiple{1,dsn}(particlenumber,cycle).TimeEnd;
                    APBin=(FitResultsMultiple{1,dsn}(particlenumber,cycle).APBin)/40;
                    if (APBin>varargin{3}) && (APBin>varargin{3}) && (APBin<varargin{4})
                        if isempty(AllTimeStarts)==1
                            AllTimeStarts=TimeStart;
                            AllTimeEnds=TimeEnd;
                            AllDurations=TimeEnd-TimeStart;
                        else
                            AllTimeStarts=cat(2,AllTimeStarts,TimeStart);
                            AllTimeEnds=cat(2,AllTimeEnds,TimeEnd);
                            AllDurations=cat(2,AllDurations,TimeEnd-TimeStart);
                        end
                    end
                end
            end
        end
    end
%This runs the code if only one bound is provided, and assumes lower bound
%is 0
elseif lowerlimitexist==1;
    for dsn=1:length(FitResultsMultiple); %dsn is short for data set number
        for particlenumber=1:length(FitResultsMultiple{1,dsn})
            for cycle=1:nc-11
                if FitResultsMultiple{1,dsn}(particlenumber,cycle).Approved==1
                    TimeStart=FitResultsMultiple{1,dsn}(particlenumber,cycle).TimeStart;
                    TimeEnd=FitResultsMultiple{1,dsn}(particlenumber,cycle).TimeEnd;
                    APBin=(FitResultsMultiple{1,dsn}(particlenumber,cycle).APBin)/40;
                    if (APBin~=0) && (APBin<varargin{3})
                        if isempty(AllTimeEnds)==1
                            AllTimeStarts=TimeStart;
                            AllTimeEnds=TimeEnd;
                            AllDurations=TimeEnd-TimeStart;
                        else
                            AllTimeStarts=cat(2,AllTimeStarts,TimeStart);
                            AllTimeEnds=cat(2,AllTimeEnds,TimeEnd);
                            AllDurations=cat(2,AllDurations,TimeEnd-TimeStart);
                        end
                    end
                end
            end
        end
    end
end

%This only plots the Histograms if a bound or bounds were provided
if lowerlimitexist==1;
    figure('units','normalized','position', [.1 .25 .35 .6]); % create new figure
    title(num2str(nc))
    subplot(3,1,1); % first subplot
    histogram(AllTimeStarts,8);
    title(['TimeStart (nc',num2str(nc),')']);
    xlabel('Time [min]','FontSize',12);
    ylabel('Number of Particles','FontSize',12);
    subplot(3,1,2);
    histogram(AllTimeEnds,8);
    title(['TimeEnd (nc',num2str(nc),')']);
    xlabel('Time [min]','FontSize',12);
    ylabel('Number of Particles','FontSize',12);
    subplot(3,1,3);
    histogram(AllDurations,8);
    title(['Duration (nc',num2str(nc),')']);
    xlabel('Time [min]','FontSize',12);
    ylabel('Number of Particles','FontSize',12);
end 

% This plots the TimeStart and TimeEnd data as a funtion of AP position
APvsTimeStart=[];
APvsTimeEnd=[];
for apbin=1:41;
    if data(1).MeanFits(apbin,cycle).TimeEnd~=0
        APvsTimeEnd=cat(2,APvsTimeEnd,data(1).MeanFits(apbin,cycle).TimeEnd);
        APvsTimeStart=cat(2,APvsTimeStart,data(1).MeanFits(apbin,cycle).TimeStart);
    else
        data(dsn).MeanFits(apbin,cycle).Approved=0;
        APvsTimeEnd=cat(2,APvsTimeEnd,NaN);
        APvsTimeStart=cat(2,APvsTimeStart,NaN);
    end
end
figure('units','normalized','position', [.3 .15 .3 .5]); % create new figure
title(num2str(nc))
subplot(2,1,1); % first subplot
x=linspace(0,1,41);
plot(x,APvsTimeStart,'.-','Color','r');
title(['TimeStart (nc',num2str(nc),')']);
xlim([0,1]);
ylabel('Time [min]','FontSize',12);
xlabel('AP Position','FontSize',12);
subplot(2,1,2); % second subplot
plot(x,APvsTimeEnd,'.-','Color','r');
xlim([0,1]);
title(['TimeEnd (nc',num2str(nc),')']);
ylabel('Time [min]','FontSize',12);
xlabel('AP Position','FontSize',12);
            
mRNA_AP=nan(length(data(1).APbinID),numberNCs(2));
SD_mRNA_AP=nan(length(data(1).APbinID),numberNCs(2));
for j=1:length(data(1).APbinID)
    for k=1:numberNCs(2)
        if (~isempty(data(1).MeanFits(j,k).RateFit))&(data(1).MeanFits(j,k).Approved==1)
            syms r
            rate = data(1).MeanFits(j,k).RateFit;
            t_off = data(1).MeanFits(j,k).TimeEnd;
            t_on = data(1).MeanFits(j,k).TimeStart;
            [mRNA_AP(j,k), SD_mRNA_AP(j,k)] = PropError(r*(t_off-t_on),[r], [rate], [data(1).MeanFits(j,k).SDRateFit]);
        end
    end
end

%Constructs arrays of integrated mRNA production and associated error,
%respectively, using formula found in SI of (Bothma, 2014).
% 
% for i=1:length(data)
%     mRNA_AP_J=nan(length(data),length(data(i).APbinID),3);
%     SD_mRNA_AP_J=nan(length(data),length(data(i).APbinID),3);
%     for j=1:length(data(i).APbinID)
%         for k=1:2
%             if (~isempty(data(i).MeanFits(j,k).RateFit))&(data(i).MeanFits(j,k).Approved==1)
%                 rate = data(i).MeanFits(j,k).RateFit;
%                 delay = GeneLength5/rate;
%                 t_off = data(i).MeanFits(j,k).TimeEnd;
%                 t_on = data(i).MeanFits(j,k).TimeStart;
%                 syms r toff ton;
%                 [ mRNA_AP_J(i,j,k), SD_mRNA_AP_J(i,j,k)] = PropError( (2 / (GeneLength5/r)) * (.5*( (toff - GeneLength5/r) + (toff - ton + GeneLength5/r) ) * GeneLength5),[toff, ton, r], [t_off, t_on, rate], [data(i).MeanFits(j,k).SDTimeEnd,data(i).MeanFits(j,k).SDTimeStart, data(i).MeanFits(j,k).SDRateFit]);
%             end
%         end
%     end
% end
 
%Generate the weighted averages for the rate
% for i=1:length(data)
% mRNA_weight=nan(length(data(1).APbinID),3);
% SD_mRNA_weight=nan(length(data(1).APbinID),3);
% SE_mRNA_weight=nan(length(data(1).APbinID),3);
% for j=1:length(data(1).APbinID)
%     for k=1:3
%         NanFilter=~(isnan(mRNA_AP(:,j,k))|isnan(SD_mRNA_AP(:,j,k)));
%         
%         if sum(NanFilter)>=MinEmbryos
%        
%             mRNA_weight(j,k)=mean(mRNA_AP(NanFilter,j,k));
%             SD_mRNA_weight(j,k)=std(mRNA_AP(NanFilter,j,k));
%             SE_mRNA_weight(j,k)=std(mRNA_AP(NanFilter,j,k))/sqrt(sum(NanFilter));
%         end
%     end
% end 
 
   
 
%Overlay single rates and mean rate
figure('units','normalized','position', [.5 .25 .4 .6])
clf
PlotHandle=[];
hold on
%Make the zeros in the data set NaN so they do not plot
dimmRNA_AP=size(mRNA_AP);
for j=1:dimmRNA_AP(1)
        for k=1:dimmRNA_AP(2)
            if mRNA_AP(j,k)==0
                mRNA_AP(j,k)=NaN;
            end
        end
end
% Plot with error bars
PlotHandle(end+1)=errorbar(data(1).APbinID,mRNA_AP(:,nc-11),SD_mRNA_AP(:,nc-11),'.-','Color','r');
% PlotHandle(end+1)=scatter(data(1).APbinID,mRNA_AP(:,nc-11))

% PlotHandle(end+1)=plot(data(1).APbinID,mRNA_AP(:,nc-11),'Color','r');
% PlotHandle(end+1)=errorbar(data(1).APbinID,mRNAWeight(:,nc-11),SDWeight(:,nc-11),'.-r');
xRange=linspace(0,1);
% plot(xRange,ones(size(xRange))*(RateNoBcdWeight(nc-11)+SERateNoBcdWeight(nc-11)),'--r')
% plot(xRange,ones(size(xRange))*(RateNoBcdWeight(nc-11)-SERateNoBcdWeight(nc-11)),'--r')
hold off
box on
xlim([0,1])
ylim([0,1.25*max(mRNA_AP(:, nc-11))])
xlabel('AP position (x/L)','FontSize',12)
ylabel('Accumulated mRNA (fu)','FontSize',12)
%StandardFigure(PlotHandle,gca)
title(['mRNA Produced in nc',num2str(nc)])
disp('There you are, Sir. I hope these results are to your liking')
 