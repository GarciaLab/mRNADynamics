%% This script requires compiled particles and the lineage
% This script will look at the number of particles expressed per nucleus
% before and after mitosis. There will also be plots looking at the
% intensity of expression before and after mitosis.
function combos=twospot(Prefix,Ellipses,schnitzcells,CompiledParticles,PreProcPath, nc14,nc13,delay)
close all
% Extract Data
% load(['/home/avaneesh/Dropbox/LivemRNAEmilia',filesep,Prefix,filesep,Prefix,'_lin.mat']);
% load(['/home/avaneesh/Dropbox/LivemRNAEmilia',filesep,Prefix,filesep,'Ellipses.mat']);
% load(['/home/avaneesh/Dropbox/LivemRNAEmilia',filesep,Prefix,filesep,'CompiledParticles.mat']);
% [SourcePath,FISHPath,DefaultDropboxFolder,MS2CodePath,PreProcPath]=...
%     DetermineLocalFolders;
% [Num,Txt,XLSRaw]=xlsread([DefaultDropboxFolder,filesep,'MovieDatabase.xlsx']);
% nc14Column=find(strcmp(XLSRaw(1,:),'nc14'));
% DataFolderColumn=find(strcmp(XLSRaw(1,:),'DataFolder'));
% Dashes=findstr(Prefix,'-');
% XLSEntry=find(strcmp(XLSRaw(:,DataFolderColumn),[Prefix(1:Dashes(3)-1),'\',Prefix(Dashes(3)+1:end)]));
% nc14=XLSRaw{XLSEntry,nc14Column};

%% Initialise Variables
firstnc14=size(schnitzcells,2); % First schnitz to appear in nc14
first14=size(Ellipses,1);
firstnc13=size(schnitzcells,2); % First schnitz to appear in nc13
first13=size(Ellipses,1);
for i=1:size(schnitzcells,2)
    if isempty(schnitzcells(1,i).frames)
        break
    end
    appears=schnitzcells(1,i).frames(1);
    if (appears<first14 & appears>nc14-2)
        firstnc14=i;
        first14=appears;
    end
    if (appears<first13 & appears>nc13-2)
        firstnc13=i;
        first13=appears;
    end
end

%% How good is the lineage?
len=0;
for i=firstnc13:firstnc14
    if schnitzcells(1,i).D>0  | schnitzcells(1,i).E>0
        len=len+1;
    end
end
display('The proportion of detected families is')
display((len/(firstnc14-firstnc13)))

%% Acquire Families
families=zeros(len,3);
count=0;

h=waitbar(0,'Acquiring Families');
for i=firstnc13:firstnc14
    waitbar(i/(-firstnc13+firstnc14));
    if (schnitzcells(1,i).D>0)
        count=count+1;
        if ~isempty(schnitzcells(1,i).E) && (schnitzcells(1,i).E~=0)
            families(count,:)=[i schnitzcells(1,i).E schnitzcells(1,i).D]';
        else
            families(count,:)=[i 0 schnitzcells(1,i).D]';
        end
    elseif (schnitzcells(1,i).E>0)
        count=count+1;
        families(count,:)=[i schnitzcells(1,i).E 0]';
    end
end
close(h);

%% Tag each family according to the combination
combs=zeros(len,1);
combos=zeros(len,3);
intens=zeros(len,3);
nucs=zeros(len,3);
h=waitbar(0,'Acquiring info about families');
for i=1:len
    waitbar(i/len);
    num1=0;
    num2=0;
    num3=0;
    for ii=1:size(CompiledParticles,2)
        if(CompiledParticles(ii).Nucleus==families(i,1))
            num1=num1+1;
            intens(i,1)=intens(i,1)+sum(CompiledParticles(ii).Fluo)/(nc14-nc13);
            nucs(i,1)=ii;
        elseif(CompiledParticles(ii).Nucleus==families(i,2))
            if any(CompiledParticles(ii).Frame==nc14+delay)
                num2=num2+1;
                intens(i,2)=intens(i,2)+sum(CompiledParticles(ii).Fluo)/(size(Ellipses,1)-nc14);
            end
            nucs(i,2)=ii;
        elseif(CompiledParticles(ii).Nucleus==families(i,3))
            if any(CompiledParticles(ii).Frame==nc14+delay)
                num3=num3+1;
                intens(i,3)=intens(i,3)+sum(CompiledParticles(ii).Fluo)/(size(Ellipses,1)-nc14);
            end
            nucs(i,3)=ii;
        end
    end
    combos(i,:)=[num1 num2 num3];
    combs(i)=100*num1+10*num2+num3;
end
spots=max(max(combos));
close(h);

%% Display Data
%try
    i=nc14+delay;
    figure(delay);
    I=imread(strcat(PreProcPath,filesep,Prefix,filesep,Prefix,'-His_',num2str(iIndex(i,3)),'.tif'));
    failed=0;
    for ii=1:len
        if(~families(ii,2)==0)
            if isempty(schnitzcells(1,families(ii,2)).frames)
                break
            end
            frame=find(schnitzcells(1,families(ii,2)).frames==i);
            if(isempty(frame))
                frame=1;
                %display('Couldnt find frame');
                failed=failed+1;
                
            end
            if ~isempty(schnitzcells(1,families(ii,2)).cenx)
                circles=[schnitzcells(1,families(ii,2)).cenx(frame) schnitzcells(1,families(ii,2)).ceny(frame) 15];
                I=insertShape(I, 'circle', circles,'Color',65535*combos(ii,:)/spots);
            end
            % display(combos(ii,:));
            % display('Inserted!')
            
        end
        if(~families(ii,3)==0)
            frame=find(schnitzcells(1,families(ii,3)).frames==i);
            if(isempty(frame))
                frame=1;
                %display('Couldnt find frame');
                failed=failed+1;
            end
            circles=[schnitzcells(1,families(ii,3)).cenx(frame) schnitzcells(1,families(ii,3)).ceny(frame) 15];
            I=insertShape(I, 'circle', circles,'Color',65535*combos(ii,:)/spots);
            % display(combos(ii,:));
            % display('Inserted!')
        end
    end
    imshow(I)
    
    filename=strcat(Prefix,delay,'_2spot.png');
    saveas(gca,filename);
% catch
%     display(['Histone pictures not present. The distribution picture cannot be shown.']);
% end

%% Intensity Plots
figure;
intens2(:,1)=intens(:,1);
intens2(:,2)=intens(:,2)+intens(:,3);
counter=0;
for i=1:length(intens)
    counter=counter+1;
    plot(1:2,intens2(i,:),'LineWidth',2);
    hold on;
end
filename=strcat(Prefix,'_2spot_intens.png');
saveas(gca,filename);

%% All Plots
count1=0;
count2=0;
count3=0;
count4=0;
for ii=1:len
    if(combs(ii)==101)
        %if count3==15
        figure(101)
        hold on
        try
            plot(CompiledParticles(nucs(ii,1)).Frame-nc14,CompiledParticles(nucs(ii,1)).Fluo);
        catch
            display(['Problem with family ',num2str(ii),' with combination ',num2str(combs(ii))]);
        end
        hold on
        try
            plot(CompiledParticles(nucs(ii,3)).Frame-nc14,CompiledParticles(nucs(ii,2)).Fluo);
        catch
            display(['Problem with family ',num2str(ii),' with combination ',num2str(combs(ii))]);
        end
        hold on
        try
            plot(CompiledParticles(nucs(ii,3)).Frame-nc14,CompiledParticles(nucs(ii,3)).Fluo);
        catch
            display(['Problem with family ',num2str(ii),' with combination ',num2str(combs(ii))]);
        end
        %end
        count3=count3+1;
    end
    if(combs(ii)==110)
        %if count4==15
        figure(110)
        hold on
        try
            plot(CompiledParticles(nucs(ii,1)).Frame-nc14,CompiledParticles(nucs(ii,1)).Fluo);
        catch
            display(['Problem with family ',num2str(ii),' with combination ',num2str(combs(ii))]);
        end
        hold on
        try
            plot(CompiledParticles(nucs(ii,3)).Frame-nc14,CompiledParticles(nucs(ii,2)).Fluo);
        catch
            display(['Problem with family ',num2str(ii),' with combination ',num2str(combs(ii))]);
        end
        hold on
        try
            plot(CompiledParticles(nucs(ii,3)).Frame-nc14,CompiledParticles(nucs(ii,3)).Fluo);
        catch
            display(['Problem with family ',num2str(ii),' with combination ',num2str(combs(ii))]);
        end
        hold on
        %end
        count4=count4+1;
    end
    
    if(combs(ii)==111)
        %if count1==30
        figure(111)
        hold on
        try
            plot(CompiledParticles(nucs(ii,1)).Frame-nc14,CompiledParticles(nucs(ii,1)).Fluo);
        catch
            display(['Problem with family ',num2str(ii),' with combination ',num2str(combs(ii))]);
        end
        hold on
        try
            plot(CompiledParticles(nucs(ii,3)).Frame-nc14,CompiledParticles(nucs(ii,2)).Fluo);
        catch
            display(['Problem with family ',num2str(ii),' with combination ',num2str(combs(ii))]);
        end
        hold on
        try
            plot(CompiledParticles(nucs(ii,3)).Frame-nc14,CompiledParticles(nucs(ii,3)).Fluo);
        catch
            display(['Problem with family ',num2str(ii),' with combination ',num2str(combs(ii))]);
        end
        
        %end
        count1=count1+1;
    end
    if(combs(ii)==100)
        %if count2==35
        figure(100)
        hold on
        try
            plot(CompiledParticles(nucs(ii,1)).Frame-nc14,CompiledParticles(nucs(ii,1)).Fluo);
        catch
            display(['Problem with family ',num2str(ii),' with combination ',num2str(combs(ii))]);
        end
        hold on
        try
            plot(CompiledParticles(nucs(ii,3)).Frame-nc14,CompiledParticles(nucs(ii,3)).Fluo);
        catch
            display(['Problem with family ',num2str(ii),' with combination ',num2str(combs(ii))]);
        end
        
        hold on
        try
            plot(CompiledParticles(nucs(ii,3)).Frame-nc14,CompiledParticles(nucs(ii,2)).Fluo);
        catch
            display(['Problem with family ',num2str(ii),' with combination ',num2str(combs(ii))]);
        end
        %end
        count2=count2+1;
    end
end
display(count1/len);
display(count2/len);
display((count3+count4)/len);
figure(314);
pie([count1,count2,count3+count4]);
legend('111','100','110');
saveas(314,'PieOfLineageBreakdown.png');
saveas(100,'DaughterCombo100.png');
saveas(111,'DaughterCombo111.png');
saveas(101,'DaughterCombo101.png');
saveas(110,'DaughterCombo110.png');
%display(failed/len);