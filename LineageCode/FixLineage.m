function  schnitzcells=FixLineage(Prefix, FrameInfo)
% Finds the times of cell division (by using nuclei radius,
% and orphan the nuclei if they survive mitosis)

% Note: This function has many nested functions, which are a bit
% complicated in Matlab.

% TrackNucleiDV(Prefix);

%% Extract all data
%Get the folders, including the default Dropbox one
[SourcePath,FISHPath,DefaultDropboxFolder,MS2CodePath,PreProcPath]=...
    DetermineLocalFolders;
%Now get the actual DropboxFolder
[SourcePath,FISHPath,DropboxFolder,MS2CodePath,PreProcPath]=...
    DetermineLocalFolders(Prefix);



%Determine division times
%Load the information about the nc from the XLS file
[Num,Txt,XLSRaw]=xlsread([DefaultDropboxFolder,filesep,'MovieDatabase.xlsx']);
XLSHeaders=Txt(1,:);
Txt=Txt(2:end,:);

ExperimentTypeColumn= strcmp(XLSRaw(1,:),'ExperimentType');
ExperimentAxisColumn=find(strcmp(XLSRaw(1,:),'ExperimentAxis'));

DataFolderColumn=find(strcmp(XLSRaw(1,:),'DataFolder'));
Dashes=findstr(Prefix,'-');

PrefixRow=find(strcmp(XLSRaw(:,DataFolderColumn),[Prefix(1:Dashes(3)-1),'\',Prefix(Dashes(3)+1:end)]));
if isempty(PrefixRow)
    PrefixRow=find(strcmp(XLSRaw(:,DataFolderColumn),[Prefix(1:Dashes(3)-1),'/',Prefix(Dashes(3)+1:end)]));
    if isempty(PrefixRow)
        error('Could not find data set in MovieDatabase.XLSX. Check if it is defined there.')
    end
end

ExperimentType=XLSRaw{PrefixRow,ExperimentTypeColumn};
ExperimentAxis=XLSRaw{PrefixRow,ExperimentAxisColumn};

%Find the different columns.
DataFolderColumn=find(strcmp(XLSRaw(1,:),'DataFolder'));
nc9Column= strcmp(XLSRaw(1,:),'nc9');
nc10Column= strcmp(XLSRaw(1,:),'nc10');
nc11Column= strcmp(XLSRaw(1,:),'nc11');
nc12Column= strcmp(XLSRaw(1,:),'nc12');
nc13Column= strcmp(XLSRaw(1,:),'nc13');
nc14Column= strcmp(XLSRaw(1,:),'nc14');
CFColumn= strcmp(XLSRaw(1,:),'CF');
Channel1Column=find(strcmp(XLSRaw(1,:),'Channel1'));
Channel2Column=find(strcmp(XLSRaw(1,:),'Channel2'));


%Find the corresponding entry in the XLS file
if (~isempty(findstr(Prefix,'Bcd')))&&(isempty(findstr(Prefix,'BcdE1')))&&...
        (isempty(findstr(Prefix,'NoBcd')))&&(isempty(findstr(Prefix,'Bcd1x')))&&(isempty(findstr(Prefix,'Bcd4x')))
    warning('This step in CheckParticleTracking will most likely have to be modified to work')
    XLSEntry=find(strcmp(XLSRaw(:,DataFolderColumn),...
        [Date,'\BcdGFP-HisRFP']));
else
    XLSEntry=find(strcmp(XLSRaw(:,DataFolderColumn),...
        [Prefix(1:Dashes(3)-1),'\',Prefix(Dashes(3)+1:end)]));
    
    if isempty(XLSEntry)
        XLSEntry=find(strcmp(XLSRaw(:,DataFolderColumn),...
            [Prefix(1:Dashes(3)-1),'/',Prefix(Dashes(3)+1:end)]));
        if isempty(XLSEntry)
            disp('%%%%%%%%%%%%%%%%%%%%%')
            error('Dateset could not be found. Check MovieDatabase.xlsx')
            disp('%%%%%%%%%%%%%%%%%%%%%')
        end
    end
end



nc9 =XLSRaw{XLSEntry,nc9Column};
nc10=XLSRaw{XLSEntry,nc10Column};
nc11=XLSRaw{XLSEntry,nc11Column};
nc12=XLSRaw{XLSEntry,nc12Column};
nc13=XLSRaw{XLSEntry,nc13Column};
nc14=XLSRaw{XLSEntry,nc14Column};
CF  =XLSRaw{XLSEntry,CFColumn};

%This checks whether all ncs have been defined
ncCheck=[nc9,nc10,nc11,nc12,nc13,nc14];
if length(ncCheck)~=6
    error('Check the nc frames in the MovieDatabase entry. Some might be missing')
end

%Do we need to convert any NaN chars into doubles?
if strcmpi(nc14,'nan')
    nc14=nan;
end
if strcmpi(nc13,'nan')
    nc13=nan;
end
if strcmpi(nc12,'nan')
    nc12=nan;
end
if strcmpi(nc11,'nan')
    nc11=nan;
end
if strcmpi(nc10,'nan')
    nc10=nan;
end
if strcmpi(nc9,'nan')
    nc9=nan;
end


%Convert the prefix into the string used in the XLS file
Dashes=findstr(Prefix,'-');

%Find the corresponding entry in the XLS file
if (~isempty(findstr(Prefix,'Bcd')))&(isempty(findstr(Prefix,'BcdE1')))&...
        (isempty(findstr(Prefix,'NoBcd')))&(isempty(findstr(Prefix,'Bcd1')))&(isempty(findstr(Prefix,'Bcd4x')))
    XLSEntry=find(strcmp(Txt(:,DataFolderColumn),...
        [Date,'\BcdGFP-HisRFP']));
else
    XLSEntry=find(strcmp(Txt(:,DataFolderColumn),...
        [Prefix(1:Dashes(3)-1),filesep,Prefix(Dashes(3)+1:end)]));
end
ncs=[nc9,nc10,nc11,nc12,nc13,nc14];

if (length(find(isnan(ncs)))==length(ncs))||(length(ncs)<6)
    error('Have the ncs been defined in MovieDatabase.XLSX?')
end

%Now do the nuclear segmentation and lineage tracking. This should be put
%into an independent function.




%% Load the lineage and save a backup
% (it is recommende to use Dropbox for version control)

load([DropboxFolder,filesep,Prefix,filesep,Prefix,'_lin.mat']);

cpexist=true; % Does Compiled Particles Exist?
try
    load([DropboxFolder,filesep,Prefix,filesep,'Particles.mat']);
    save([DropboxFolder,filesep,Prefix,filesep,'Particles_backup.mat'],Particles);
    cpexist=true;
catch
    cpexist=false;
end
load([DropboxFolder,filesep,Prefix,filesep,'FrameInfo.mat']);

% The following are counter variables that are used to track how many times
% parts of the code have been run. So basically, they serve no utility

DaughtersOlderThanParents=0;  % Number of daughter cells with associated parents that appear in a frame later than themselves. Obvious error, should be 0
NumberOfStitchesPerformed=0; % Number of times pairs of traces of nuclei and associated particles have been stitched
SortErrors=0; % The number of times the schnitzcells after sort are not the same size as the schnitzcells before
NumberOfStitchMismatches=0; % The number of times the traces being stitched have different parent cells. Some of the parent information will be lost

% schnitz_stitch is a nested function that stitches traces that are close
% by and such that the beginning of the second comes right after the end of
% the first
schnitzcells=schnitz_stitch(schnitzcells,ncs);

% schnitz_sort is a nested function that sorts the lineage data according
% to the first frame in which they appear
schnitzcells=schnitz_sort(schnitzcells);
AdditionalSchnitzes=0; % Number of new schnitzcell elements introduced
NumberOfSchnitzes=size(schnitzcells,2); % Original Number of schnitzcell elements


%% Remove nuclei that appear and disappear
% This section removes the nuclei that appear only for a short while at the
% edges, before disappearing beyond the frame of view

margins=5; % Choose (in number of pixels) the space to be ignored when fixing lineage

% Instead of a for loop, I use a while loop iterating over the index
% variable i
i=1;
while(i<NumberOfSchnitzes)
    
    % We can choose how short a trace needs to be before it is suspect
    if length(schnitzcells(i).frames)<5
        
        % We only look at the schnitzes that disappear at the margins
        if (any(schnitzcells(i).cenx<margins) || ...
                any(schnitzcells(i).cenx>(FrameInfo(1).LinesPerFrame-margins)) ||...
                any(schnitzcells(i).ceny<margins) ||...
                any(schnitzcells(i).ceny>(FrameInfo(1).LinesPerFrame-margins)))
            
            % Instead of deleting, we make the frames an Inf value.
            % We later remove all such schnitzes. This is done so as not
            % to mess up the indices
            schnitzcells(i).frames=Inf;
            
            % If the schnitz has associated parents or daughters, an empty
            % matrix is put in its place
            if ~isempty(schnitzcells(i).P)
                parent_del=schnitzcells(i).P;
                if schnitzcells(parent_del).D==i
                    schnitzcells(parent_del).D=[];
                elseif schnitzcells(parent_del).E==i
                    schnitzcells(parent_del).E=[];
                end
            end
            if ~isempty(schnitzcells(i).D)
                child_del=schnitzcells(i).D;
                schnitzcells(child_del).P=[];
            end
            if ~isempty(schnitzcells(i).E)
                child_del=schnitzcells(i).E;
                schnitzcells(child_del).P=[];
            end
        end
    end
    i=i+1;
end

%% Schnitz Splits
% The schnitz traces are now split if they go across a nuclear cycle
% There are a few parameters to tweak here that are dependent on
% TimeResolution

TimeResolution=20;
i=1;
while(i<NumberOfSchnitzes)
    % Check ncs
    t=i;
    for nc=1:length(ncs)
        if isempty(schnitzcells(i).frames)
            schnitzcells(i).frames=Inf;
        end
        if (min(schnitzcells(i).frames) <ncs(nc)-(2*20/TimeResolution)) && ...
                (max(schnitzcells(i).frames) >ncs(nc)+(8*20/TimeResolution))
            
            AdditionalSchnitzes=AdditionalSchnitzes+1; % This is only for counting
            schnitzcells(NumberOfSchnitzes+1)=schnitzcells(i);
            schnitzcells(i).E=NumberOfSchnitzes+1;
            schnitzcells(i).D=[];
            
            schnitzcells(NumberOfSchnitzes+1).P=i;
            % The schnitz trace before the split
            FramesBeforeSplit=schnitzcells(i).frames(schnitzcells(i).frames<ncs(nc));
            [~,SplitFrame]=max(FramesBeforeSplit); % SplitFrame is the exact frame where the split takes place
            ParticleIndices=[]; % Particles associated with the current nucleus
            if cpexist
                for ii=1:size(Particles,2)
                    if  Particles(ii).Nucleus==i
                        ParticleIndices=[ParticleIndices ii];
                    end
                end
            end
            
            temp=schnitzcells(i);
            
            
            % The splitting, with corner cases for splits at beginning or
            % end
            if ~isempty(temp.frames(1:SplitFrame))
                schnitzcells(i).frames=temp.frames(1:SplitFrame);
                schnitzcells(i).cenx=temp.cenx(1:SplitFrame);
                schnitzcells(i).ceny=temp.ceny(1:SplitFrame);
                schnitzcells(i).len=temp.len(1:SplitFrame);
                schnitzcells(i).cellno=temp.cellno(1:SplitFrame);
                
                
                
                if ~isempty(temp.frames(SplitFrame+1:end))
                    schnitzcells(NumberOfSchnitzes+1).frames=temp.frames(SplitFrame+1:end);
                    schnitzcells(NumberOfSchnitzes+1).cenx=temp.cenx(SplitFrame+1:end);
                    schnitzcells(NumberOfSchnitzes+1).ceny=temp.ceny(SplitFrame+1:end);
                    schnitzcells(NumberOfSchnitzes+1).len=temp.len(SplitFrame+1:end);
                    schnitzcells(NumberOfSchnitzes+1).cellno=temp.cellno(SplitFrame+1:end);
                    for jj=1:length(ParticleIndices)
                        if any(Particles(ParticleIndices(jj)).Frame>SplitFrame)
                            Particles(ParticleIndices(jj)).Nucleus=NumberOfSchnitzes+1;
                        end
                    end
                end
                
            else
                schnitzcells(i).frames=temp.frames(SplitFrame+1:end);
                schnitzcells(i).cenx=temp.cenx(SplitFrame+1:end);
                schnitzcells(i).ceny=temp.ceny(SplitFrame+1:end);
                schnitzcells(i).len=temp.len(SplitFrame+1:end);
                schnitzcells(i).cellno=temp.cellno(SplitFrame+1:end);
                for jj=1:length(ParticleIndices)
                    if any(Particles(ParticleIndices(jj)).Frame>SplitFrame)
                        Particles(ParticleIndices(jj)).Nucleus=NumberOfSchnitzes+1;
                    end
                end
            end
            % If split, move the index back one schnitz, in order to do
            % additional splits on some schnitzes if needed
            t=i-1;
        end
    end
    
    NumberOfSchnitzes=size(schnitzcells,2);
    i=max(t+1,0);
    
end

%% Sort
schnitzcells=schnitz_sort(schnitzcells);
i=1;
while i<size(schnitzcells,2)
    if isinf(schnitzcells(i).frames(1))
        schnitzcells(i)=[];
        i=i-1;
    end
    i=i+1;
end
display(NumberOfStitchMismatches);
display(DaughtersOlderThanParents);
display(NumberOfStitchesPerformed);
display(SortErrors);

    %% Nested function to stitch schnitzes that are too short
    function schnitz=schnitz_stitch(schnitz,ncs)
        NumberOfSchnitzes=size(schnitz,2);
        k=1;
        while k<NumberOfSchnitzes
            if ~isempty(schnitz(k).frames)
                nc=0; % Current nc
                for kk=1:length(ncs)
                    if schnitz(k).frames(1)>ncs(kk)
                        nc=kk;
                    end
                end
                Margins=30;
                try
                if (schnitz(k).cenx(1)>Margins && schnitz(k).cenx(1)<(FrameInfo(1).PixelsPerLine-Margins) ...
                    && schnitz(k).ceny(1)>Margins && schnitz(k).ceny(1)<(FrameInfo(1).LinesPerFrame-Margins))
                
                    %  Determine the minimum number of frames in which the schnitz should appear 
                    if nc<length(ncs)
                        mx_frame=ncs(nc+1)-2;
                    else
                        mx_frame=size(FrameInfo,2);
                    end
                    mn_frame=ncs(nc)+8;
                    if length(schnitz(k).frames)<(mx_frame-mn_frame)
                        
                        % Find the closest matching nuclei in the immediate
                        % next frame
                        Closest=0;
                        min_distance=Inf;
                        kkk=1;
                        
                        % Go through all the schnitzes. This could be
                        % made faster if needed
                        while(kkk<NumberOfSchnitzes)
                            if (schnitz(kkk).frames(1)==schnitz(k).frames(end)+1) &&...
                                    (schnitz(kkk).frames(1)-mx_frame)>4
                                
                                dist=(schnitz(kkk).cenx(1)-schnitz(k).cenx(end))^2+...
                                    (schnitz(kkk).ceny(1)-schnitz(k).ceny(end))^2;
                                if dist<min_distance;
                                    min_distance=dist;
                                    Closest=kkk;
                                end
                                
                            end
                            kkk=kkk+1;
                            
                            % Do no consider proximity to yourself
                            if kkk==k
                                kkk=kkk+1;
                            end
                        end
                        
                        % If a closest nuclei has been found
                        if Closest~=0 && min_distance<2000
                            
                            schnitz(k).E=schnitz(Closest).E;
                            schnitz(Closest).E=[];
                            
                            schnitz(k).D=schnitz(Closest).D;
                            schnitz(Closest).D=[];
                            
                            schnitz(k).frames=[schnitz(k).frames; schnitz(Closest).frames];
                            schnitz(Closest).frames=Inf;
                            
                            schnitz(k).cenx=[schnitz(k).cenx schnitz(Closest).cenx];
                            schnitz(Closest).cenx=[];
                            
                            schnitz(k).ceny=[schnitz(k).ceny schnitz(Closest).ceny];
                            schnitz(Closest).ceny=[];
                            
                            schnitz(k).len=[schnitz(k).len schnitz(Closest).len];
                            schnitz(Closest).len=[];
                            
                            schnitz(k).cellno=[schnitz(k).cellno schnitz(Closest).cellno];
                            schnitz(Closest).cellno=[];
                            
                            if ~isempty(schnitz(Closest).P)
                                NumberOfStitchMismatches=NumberOfStitchMismatches+1;
                            end
                            if ~isempty(schnitz(Closest).D)
                                child=schnitz(Closest).D;
                                schnitz(child).P=Closest;
                            end
                            if ~isempty(schnitz(Closest).E)
                                child=schnitz(Closest).E;
                                schnitz(child).P=Closest;
                            end
                            if cpexist
                                for ii=1:length(Particles)
                                    if Particles(ii).Nucleus==Closest
                                        Particles(ii).Nucleus=k;
                                    end
                                end
                            end
                            NumberOfStitchesPerformed=NumberOfStitchesPerformed+1;
                        end
                    end
                    
                end
                catch
                end
            end
            k=k+1;
        end
    end

    %% Sorts the schnitzes using quicksort and changes the associated
    %  nuclei as necessary
    function schnitz=schnitz_sort(schnitz)
        lengt=size(schnitz,2);
        if lengt>2
            % Choose pivot
            pivot=ceil(lengt*rand);
            schnitz_swap(lengt,pivot);
            
            % 3-way partition
            ii = 1;
            k = 1;
            p = lengt;
            if isempty(schnitz(ii).frames)
                schnitz(ii).frames=Inf;
            end
            if isempty(schnitz(lengt).frames)
                schnitz(lengt).frames=Inf;
            end
            while ii < p
                if min(schnitz(ii).frames) < min(schnitz(lengt).frames)
                    schnitz_swap(ii,k);
                    ii=ii+1;
                    k=k+1;
                    
                elseif min(schnitz(ii).frames) == min(schnitz(lengt).frames)
                    p=p-1;
                    schnitz_swap(p,ii);
                    
                else ii=ii+1;
                end
            end
            % move pivots to center
            n=lengt;
            m = min(p-k,n-p+1);
            for ii=0:m-1
                schnitz_swap(k+ii,n-m+1+ii);
            end
            % recursive sorts
            schnitz_1=schnitz_sort(schnitz(1:k-1));
            schnitz_2=schnitz_sort(schnitz(n-p+k+1:n));
            schnitz(1:k-1)=schnitz_1;
            schnitz(k:n-p+k)=schnitz(k:n-p+k);
            schnitz(n-p+k+1:n)=schnitz_2;
        end
        
        if size(schnitz,2)~=lengt
            SortErrors=SortErrors+1;
            display('Loss in sort!')
        end
        
        function schnitz_swap(a,b)
            if a==b
                return;
            end
            swap=schnitz(b);
            schnitz(b)=schnitz(a);
            if ~isempty(schnitzcells(b).P)
                if(schnitzcells(b).P~=0)
                    parent=schnitzcells(b).P;
                    if schnitzcells(parent).E==b;
                        schnitzcells(parent).E=a;
                    elseif schnitzcells(parent).D==b;
                        schnitzcells(parent).D=a;
                    end
                end
            end
            if ~isempty(schnitzcells(b).E)
                if schnitzcells(b).E~=0
                    child=schnitzcells(b).E;
                    if schnitzcells(child).P==b;
                        if child>a
                            % display(['Downs everywhere!']);
                            DaughtersOlderThanParents=DaughtersOlderThanParents+1;
                        end
                        schnitzcells(child).P=a;
                    end
                end
            end
            if ~isempty(schnitzcells(b).D)
                if schnitzcells(b).D~=0
                    child=schnitzcells(b).D;
                    if schnitzcells(child).P==b;
                        if child>b
                            % display(['Downs everywhere!']);
                            DaughtersOlderThanParents=DaughtersOlderThanParents+1;
                        end
                        schnitzcells(child).P=a;
                    end
                end
            end
            schnitz(a)=swap;
            if ~isempty(schnitzcells(a).P)
                if    schnitzcells(a).P~=0
                    parent=schnitzcells(a).P;
                    if schnitzcells(parent).E==a;
                        schnitzcells(parent).E=b;
                    elseif schnitzcells(parent).D==a;
                        schnitzcells(parent).D=b;
                    end
                end
            end
            if ~isempty(schnitzcells(a).E)
                if schnitzcells(a).E~=0
                    child=schnitzcells(a).E;
                    if schnitzcells(child).P==a;
                        if child>b
                            % display(['Downs everywhere!']);
                            DaughtersOlderThanParents=DaughtersOlderThanParents+1;
                        end
                        schnitzcells(child).P=b;
                    end
                end
            end
            if ~isempty(schnitzcells(a).D)
                if schnitzcells(a).D~=0
                    child=schnitzcells(a).D;
                    if schnitzcells(child).P==a;
                        if child>b
                            % display(['Downs everywhere!']);
                            DaughtersOlderThanParents=DaughtersOlderThanParents+1;
                        end
                        schnitzcells(child).P=b;
                    end
                end
            end
            
            % Particles
            if(cpexist)
                ParticleIndicesForSwap= Particles.Nucleus==a;
                Particles(ParticleIndicesForSwap).Nucleus=b;
                ParticleIndicesForSwap= Particles.Nucleus==b;
                Particles(ParticleIndicesForSwap).Nucleus=a;
            end
        end
        
    end
end