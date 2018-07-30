function Values=CreateMovies(Prefix)

%The idea here is to create different types of movies

%% Determine folders and load the data

% Prefix = '2012-06-26-MCP(10)TM3-X1_v2';     %This was posterior in the paper
% Prefix = '2012-06-26-MCP(10)TM3-X1';       %This was anterior in the paper

[SourcePath,FISHPath,DefaultDropboxFolder,MS2CodePath,SchnitzcellsFolder]=...
    DetermineLocalFolders;
%Now get the dropbox folders for the different data sets.
[SourcePath,FISHPath,DropboxFolder,MS2CodePath,SchnitzcellsFolder]=...
    DetermineLocalFolders(Prefix);



%Load the data
load([DropboxFolder,filesep,Prefix,'\CompiledParticles.mat'])
load([DropboxFolder,filesep,Prefix,filesep,Prefix,'_lin.mat'])
load([DropboxFolder,filesep,Prefix,'\Ellipses.mat'])


%We need to figure out whether this was done with the Elowitz or Laurent
%versions of the tracking code. We do this by looking at Ellipses.

if ~(exist([FISHPath,filesep,'Analysis',filesep,Prefix,filesep,'dataStructure.mat']))|...
    ~(exist([FISHPath,filesep,'Analysis',filesep,Prefix,filesep,'TrackingDataStructure.mat']))
    TrackingVersion=0;      %Elowitz tracking
    
    %We'll remove one frame from all schnitzcell elements
    for i=1:length(schnitzcells)
        schnitzcells(i).frames=schnitzcells(i).frames-1;
    end
else
    TrackingVersion=1;      %Laurent tracking
end




%% Movie of probability of turning on


%Make folder for the movie output
mkdir([DropboxFolder,filesep,Prefix,filesep,'OnMovie']);
mkdir([DropboxFolder,filesep,Prefix,filesep,'OnMovie\MCPChannel']);
mkdir([DropboxFolder,filesep,Prefix,filesep,'OnMovie\HisChannel']);
mkdir([DropboxFolder,filesep,Prefix,filesep,'OnMovie\MCP+HisChannel']);
mkdir([DropboxFolder,filesep,Prefix,filesep,'OnMovie\Ellipses']);




CurrentNC=14;



if CurrentNC==14
    FrameRange=nc14:length(ElapsedTime);
elseif CurrentNC==13
    FrameRange=nc13:nc14;
else
    error('Only for nc13 and nc14');
end


%Choose which channels to output
OutputHistone=0;
OutputMCP=0;
OutputOverlay=1;
OutputEllipses=0;



%Add nuclei as detected
%This is something I had to do for one particular nucleus that had a
%particle, but no frame within it was quantifiable
if strcmp(Prefix,'2012-06-26-MCP(10)TM3-X1_v2')
    OriginalParticles=load([DropboxFolder,filesep,Prefix,'\Particles.mat']);
    OriginalParticles.Particles(253)
    
    CompiledParticles(end+1).OriginalParticle=253;
    CompiledParticles(end).Frame=OriginalParticles.Particles(253).Frame;
    CompiledParticles(end).Index=OriginalParticles.Particles(253).Index;
    CompiledParticles(end).Approved=1;
    CompiledParticles(end).xPos=OriginalParticles.Particles(253).xPos;
    CompiledParticles(end).yPos=OriginalParticles.Particles(253).yPos;
    CompiledParticles(end).APpos=OriginalParticles.Particles(253).APpos;
    CompiledParticles(end).APPos=OriginalParticles.Particles(253).APpos;
    CompiledParticles(end).FirstFrame=OriginalParticles.Particles(253).Frame(1);
    CompiledParticles(end).nc=14;
    CompiledParticles(end).Nucleus=OriginalParticles.Particles(253).Nucleus;
end
    
    




Values=[];
for CurrentFrame=FrameRange

    
    %Make a maximum projection of the mRNA channel
    D=dir([FISHPath,'\Data\',Prefix,filesep,Prefix,'_',iIndex(CurrentFrame,3),'_z*.tif']);
    %Do not load the first and last frame as they are black
    ImageTemp=[];
    for i=2:(length(D)-1)
        ImageTemp(:,:,i-1)=imread([FISHPath,'\Data\',Prefix,filesep,D(i).name]);
    end
    mRNAImage=max(ImageTemp,[],3);

  
    
    %Load the corresponding histone image
    HistoneImage=imread([FISHPath,'\Data\',Prefix,filesep,...
        Prefix,'-His_',iIndex(CurrentFrame,3),'.tif']);

    
    if OutputHistone

        figure(1)
        clf

        %Show the histone image
        imshow(HistoneImage,[],'Border','Tight')


        %Plot the ellipses on top only if they have been detected

        %Which nuclei had already a particle detected by now?
        ParticleFilter=([CompiledParticles.FirstFrame]<=CurrentFrame)&...
            ([CompiledParticles.nc]==CurrentNC);
        ParticlesToShow=find(ParticleFilter);

    %     hold on       
    % 
    %     EllipseHandle=[];
    %     for i=1:length(ParticlesToShow)
    %         try
    %             CurrentEllipse=...
    %                 schnitzcells(CompiledParticles(ParticlesToShow(i)).Nucleus).cellno(...
    %                 schnitzcells(CompiledParticles(ParticlesToShow(i)).Nucleus).frames==...
    %                 CurrentFrame);
    %                 EllipseHandle(end+1)=notEllipse(Ellipses{CurrentFrame}(CurrentEllipse,3),...
    %                     Ellipses{CurrentFrame}(CurrentEllipse,4),...
    %                     Ellipses{CurrentFrame}(CurrentEllipse,5),...
    %                     Ellipses{CurrentFrame}(CurrentEllipse,1)+1,...
    %                     Ellipses{CurrentFrame}(CurrentEllipse,2)+1,'r',50);
    %         catch
    %             display(['Error in frame ',num2str(CurrentFrame)])
    %         end
    %     end
    %     hold off
    %     set(EllipseHandle,'LineWidth',5)

    %     text(375,45,[iIndex(round(ElapsedTime(CurrentFrame)-ElapsedTime(FrameRange(1))),2),...
    %         ' min'],'Color','k','FontSize',30,'BackgroundColor',[1,1,1])
        xlim([25,487])
        ylim([25,223])
        drawnow

        SavedFlag=0;
        while ~SavedFlag
            try
%                 saveas(gcf,[DropboxFolder,filesep,Prefix,filesep,'OnMovie\HisChannel\nc',num2str(CurrentNC),...
%                     '-',num2str(iIndex(CurrentFrame,3)),'.tif']); 
                hgexport(gcf,[DropboxFolder,filesep,Prefix,filesep,'OnMovie\HisChannel\nc',num2str(CurrentNC),...
                     '-',num2str(iIndex(CurrentFrame,3)),'.tif'],...
                        hgexport('factorystyle'), 'Format', 'tiff')
                SavedFlag=1;
            catch
                display('Error saving. Retrying...')
            end
        end
    end
    
    if OutputMCP

        %Show the MCP channel
        figure(2)
        clf

        Values=[Values;min(min(mRNAImage)),max(max(mRNAImage))];

        imshow(mRNAImage,[145,1801/2],'Border','Tight')

    %     text(375,45,[iIndex(round(ElapsedTime(CurrentFrame)-ElapsedTime(FrameRange(1))),2),...
    %         ' min'],'Color','k','FontSize',30,'BackgroundColor',[1,1,1])
        xlim([25,487])
        ylim([25,223])
        drawnow

        SavedFlag=0;
        while ~SavedFlag
            try
                saveas(gcf,[DropboxFolder,filesep,Prefix,filesep,'OnMovie\MCPChannel\nc',num2str(CurrentNC),...
                    '-',num2str(iIndex(CurrentFrame,3)),'.tif']);   
                SavedFlag=1;
            catch
                display('Error saving. Retrying...')
            end
        end
    end
    
    if OutputOverlay


        %Overlay all channels
        MCPChannel=mat2gray(mRNAImage,[145,1801/2]);
        HistoneChannel=mat2gray(HistoneImage);

        %ImOverlay=cat(3,HistoneChannel,MCPChannel,zeros(size(MCPChannel)));
        ImOverlay=cat(3,HistoneChannel,HistoneChannel,HistoneChannel);

        figure(3)
        clf
        imshow(ImOverlay)

        %Plot the ellipses on top only if they have been detected

        %Which nuclei had already a particle detected by now?
        ParticleFilter=([CompiledParticles.FirstFrame]<=CurrentFrame)&...
            ([CompiledParticles.nc]==CurrentNC);
        ParticlesToShow=find(ParticleFilter);

        hold on       
        EllipseHandle=[];
        for i=1:length(ParticlesToShow)
            try
                EllipseIncrease=1.2;
                CurrentEllipse=...
                    schnitzcells(CompiledParticles(ParticlesToShow(i)).Nucleus).cellno(...
                    schnitzcells(CompiledParticles(ParticlesToShow(i)).Nucleus).frames==...
                    CurrentFrame);
                    EllipseHandle(end+1)=notEllipse(Ellipses{CurrentFrame}(CurrentEllipse,3)*EllipseIncrease,...
                        Ellipses{CurrentFrame}(CurrentEllipse,4)*EllipseIncrease,...
                        Ellipses{CurrentFrame}(CurrentEllipse,5),...
                        Ellipses{CurrentFrame}(CurrentEllipse,1)+1,...
                        Ellipses{CurrentFrame}(CurrentEllipse,2)+1,'w',50, gca);
            catch
                display(['Error in frame ',num2str(CurrentFrame)])
            end
        end
        hold off
        set(EllipseHandle,'LineWidth',2)

    %     text(375,45,[iIndex(round(ElapsedTime(CurrentFrame)-ElapsedTime(FrameRange(1))),2),...
    %         ' min'],'Color','k','FontSize',30,'BackgroundColor',[1,1,1])
        xlim([25,487])
        ylim([25,223])
        drawnow

        SavedFlag=0;
        while ~SavedFlag
            try
                saveas(gcf,[DropboxFolder,filesep,Prefix,filesep,'OnMovie\MCP+HisChannel\nc',num2str(CurrentNC),...
                    '-',num2str(iIndex(CurrentFrame,3)),'.tif']);   
                SavedFlag=1;
            catch
                display('Error saving. Retrying...')
            end
        end
    end
    
    if OutputEllipses
    
        %Ellipse channel
        EmptyImage=zeros(size(mRNAImage));


        figure(4)
        clf
        imshow(EmptyImage)

        %Plot the ellipses on top only if they have been detected

        %Which nuclei had already a particle detected by now?
        ParticleFilter=([CompiledParticles.FirstFrame]<=CurrentFrame)&...
            ([CompiledParticles.nc]==CurrentNC);
        ParticlesToShow=find(ParticleFilter);

        hold on       
        EllipseHandle=[];
        for i=1:length(ParticlesToShow)
            try
                EllipseIncrease=1.2;
                CurrentEllipse=...
                    schnitzcells(CompiledParticles(ParticlesToShow(i)).Nucleus).cellno(...
                    schnitzcells(CompiledParticles(ParticlesToShow(i)).Nucleus).frames==...
                    CurrentFrame);
                    EllipseHandle(end+1)=notEllipse(Ellipses{CurrentFrame}(CurrentEllipse,3)*EllipseIncrease,...
                        Ellipses{CurrentFrame}(CurrentEllipse,4)*EllipseIncrease,...
                        Ellipses{CurrentFrame}(CurrentEllipse,5),...
                        Ellipses{CurrentFrame}(CurrentEllipse,1)+1,...
                        Ellipses{CurrentFrame}(CurrentEllipse,2)+1,'b',50, gca);
            catch
                display(['Error in frame ',num2str(CurrentFrame)])
            end
        end
        hold off
        set(EllipseHandle,'LineWidth',2)

    %     text(375,45,[iIndex(round(ElapsedTime(CurrentFrame)-ElapsedTime(FrameRange(1))),2),...
    %         ' min'],'Color','k','FontSize',30,'BackgroundColor',[1,1,1])
        xlim([25,487])
        ylim([25,223])
        drawnow

        SavedFlag=0;
        while ~SavedFlag
            try
                saveas(gcf,[DropboxFolder,filesep,Prefix,filesep,'OnMovie\Ellipses\nc',num2str(CurrentNC),...
                    '-',num2str(iIndex(CurrentFrame,3)),'.tif']);   
                SavedFlag=1;
            catch
                display('Error saving. Retrying...')
            end
        end
    end
    
end
    

    












