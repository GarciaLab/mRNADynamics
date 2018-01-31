% Script to Build Classifier for Sister Detection 

Prefix = '2017-07-10-eve2_20sec_20';
DropboxFolder = 'D:\Data\Nick\LivemRNA\LivemRNAFISH\Dropbox (Garcia Lab)\mHMM\weka';
FISHPath = 'D:\Data\Nick\LivemRNA\LivemRNAFISH\';
%Load the data
load([DropboxFolder,filesep,Prefix,'\Spots.mat'])
load([DropboxFolder,filesep,Prefix,filesep,'Particles.mat'])
load([DropboxFolder,filesep,Prefix,'\Ellipses.mat'])
%%
%%%Project Path 
data_path = 'D:\Data\Nick\projects\hmmm\dat\mHMMeve2_weka_inf_2018_01_23\';
sister_name = [data_path 'sister_training_data.mat'];
sister_file_flag = 0;
if exist(sister_name) == 2
    sister_file_flag = 1;
end
%%% Analysis Parameters
snippet_size = 21;
train = 1;
classify = 0;
% ---------- Load Snippets on the fly to generate classifications --------%
if train
    particle_frames = [Particles.Frame];
    particle_nuclei = [];
    for i = 1:length(Particles)
        particle_nuclei = [particle_nuclei repelem(Particles(i).Nucleus,length(Particles(i).Frame))];
    end    
    particle_indices = [Particles.Index];
    particle_id = [];
    for i = 1:length(Particles)
        particle_id = [particle_id repelem(i,length(Particles(i).Frame))];
    end
    fn = fieldnames(Particles);
    if sister_file_flag
        load(sister_name);
        PrefixList = sister_struct.PrefixList;
        IDvec = sister_struct.IDvec;
        SetID = find(strcmp(Prefix,PrefixList)); % check to see if project is already in struct
        particle_labels = sister_struct.particle_labels;
        if isempty(SetID)
            disp('Current set not found in PefixList. Adding')
            sister_struct.PrefixList = [PrefixList{:} {Prefix}];
            SetID = length(PrefixList);
            sister_struct.IDvec = [IDvec SetID];            
        end
        labeled_nc = particle_labels([particle_labels.SetID]==SetID).Nucleus;
        labeled_frames = particle_labels([particle_labels.SetID]==SetID).Frame;
    else % if no pre-existing file, generate one
        disp('No sister data set detected. New one will be generated')
        sister_struct.IDvec = [1];
        sister_struct.PrefixList = {Prefix};
        particle_labels = [];
        labeled_nc = [];
        labeled_frames = [];
    end
    
    % limit search to unlabeled spots
    label_vec = ~(ismember(particle_nuclei,labeled_nc)&ismember(particle_frames,labeled_frames));
    selection_vec = find(label_vec);    
    % randomize options
    shuffled_vec = randsample(selection_vec,length(selection_vec),false);
    exit_flag = 0;
    iter = 1;
    cm = jet(64);
    for i = shuffled_vec
        index = shuffled_vec(i);
        % create temp struct
        temp = struct;
        temp.SetID = SetID;
        temp.selection_index = index;                
        Frame = particle_frames(index);
        temp.Frame = Frame;
        ParticleIndex = index;
        temp.ParticleIndex = ParticleIndex;
        SpotIndex = particle_indices(index);
        temp.SpotIndex = SpotIndex;        
        Nucleus = particle_nuclei(index);
        temp.Nucleus = Nucleus;
        % get spot location
        xSpot = Spots(Frame).Fits(SpotIndex).xDoG;
        ySpot = Spots(Frame).Fits(SpotIndex).yDoG;
        zSpot = Spots(Frame).Fits(SpotIndex).z;
        temp.x = xSpot; 
        temp.y = ySpot;
        temp.z = zSpot;
        bz = Spots(Frame).Fits(SpotIndex).brightestZ;
        temp.brightestZ = bz;
        zSpotGrid = min(zSpot)-1:max(zSpot)+1;
        % load appropriate image stack
        D=dir([FISHPath,'\Data\PreProcessedData\',Prefix,filesep,Prefix,'_',iIndex(Frame,3),'_z*.tif']);
        % Do not load the first and last frame as they are black    
        ImageTemp=[];    
        for j=2:(length(D)-1)        
            ImageTemp(:,:,j-1)=imread([FISHPath,'\Data\PreProcessedData\',Prefix,filesep,D(j).name]);                    
        end
        sd = floor(snippet_size/2);
        mRNAImage = ImageTemp(max(1,ySpot(zSpot==bz)-sd):min(256,ySpot(zSpot==bz)+sd),...
                    max(1,xSpot(zSpot==bz)-sd):min(512,xSpot(zSpot==bz)+sd),:);
        MaxMCP = max(mRNAImage,[],3);
        HistoneImage=imread([FISHPath,'\Data\PreProcessedData\',Prefix,filesep,...
            Prefix,'-His_',iIndex(Frame,3),'.tif']);
        HistoneImage = HistoneImage(max(1,ySpot(zSpot==bz)-sd):min(256,ySpot(zSpot==bz)+sd),...
                    max(1,xSpot(zSpot==bz)-sd):min(512,xSpot(zSpot==bz)+sd));
        ImOverlay=cat(3,mat2gray(HistoneImage),mat2gray(MaxMCP),zeros(size(MaxMCP)));
    %     mRNAImage=max(ImageTemp,[],3); % Max project 
        nc = bz; % initially look at brightest frame
        cc = '';
        while ~strcmp(cc,'0')&&~strcmp(cc,'1')&&~strcmp(cc,'2')&&~strcmp(cc,'3')&&~strcmp(cc,'4')        
            HistFig = figure;
            hold on
            imshow(ImOverlay,'InitialMagnification','fit');                        
            title('Histone Channel')
            HistFig.Position = [950,350,256,256];
            axis([1 snippet_size 1 21])

            ZYFig = figure;
            hold on
            imshow(mat2gray(reshape(max(mRNAImage,[],2),snippet_size,21)'),'InitialMagnification','fit');
            plot(1:snippet_size,repelem(nc-1,snippet_size),'Color',cm(35,:));
            axis on
            ylabel('z slice')
            title('Z Profile (max X projection)')
            ZYFig.Position = [700,350,256,256];
            axis([1 snippet_size 1 21])

            ZXFig = figure;
            hold on
            imshow(mat2gray(reshape(max(mRNAImage,[],1),snippet_size,21)'),'InitialMagnification','fit');
            plot(1:snippet_size,repelem(nc-1,snippet_size),'Color',cm(35,:));
            axis on
            ylabel('z slice')
            title('Z Profile (max Y projection)')
            ZXFig.Position = [700,50,256,256];
            axis([1 snippet_size 1 21])
            SpotFig = figure('Position',[0 0 512 512]);
            slice = repelem(mRNAImage(:,:,nc-1),25,25); 
            imshow(mat2gray(slice));
            title(['Particle: ' num2str(ParticleID) ' (' num2str(iter) ...
                ' of ' num2str(length(shuffled_vec)) ') Frame: ' num2str(Frame) ...
                    ' Z Slice: ' num2str(nc)])

            ct=waitforbuttonpress;
            cc=get(SpotFig,'currentcharacter');
            % 4 unusable (QC); 3 interference from neighboring nucleus
            % 2 ambiguous; 1 resolved sisters; 0 no resolved sisters
            if strcmp(cc,'3')||strcmp(cc,'2')||strcmp(cc,'1')||strcmp(cc,'0')
                temp.SisterClass = eval(cc);
                particle_labels = [particle_labels temp];                
            elseif strcmp(cc,'x')
                exit_flag = 1;
                break
            elseif strcmp(cc,'a')
                nc = min(21,nc-1)+2;
            elseif strcmp(cc,'z')
                nc = max(1,nc-2)+1;            
            end       
        end 
        close all
        if exit_flag
            disp('Exiting')
            break
        end
        iter = iter + 1;
    end
    % add labels to struct
    sister_struct.particle_labels = particle_labels;
end    
save(sister_name,'sister_struct')
% %% Movie of trace duration
% 
% [px, py] = meshgrid(1:512,1:256); % reference coordinate grid
% Frames = CompiledParticles(trace_id).Frame; % get frames for relevant particle
% FrameRange = min(CompiledParticles(trace_id).Frame):max(CompiledParticles(trace_id).Frame);
% Time = ElapsedTime(FrameRange); 
% Fluo = CompiledParticles(trace_id).Fluo;
% x_pt = CompiledParticles(trace_id).xPos; % particle x position
% y_pt = CompiledParticles(trace_id).yPos; % particle y position
% FluoInterp = interp1(ElapsedTime(Frames),Fluo,Time); % Interpolate to fill gaps
% x_pt = interp1(ElapsedTime(Frames),x_pt,Time); % Interpolate to fill gaps
% y_pt = interp1(ElapsedTime(Frames),y_pt,Time); % Interpolate to fill gaps
% start = min(FrameRange);
% Time = Time - ElapsedTime(nc14);
% % Display Params
% x_size = 20; % x dimension in pixels
% y_size = 20; % y dimension in pixels
% zoomFactor = 1; % Factor by which to expand image
% radius = 4;
% 
% iter = 0;
% for CurrentFrame = FrameRange 
%     iter = iter + 1;
%     % Get Nucleus Info    
%     CurrentEllipse=...
%         schnitzcells(CompiledParticles(trace_id).Nucleus).cellno(...
%         schnitzcells(CompiledParticles(trace_id).Nucleus).frames==...
%         CurrentFrame);
%     x_nc =  Ellipses{CurrentFrame}(CurrentEllipse,1)+1; % x position of nucleus
%     y_nc =  Ellipses{CurrentFrame}(CurrentEllipse,2)+1; % y position of nucleus
%     % Draw a Circle Around Spot
%     circle_mat = zeros(256,512);
%     xp = x_pt(iter);
%     yp = y_pt(iter);
%     y_dist = abs(py - yp);
%     x_dist = abs(px - xp);
%     r_mat = round(sqrt(x_dist.^2 + y_dist.^2));
%     circle_mat(r_mat==radius) = .7; 
% %     circle_mat(r_mat<radius) = .2; 
%     %Make a maximum projection of the mRNA channel
%     D=dir([FISHPath,'\Data\PreProcessedData\',Prefix,filesep,Prefix,'_',iIndex(CurrentFrame,3),'_z*.tif']);
%     % Do not load the first and last frame as they are black
%     ImageTemp=[];
%     for i=2:(length(D)-1)
%         ImageTemp(:,:,i-1)=imread([FISHPath,'\Data\PreProcessedData\',Prefix,filesep,D(i).name]);
%     end
%     mRNAImage=max(ImageTemp,[],3); % Max project    
%     mRNAImage = mat2gray(mRNAImage);
%     %Load the corresponding histone image
%     HistoneImage=imread([FISHPath,'\Data\PreProcessedData\',Prefix,filesep,...
%         Prefix,'-His_',iIndex(CurrentFrame,3),'.tif']);        
%     
%     %Overlay all channels
%     MCPChannel=mat2gray(mRNAImage);
%     HistoneChannel=mat2gray(HistoneImage);    
%     BlueChannel  = mat2gray(circle_mat,[0,1]);
%     ImOverlay=cat(3,HistoneChannel,MCPChannel,BlueChannel);
% 
%     ImCropped = ImOverlay(max(1,y_nc-y_size):min(256,y_nc+y_size),...
%                 max(1,x_nc-x_size):min(512,x_nc+x_size),:);
%     ImCropped = repelem(ImCropped,zoomFactor,zoomFactor); % Blow it up a bit...     
%     OverlayFig = figure(1);
%     OverlayFig.Position = [0 0 1024 512];%('Visible','off');
%     OverlayFig.Visible = 'off';
%     clf
%     subplot(1,2,1);
%     imshow(ImCropped)   
%     
%     text(1,2*y_size*zoomFactor,[iIndex(round(ElapsedTime(CurrentFrame)-ElapsedTime(nc14)),2),...
%         ' min'],'Color','k','FontSize',10,'BackgroundColor',[1,1,1,.5])    
%     
%     subplot(1,2,2);
%     hold on
%     plot(Time(1:(CurrentFrame-start+1)),FluoInterp(1:(CurrentFrame-start+1)),'LineWidth',2,'Color','black');
%     scatter(Time(CurrentFrame-start+1),FluoInterp(CurrentFrame-start+1),50,'black','filled')
%     axis([0 1.2*Time(end) 0 1.2*max(FluoInterp)])
%     xlabel('time (min)')
%     ylabel('AU')
% %     pause(.05)
%     zs = '000';
%     f_string = [zs(1:3-length(num2str(CurrentFrame))) num2str(CurrentFrame)];
%     saveas(figure(1),[OutPath 'Trace_' num2str(trace_id) '_Frame_' ...
%         num2str(CurrentFrame) '.tiff'],'tiff');
% end 