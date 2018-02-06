% Script to Build Classifier for Sister Detection 

Prefix = '2017-07-10-eve2_20sec_20';
DropboxFolder = 'D:\Data\Nick\LivemRNA\LivemRNAFISH\Dropbox (Garcia Lab)\mHMM\weka';
FISHPath = 'D:\Data\Nick\LivemRNA\LivemRNAFISH\';
%Load the data
load([DropboxFolder,filesep,Prefix,'\Spots.mat'])
load([DropboxFolder,filesep,Prefix,filesep,'Particles.mat'])
load([DropboxFolder,filesep,Prefix,'\Ellipses.mat'])
load([DropboxFolder,filesep,Prefix,'\FrameInfo.mat'])
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
pixel_size = FrameInfo.PixelSize;
z_step_size = FrameInfo.ZStep;
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
        setID = find(strcmp(Prefix,PrefixList)); % check to see if project is already in struct
        particle_labels = sister_struct.particle_labels;
        if isempty(setID)
            disp('Current set not found in PefixList. Adding')
            sister_struct.PrefixList = [PrefixList{:} {Prefix}];
            setID = length(PrefixList);
            sister_struct.IDvec = [IDvec setID];            
        end
        labeled_nc = particle_labels([particle_labels.setID]==setID).Nucleus;
        labeled_frames = particle_labels([particle_labels.setID]==setID).Frame;
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
        si = strfind(Prefix,'_');
        setID = Prefix(si(end)+1:end);
        % create temp struct
        temp = struct;
        temp.setID = setID;
        temp.selection_index = index;                
        Frame = particle_frames(index);
        ParticleID = particle_id(index);
        temp.Frame = Frame;
        temp.ParticleID = ParticleID;
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

%% Use Classified Images to Train Simple Classifier
load(sister_name)
particle_labels = sister_struct.particle_labels;
% generate relevant statistics for each particle
for i = 1:length(particle_labels)
    Frame = particle_labels(i).Frame;
    brightestZ = particle_labels(i).brightestZ;
    xSpot = particle_labels(i).x;
    ySpot = particle_labels(i).y;
    zSpot = particle_labels(i).z;    
    % first load image
    D=dir([FISHPath,'\Data\PreProcessedData\',Prefix,filesep,Prefix,'_',iIndex(Frame,3),'_z*.tif']);
    % Do not load the first and last frame as they are black    
    ImageTemp=[];    
    for j=2:(length(D)-1)        
        ImageTemp(:,:,j-1)=imread([FISHPath,'\Data\PreProcessedData\',Prefix,filesep,D(j).name]);                    
    end
    [X1,X2,X3] = meshgrid(1:snippet_size,1:snippet_size,1:snippet_size);

    for n = 1:n_sim
        snippet = zeros(snippet_size,snippet_size,snippet_size);
        n_spots = 2 - round(rand());
        for s = 1:n_spots       
            mu = [randsample(mean_r:snippet_size-mean_r,1),randsample(mean_r:snippet_size-mean_r,1), randsample(mean_r:snippet_size-mean_r,1)]; 
            xr = mean_r - sigma_r*norminv(rand()*lower_r);  
            yr = mean_r - sigma_r*norminv(rand()*lower_r); 
            zr = mean_r - sigma_r*norminv(rand()*lower_r); 
            F = mvnpdf([X1(:) X2(:) X3(:)],mu,[xr,yr,zr]);
            F = reshape(F,snippet_size,snippet_size,snippet_size);
            f = mean_f - sigma_f*norminv(rand()*lower_f);        
            snippet = snippet + (1/n_spots)*f*F;        
        end
        x_center  = sum((1:snippet_size).*sum(sum(snippet),3))/sum(snippet(:));
        y_center  = sum((1:snippet_size).*reshape(sum(sum(snippet,2),3),1,[]))/sum(snippet(:));
        z_center  = sum((1:snippet_size).*reshape(sum(sum(snippet,1),2),1,[]))/sum(snippet(:));
    %     r_sq_mat = (X1-x_center).^2 + (X2 - y_center).^2;
        r_sq_mat = (X1-x_center).^2 + (X2 - y_center).^2 + (X3 - z_center).^2;
        moment2 = [moment2 sum(sum(sum(r_sq_mat.*snippet)))/sum(snippet(:))];
        com_register = [com_register {[x_center,y_center]}];
        snippet_register{n} = snippet;
        spot_counts = [spot_counts n_spots];
    end

    mc = jet(64);

    metric_fig = figure;
    hold on
    g_size = ceil(max(moment2)/100);
    grid = 0:g_size:ceil(max(moment2));
    s1 = histc(moment2(spot_counts==1),grid);
    s1 = s1/sum(s1);
    s2 = histc(moment2(spot_counts==2),grid);
    s2 = s2/sum(s2);
    bar(grid,s1,'FaceColor',cm(15,:),'BarWidth',1,'FaceAlpha',.5,'EdgeAlpha',0);
    bar(grid,s2,'FaceColor',cm(30,:),'BarWidth',1,'FaceAlpha',.5,'EdgeAlpha',0);

end
    
    