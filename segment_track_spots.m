function segment_track_spots(Prefix, thresh, show_status,save_status, num_frames)

%show_status takes 0 or 1 depending on whether you want to display plots
%and images

%save_status takes 0 or 1 depending on whether you want to save your files

%num_frames for debugging should be kept at 5-20

%thresh should be kept at ~90-200

% try
%     parpool;
% catch
% end

%%    
%Get the relevant folders now:

tic;
[SourcePath,FISHPath,DropboxFolder,MS2CodePath,PreProcPath]=...
    DetermineLocalFolders(Prefix);

%Figure out what type of experiment we have
[XLSNum,XLSTxt]=xlsread([DropboxFolder,filesep,'MovieDatabase.xlsx']);
DataFolderColumn=find(strcmp(XLSTxt(1,:),'DataFolder'));
ExperimentTypeColumn=find(strcmp(XLSTxt(1,:),'ExperimentType'));
Channel1Column=find(strcmp(XLSTxt(1,:),'Channel1'));
Channel2Column=find(strcmp(XLSTxt(1,:),'Channel2'));

% Convert the prefix into the string used in the XLS file
Dashes = strfind(Prefix, '-');
PrefixRow = find(strcmp(XLSTxt(:, DataFolderColumn),...
    [Prefix(1:Dashes(3)-1), '\', Prefix(Dashes(3)+1:end)]));
if isempty(PrefixRow)
    PrefixRow = find(strcmp(XLSTxt(:, DataFolderColumn),...
        [Prefix(1:Dashes(3)-1), '/', Prefix(Dashes(3)+1:end)]));
    if isempty(PrefixRow)
        error('Could not find data set in MovieDatabase.XLSX. Check if it is defined there.')
    end
end

ExperimentType=XLSTxt(PrefixRow,ExperimentTypeColumn);

if strcmp(ExperimentType,'1spot')
    NChannels=1;
elseif strcmp(ExperimentType,'2spot')
    NChannels=1;
elseif strcmp(ExperimentType,'2spot1color')
    NChannels=1;
elseif strcmp(ExperimentType,'2spot2color')
    NChannels=2;
elseif strcmp(ExperimentType,'inputoutput')
    NChannels=1;
else
    error('Experiment type not recognized in MovieDatabase.XLSX')
end

if ~exist('Thresholds')
    Thresholds=ones(1,NChannels)*inf;
else
    if length(Thresholds)~=NChannels
        error('Number of channels in movie does not match number of thresholds input')
    end
end
HyphenPositionR = find(Prefix == '-');
DateFolderS = Prefix(1 : HyphenPositionR(3)-1);
LineFolderS = Prefix(HyphenPositionR(3)+1 : end);
Folder = [SourcePath, filesep, DateFolderS, filesep, LineFolderS];
%Determine whether we're dealing with 2-photon data from Princeton or LSM
%data. 2-photon data uses TIF files. In LSM mode multiple files will be
%combined into one.
DTIF=dir([Folder,filesep,'*.tif']);
DLSM=dir([Folder,filesep,'*.lsm']);
DLIF=dir([Folder,filesep,'*.lif']);
DLAT=dir([Folder,filesep,'*_Settings.txt']);

if (length(DTIF)>0)&(isempty(DLSM))
    if length(DLIF)==0
        if length(DLAT)==0
            display('2-photon @ Princeton data mode')
            FileMode='TIF';
        else
            display('Lattice Light Sheet data mode')
            FileMode='LAT';
        end
    else
        display('LIF export mode')
        FileMode='LIFExport';
    end
elseif (length(DTIF)==0)&(length(DLSM)>0)
    display('LSM mode')
    FileMode='LSM';
else
    error('File type not recognized. For LIF files, were they exported to TIF?')
end

%%
%Here I want to load each of the tifs and iterate through them, generating
%dog images in each iteration

% Load data generated with ExportDataForFISH


if strcmp(FileMode,'LIFExport')
    load([DropboxFolder,filesep,Prefix,filesep,'FrameInfo.mat']);
    im_stack = {};
    if num_frames == 0
        num_frames = length(FrameInfo);
    end
    
    for current_frame = 1:num_frames
        for current_z = 1:FrameInfo(1).NumberSlices
            im_stack{current_frame,current_z} = imread([PreProcPath,filesep,Prefix,filesep,Prefix,'_',iIndex(current_frame,3),'_z',iIndex(current_z,2),'.tif']);
        end
    end

elseif strcmp(FileMode,'LAT')
    im_stack = {};
    
    if num_frames == 0
        num_frames = length(DTIF);
    end
    
    for j = 1:num_frames-1 %For a full analysis, i'll run this to length(DTIF)
        fname = [Folder, filesep, DTIF(j).name];
        info = imfinfo(fname);
        num_images = numel(info);
        for i = 1:num_images
            im_stack{j, i} = imread(fname, i, 'Info', info);
        end
    end
end

OutputFolder1=[FISHPath,filesep,Prefix,'_',filesep,'MYCODEdogsMYCODE'];
OutputFolder2=[FISHPath,filesep,Prefix,'_',filesep,'MYCODEsegsMYCODE'];

mkdir(OutputFolder1)
mkdir(OutputFolder2)
%%
%DoG Stuff
%filterSize >> sigma 2 > sigma 1. these values should be good for a first
%pass.
pixelSize = 100; %nm
% pixelSize = 200;
sigma1 = 150 / pixelSize;
sigma2 = 250 / pixelSize;
filterSize = 1500 /pixelSize; 
neighb = 500 / pixelSize; %This should work for a first pass and shouldn't fail on sisters.
thr = thresh; 
dog_stack  = {}; 
all_frames = {}; 
for current_frame = 1:num_frames-1 
    for current_z = 1:size(im_stack,2) %z-slices
        im = im_stack{current_frame,current_z};
        %filterSize >> sigma 2 > sigma 1. these values should be good for a first pass.
        dog_stack{current_frame,current_z} = conv2(single(im), single(fspecial('gaussian',filterSize, sigma1) - fspecial('gaussian',filterSize, sigma2)),'same');
        if save_status
            dog_name = ['DOG_',Prefix,'_',iIndex(current_frame,3),'_z',iIndex(current_z,2),'.tif'];
            imwrite(uint16(dog_stack{current_frame,current_z}), [OutputFolder1,filesep,dog_name])
        end
        dog = dog_stack{current_frame,current_z}(10:end-10, 10:end-10);
        im = im(10:end-10, 10:end-10);
        if show_status
            figure(1)
            imshow(im,[]);
        end
        thrim = dog > thr;
        [im_label, n_spots] = bwlabel(thrim); 
        temp_particles = {};
        rad = 600 / pixelSize; %500nm is roughly the size of a sister chromatid diffraction limited spot.
%         rad = 500/pixelSize;
        for k = 1:n_spots
            [r,c] = find(im_label == k);

            %Find spot centroids in the actual image by hunting for absolute maxima in
            %neighborhoods around spots that were just located

            possible_cent = [];
            pcentloc = {};
            cent = [];
            cent_intensity = 0;
            for o = 1:2*neighb
                for p = 1:2*neighb
                    if r(1) - neighb + o > 0 && c(1) - neighb + p > 0 ... 
                            && r(1) - neighb + o < size(im,1)  && c(1) - neighb + p < size(im,2)
                        possible_cent(o,p) = im(r(1)-neighb+o, c(1)-neighb+p);
                        pcentloc{o,p} = [r(1)-neighb+o, c(1)-neighb+p];
                    end
                end
            end
            if ~isempty(possible_cent)
                [inten, index] = max(possible_cent(:));
                [row, col] = ind2sub(size(possible_cent),index);
                cent_y = pcentloc{row,col}(1); 
                cent_x = pcentloc{row,col}(2);
               % temp_particles = [temp_particles,[0, cent_x, cent_y, 0, 0]];
               temp_particles = {};
               if show_status
                    ellipse(neighb/2,neighb/2,0,cent_y,cent_x,'r');
               end

                if cent_y - rad > 1 && cent_x - rad > 1 && cent_y + rad < size(im, 1) && cent_x + rad < size(im,2)
                    snip = im(cent_y-rad:cent_y+rad, cent_x-rad:cent_x+rad);
%                      [f1, res1, f2, res2] = fitGausses(snip);

                    % Set parameters to use as initial guess in the
                    % fitting. For the lattice data, try NeighborhoodSize =
                    % [9 9], MaxThreshold = 2000, WidthGuess = 5,
                    % OFfsetGuess = 1000.
                    
                    NeighborhoodSize = [5 5];
                    MaxThreshold = 30;
                    WidthGuess = 1;
                    OffsetGuess = 10;
                    [f1, res1, residual, exitflag, output, lambda, jacobian] =  ...
                        fitTwoGausses(snip, NeighborhoodSize, MaxThreshold, ...
                        WidthGuess, OffsetGuess, show_status);
                    
                    ci = nlparci(f1,residual,'jacobian',jacobian);
                    errors = zeros(1, length(f1));
                    for ndx = 1:length(ci)
                        errors(ndx) = abs((abs(ci(ndx, 1)) - abs(ci(ndx, 2)))/2);
                    end
                    rel_errors = abs(errors./f1);
                    disp(rel_errors);
                    
                    % Quality control.
                    % TODO: make some quality control using the errors in
                    % the fitting. Using any(rel_errors > 0.3) (for
                    % example) is not the best thing to do, because
                    % sometimes the second gaussian doesn't get a good fit
                    % but the first one does, and the second one is good
                    % enough to position its center.
                    
%                     if f1(3) > rad+3 || f1(5) > rad+3
                    if 1
                        c_x = f1(2) - rad + cent_x;
                        c_y = f1(4) - rad + cent_y;
                        int_x = [round(c_x - f1(3)), round(c_x + f1(3))];
                        int_y = [round(c_y - f1(5)), round(c_y + f1(5))];
                        area = (int_x(2) - int_x(1)) * (int_y(2) - int_y(1));
                        fluor = 0;
                        if int_x(1) > 1 && int_y(1) > 1 && int_x(2) < size(im,2) && int_y(2) < size(im,1)
                            for w = int_x(1):int_x(2)
                                for v = int_y(1): int_y(2)
                                    fluor = fluor + double(im(v,w));
                                end
                            end
                            temp = {{fluor, c_x, c_y, f1(6), snip, area}};
                            temp_particles = [temp_particles,temp];
                        end
                    end
                end
            end
            if k == n_spots && save_status
                seg_name = ['SEG_',Prefix,'_',iIndex(current_frame,3),'_z',iIndex(current_z,2),'.tif'];
                saveas(gcf,[OutputFolder2,filesep,seg_name]);
            end
        end
        all_frames{current_frame,current_z} = temp_particles;
    end
end

%%

clear im_stack
clear dog_stack
n = 1;
nframes = size(all_frames,1);
nz = size(all_frames,2);
for i = 1:nframes 
    for j = 1:nz 
         for spot = 1:length(all_frames{i,j}) %spots within particular image
             if ~isempty(all_frames{i,j}{spot})
                 Particles(n).Intensity(1) = cell2mat(all_frames{i,j}{spot}(1));
                 Particles(n).x(1) = cell2mat(all_frames{i,j}{spot}(2));
                 Particles(n).y(1) = cell2mat(all_frames{i,j}{spot}(3));
                 Particles(n).Offset(1) = cell2mat(all_frames{i,j}{spot}(4));
%                  Particles(n).Sister(1) = all_frames{i,j}{spot}(5);
                 Particles(n).Snippet{1} = cell2mat(all_frames{i,j}{spot}(5));
                 Particles(n).Area{1} = cell2mat(all_frames{i,j}{spot}(6));
                 Particles(n).z(1) = j;
                 Particles(n).t(1) = i;
                 Particles(n).r = 0;
                 n = n + 1;
             end
         end
    end
end

fields = fieldnames(Particles);

% z tracking
changes = 1;
while changes ~= 0
    changes = 0;
    i = 1;
    for n = 1:nframes   
        i = i + length(Particles([Particles.t] == (n - 1) ));
        for j = i:i+length(Particles([Particles.t] == n)) - 1
            for k = j+1:i+length(Particles([Particles.t] == n)) - 1
                dist = sqrt( (Particles(j).x(end) - Particles(k).x(end))^2 + (Particles(j).y(end) - Particles(k).y(end))^2); 
                if dist < neighb && Particles(j).z(end) ~= Particles(k).z(end)
                    for m = 1:numel(fields)-2 %do not include fields 'r' or 't'
                        Particles(j).(fields{m}) = [Particles(j).(fields{m}), Particles(k).(fields{m})];
                    end
                    Particles(k).r = 1;
                    changes = changes + 1;
                end
            end
        end
    end
    Particles = Particles([Particles.r]~=1);
end


%pick the brightest z-slice
for i = 1:length(Particles)
    [~, max_index] = max(Particles(i).Intensity);
    for j = 1:numel(fields)-2 %do not include fields 'r' or 't'
        Particles(i).(fields{j}) = Particles(i).(fields{j})(max_index);
    end
end

% %time tracking
% changes = 1;
% while changes ~= 0
%     changes = 0;
%     for n = 1:length(Particles)-1 %particle of interest
%         for j = n+1:length(Particles) %particle to compare to
%             if Particles(n).t(end) == (Particles(j).t(end) -  1)
%                 dist = sqrt( (Particles(n).x(end) - Particles(j).x(end))^2 + (Particles(n).y(end) - Particles(j).y(end))^2); 
%                 if dist < neighb
%                     for m = 1:numel(fields)-1 %do not include fields 'r'
%                         Particles(n).(fields{m}) = [Particles(n).(fields{m}), Particles(j).(fields{m})];
%                     end
%                     Particles(j).r = 1;
%                     changes = changes + 1;
%                 end
%             end
%         end
%     end
% Particles = Particles([Particles.r]~=1);
% end

mkdir([DropboxFolder,filesep,Prefix,filesep,'mycode']);
save([DropboxFolder,filesep,Prefix,filesep,'mycode',filesep,'Particles.mat'], 'Particles');
for i = 1:length(Particles)
    if length(Particles(i).t) > 50
        plot(Particles(i).t, Particles(i).Intensity)
        i
        hold on
    end
end
t = toc
end