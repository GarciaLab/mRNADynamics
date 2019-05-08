function combineAndCleanTraces(project, varargin)
%
% DESCRIPTION
% Funcion to compile relevant outputs from image analysis pipeline across
% multiple experiments
%
%
% ARGUMENTS
% project: master ID variable 
%
% keyword: String contained in all folder names for projects one wishes to 
%          compile. For instance: 'Eve2MS2', will full all projects
%          containing this string
%
% OPTIONS
% include_vec:  Vector specifying IDs within list of projects
%               pulled by 'keyword' to keep. If passed as empty vector, all
%               matching projects will be taken
% dropboxFolder: Pass this option, followed by the path to data folder 
%                where you wish to save
%                pipeline-generated data sets and figures. If this
%                var is not specified, output will be saved one level
%                above git repo in the folder structure
% first_nc: script defaults to taking only nc14 traces. If you want
%           earlier traces, pass 'first_nc', followed by desired nuclear cycle
%           number
%
% OUTPUT: nucleus_struct: compiled data set contain key nucleus and
% particle attributes

%% ..............Process compileTraces Options..........%

[fnames, ncs_used, DropboxFolder, Prefixes] = processCompileTracesOptions(varargin);
minDP = 5;
nan_buffer = 2;

%% ...............Extract CompiledParticles data from each file.......%

% initialize filename vectors
cp_filenames = cell(1, numel(fnames)); % particles
cn_filenames = cell(1, numel(fnames)); % protein data
ap_filenames = cell(1, numel(fnames)); % ap info
nc_filenames = cell(1, numel(fnames)); % nuclei
fov_filenames = cell(1, numel(fnames)); % fov info
sp_filenames = cell(1, numel(fnames)); % spot info
pa_filenames = cell(1, numel(fnames)); % particle info that contains z info
for f = 1:numel(fnames)
    thisdir = fnames{f};          
    % append file paths
    cn_filenames{f} = [thisdir '/CompiledNuclei.mat'];
    cp_filenames{f} = [thisdir '/CompiledParticles.mat'];    
    ap_filenames{f} = [thisdir '/APDetection.mat'];    
    nc_filenames{f} = [thisdir '/' Prefixes{f} '_lin.mat'];           
    fov_filenames{f} = [thisdir '/FrameInfo.mat'];
    sp_filenames{f} = [thisdir '/Spots.mat'];
    pa_filenames{f} = [thisdir '/Particles.mat'];
end

%% .......Generate Structure for info about every nucleus...........%

nucleus_struct = [];
% Loop through filenames    
for i = 1:length(cp_filenames) 
    
    % read in raw files
    try
        load(nc_filenames{i}) % Ellipse Info
        load(fov_filenames{i}) % FrameInfo Info        
        load(pa_filenames{i}); % Raw Particles
        if ~iscell(Particles)
            Particles = {Particles};
        end
    catch
        continue
    end
    
    [Date, ExperimentType, ExperimentAxis, CoatProtein, StemLoopEnd, APResolution,...
    Channel1, Channel2, Objective, Power, DataFolder, DropboxFolderName, Comments,...
    nc9, nc10, nc11, nc12, nc13, nc14, CF,Channel3,prophase,metaphase, anaphase] = getExperimentDataFromMovieDatabase(Prefixes{i}, DropboxFolder);
    
    % loads AP information
    ap_flag = 1;
    try
        load(ap_filenames{i}) % AP Info 
        % get angle between the x-axis and the AP-axis 
        APAngle = round(atan((coordPZoom(2)-coordAZoom(2))/(coordPZoom(1)-coordAZoom(1)))*360 / (2*pi));    
        %Correction for if APAngle is in quadrants II or III
        if coordPZoom(1)-coordAZoom(1) < 0       
            APAngle = APAngle + 180;
        end  
    catch
        warning('No AP Data detected. Proceeding without AP info')
        ap_flag = 0;
        APAngle = 0;
    end   
    
     
    
    % extract data structures
    processed_data = load(cp_filenames{i}); % processed particles  

    % Extract compiled particles structure         
    cp_particles = processed_data.CompiledParticles;    
    if ~iscell(cp_particles)
        cp_particles = {cp_particles};
    end   
    cp_z_flag = 0;
    if isfield(cp_particles,'zPos')
        cp_z_flag = 1;
    else        
        load(sp_filenames{i}) % Spots info
        if ~iscell(Spots)
            Spots = {Spots};
        end
    end
    
    cp_protein = {}; % empty iff there is no protein data
    if contains(ExperimentType, 'input')
        protein_data = load(cn_filenames{i}); % processed nulcei 
        % extract protein data
        cp_protein = protein_data.CompiledNuclei;   
        if ~iscell(cp_protein)
            cp_protein = {cp_protein};
        end
    end
    
    % sets identifier (ordered by date)
    setID = i;            
    % pull trace and nucleus variables
    time_raw = processed_data.ElapsedTime*60; % time vector   

    traces_raw = processed_data.AllTracesVector; % array with a column for each trace 
    % check to see if traces are stored in cell array
    if ~iscell(traces_raw)
        traces_raw = {traces_raw};
    end
    
    frames_raw = 1:length(time_raw); % Frame list   
    num_outputs = length(traces_raw);
    
    % loops through nuclear cycles
    for nc = ncs_used
        
        % check if we have data taken during this nuclear cycle; skips it
        % if we don't
        try
            if ~(processed_data.(['nc' num2str(nc)]) || processed_data.(['nc' num2str(nc+1)]) > 1)
                continue
            end
        catch
            continue
        end
        
        % finds beginning of nuclear cycle
        first_frame = max([1, processed_data.(['nc' num2str(nc)])]);
        
        % finds end (or last frame we have) of the nuclear cycle
        try
            last_frame = processed_data.(['nc' num2str(nc + 1)]) - 1;
            if isnan(last_frame) || last_frame <= 1
                last_frame = length(time_raw);
            end
        catch
            last_frame = length(time_raw);
        end

        % takes times within nuclear cycle
        time_clean = time_raw(first_frame:last_frame);    
        time_clean = time_clean - min(time_clean); % Normalize to start of nc    
        frames_clean = frames_raw(first_frame:last_frame);    

        [s_cells] = compileSchnitz(schnitzcells, frames_clean, setID, i, ...
            cp_filenames, cp_protein, nc, ap_flag, FrameInfo, num_outputs, time_clean);

        % now add particle info

        % Index vector to cross-ref w/ particles            
        e_index = [s_cells.Nucleus]; 
        nc_index = [s_cells.ncStart];
        
        % iterates through multiple colors
        for cidx = 1:length(traces_raw)
            traces_clean = traces_raw{cidx}(first_frame:last_frame,:);

            % iterate through traces 
            for j = 1:size(traces_clean,2)  
                % Raw fluo trace
                raw_trace = traces_clean(:,j); 
                % Get nucleus ID
                schnitz = cp_particles{cidx}(j).schnitz;
                
                % find corresponding nucleus index        
                nc_ind = find(e_index==schnitz & nc_index == nc);        
                if length(nc_ind) ~= 1
                    warning('Problem with Particle-Nucleus Crossref')
                    continue
                end 
                
                % skip particles not in nc range
                if sum(~isnan(raw_trace)) == 0
                    continue
                end
                
                trace_start = find(~isnan(raw_trace),1);
                trace_stop = find(~isnan(raw_trace),1,'last');        
                %Creat versions with all intervening frames present (missing frames
                %appear as NaNs)
                trace_full = raw_trace(trace_start:trace_stop)';
                frames_full = frames_clean(trace_start:trace_stop); 
                
                % Find intersection btw full frame range and CP frames        
                raw_pt_frames = cp_particles{cidx}(j).Frame;
                cp_frames = raw_pt_frames(ismember(raw_pt_frames,frames_clean));        
                qc_flag = 1;
                if numel(cp_frames) < minDP
                    qc_flag = 0;
                end
                
                % gets other fluorescence measures
                raw_trace3 = nan(1, length(frames_clean));
                raw_trace5 = nan(1, length(frames_clean));
                raw_trace3D = nan(1, length(frames_clean));
                if isfield(cp_particles{cidx}(j), 'Fluo3')
                    raw3 = nan(1, length(frames_raw));
                    raw3(raw_pt_frames) = cp_particles{cidx}(j).Fluo3;
                    raw_trace3 = raw3(first_frame:last_frame);
                end
                if isfield(cp_particles{cidx}(j), 'Fluo5')
                    raw5 = nan(1, length(frames_raw));
                    raw5(raw_pt_frames) = cp_particles{cidx}(j).Fluo5;
                    raw_trace5 = raw5(first_frame:last_frame);
                end
                if isfield(cp_particles{cidx}(j), 'FluoGauss')
                    raw3D = nan(1, length(frames_raw));
                    raw3D(raw_pt_frames) = cp_particles{cidx}(j).FluoGauss;
                    raw_trace3D = raw3D(first_frame:last_frame);
                end
                
                % stores z information if necessary
                if ~cp_z_flag
                    % Initialize arrays to store Z info
                    bZ = NaN(1, length(frames_full));        
                    % Get spot Z location info
                    part_id = cp_particles{cidx}(j).OriginalParticle;        
                    particle_frames_raw = Particles{cidx}(part_id).Frame;
                    for id = 1:length(cp_frames)            
                        abs_idx = cp_frames(id) == frames_full;
                        rel_id = raw_pt_frames==cp_frames(id);
                        if Particles{cidx}(part_id).Index(particle_frames_raw==cp_frames(id)) > length(Spots{cidx}(cp_frames(id)) ...
                                .Fits)
                            error('Mismatch between Particles and Spots dimensions')                                
                        else
                            bZ(abs_idx) = Spots{cidx}(cp_frames(id)) ...
                                    .Fits(Particles{cidx}(part_id).Index(rel_id)).brightestZ;                
                        end
                    end
                end
                                                        
                % find overlap between nucleus and trace
                nc_frames = s_cells(nc_ind).frames;         
                spot_filter = ismember(nc_frames,frames_full);            
                s_cells(nc_ind).(['spot_frames' num2str(cidx)]) = spot_filter;
                if sum(spot_filter) < numel(frames_full)

                    error('Inconsistent particle and nucleus frames')
                end
                % record fluorescence info             
                s_cells(nc_ind).fluo(cidx,:) = raw_trace';
                s_cells(nc_ind).fluo3(cidx,:) = raw_trace3;
                s_cells(nc_ind).fluo5(cidx,:) = raw_trace5;
                s_cells(nc_ind).fluo3D(cidx,:)= raw_trace3D;
                s_cells(nc_ind).APAngle = APAngle;
                % x and y info                                
                s_cells(nc_ind).xPosParticle(cidx,ismember(nc_frames,cp_frames)) = ...
                    cp_particles{cidx}(j).xPos(ismember(cp_frames,nc_frames));
                s_cells(nc_ind).yPosParticle(cidx,ismember(nc_frames,cp_frames)) = ...
                    cp_particles{cidx}(j).yPos(ismember(cp_frames,nc_frames));
                % add z info       
                if ~cp_z_flag
                    s_cells(nc_ind).zPosParticle(cidx,spot_filter) = bZ;   
                else
                    s_cells(nc_ind).zPosParticle(cidx,ismember(nc_frames,cp_frames)) = ...
                    cp_particles{cidx}(j).zPos(ismember(cp_frames,nc_frames));
                end
                % add ap info
                if ap_flag
                    s_cells(nc_ind).apPosParticle(cidx,ismember(nc_frames,cp_frames)) = ...
                    cp_particles{cidx}(j).APpos(ismember(cp_frames,nc_frames));
                end
                % Identifier variables                        
                particle = cp_particles{cidx}(j).OriginalParticle;                        
                s_cells(nc_ind).ParticleID(cidx) = eval([num2str(setID) '.' sprintf('%04d',particle)]);   
                s_cells(nc_ind).qc_flag{cidx} = qc_flag; 
            end      
        end
        nucleus_struct = [nucleus_struct  s_cells]; 
    end
end

%% ..................Interpolate and Clean Traces........................
% estimate interpolation time
med_time = nanmedian(diff([nucleus_struct.time]));
TresInterp = round(med_time);
[nucleus_struct] = interpolateTraces(nucleus_struct, minDP, ...
    TresInterp, num_outputs, nan_buffer); 

%% ..................Save Trace and Nucleus Data .........................
mkdir([DropboxFolder,filesep,project]);
save([DropboxFolder,filesep,project,filesep,'nucleus_struct.mat'],...
    'nucleus_struct','-v7.3');
end

