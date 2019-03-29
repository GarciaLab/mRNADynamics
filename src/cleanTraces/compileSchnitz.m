function [s_cells] = compileSchnitz(schnitzcells, frames_clean, setID, i, ...
    cp_filenames, cp_protein, nc, ap_flag, FrameInfo, num_outputs, time_clean)
%COMPILESCHNITZ Summary of this function goes here
%   Detailed explanation goes here

    % compile schnitz info
    s_cells = struct;
    e_pass = 1;    
    for e = 1:length(schnitzcells)
        e_frames = schnitzcells(e).frames;
        nc_filter = ismember(e_frames,frames_clean);
        nc_frames = e_frames(nc_filter);
        if length(nc_frames) >= 1 % skip nuclei not desired nc range                     
            %Will be set to particle real values for nuclei with matching
            %particle
            for cidx = 1:num_outputs
                s_cells(e_pass).(['ParticleID' num2str(cidx)]) = NaN;
                s_cells(e_pass).(['xPosParticle' num2str(cidx)]) = NaN(1,sum(nc_filter));
                s_cells(e_pass).(['yPosParticle' num2str(cidx)]) = NaN(1,sum(nc_filter));
                s_cells(e_pass).(['zPosParticle' num2str(cidx)]) = NaN(1,sum(nc_filter));            
                s_cells(e_pass).(['fluo' num2str(cidx)]) = NaN(1,sum(nc_filter));
            end
            % add core nucleus info
            x = schnitzcells(e).cenx;            
            y = schnitzcells(e).ceny;                           
            s_cells(e_pass).xPos = x(nc_filter);
            s_cells(e_pass).yPos = y(nc_filter); 
            s_cells(e_pass).frames = nc_frames';            
            s_cells(e_pass).Nucleus = e; 

            s_cells(e_pass).ncID = eval([num2str(setID) '.' sprintf('%04d',e)]);
            s_cells(e_pass).xMean = mean(x(nc_filter));
            s_cells(e_pass).yMean = mean(y(nc_filter));
            s_cells(e_pass).ncStart = nc;
            % time and set info
            s_cells(e_pass).time = time_clean(ismember(frames_clean,nc_frames));
            if numel(s_cells(e_pass).time) > 1
                if s_cells(e_pass).time(2) == 0
                    error('wtf')
                end
            end
            s_cells(e_pass).setID = setID;
            fn = cp_filenames{i}; % Get filename to store in struct  
            fn = fn(1:strfind(fn,'/')-1);
            s_cells(e_pass).source_path = fn;                        
            s_cells(e_pass).PixelSize = FrameInfo(1).PixelSize;       
            s_cells(e_pass).ap_flag = ap_flag;
            % add protein info      
            if ~isempty(cp_protein)
                for j = 1:length(cp_protein)
                    prot_nuc = cp_protein{j}([cp_protein{j}.schnitz] == e);                     
                    pt_vec = prot_nuc.FluoMax(nc_filter);
                    s_cells(e_pass).(['protein' j]) = pt_vec';
                end
            end
            e_pass = e_pass + 1;
        end
    end
end

