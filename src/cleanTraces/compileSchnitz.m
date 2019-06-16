function [s_cells] = compileSchnitz(schnitzcells, frames_clean, setID, i, ...
    cp_filenames, cp_protein, nc, ap_flag, FrameInfo, num_outputs, time_clean, ap_ref_grid)
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
            s_cells(e_pass).ParticleID = NaN(num_outputs,1);            
            s_cells(e_pass).xPosParticle = NaN(num_outputs,sum(nc_filter));
            s_cells(e_pass).yPosParticle= NaN(num_outputs,sum(nc_filter));
            s_cells(e_pass).zPosParticle = NaN(num_outputs,sum(nc_filter));             
            s_cells(e_pass).fluo = NaN(num_outputs,sum(nc_filter));
            s_cells(e_pass).fluo3D = NaN(num_outputs,sum(nc_filter));

            % add core nucleus info
            x = schnitzcells(e).cenx;            
            y = schnitzcells(e).ceny;  
            s_cells(e_pass).xPos = round(x(nc_filter));
            s_cells(e_pass).yPos = round(y(nc_filter)); 
            apPos = NaN(size(x(nc_filter)));
            for m = 1:numel(apPos)
                apPos(m) = ap_ref_grid(s_cells(e_pass).yPos(m),s_cells(e_pass).xPos(m));
            end
            s_cells(e_pass).apPos = apPos;            
            s_cells(e_pass).frames = nc_frames';            
            s_cells(e_pass).Nucleus = e; 

            s_cells(e_pass).ncID = eval([num2str(setID) '.' sprintf('%04d',e)]);
            s_cells(e_pass).xMean = mean(x(nc_filter));
            s_cells(e_pass).yMean = mean(y(nc_filter));
            s_cells(e_pass).ncStart = nc;
            % time and set info 
            s_cells(e_pass).time = time_clean(ismember(frames_clean,e_frames));

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
                    s_cells(e_pass).protein(j,:) = pt_vec';
                end
            end
            e_pass = e_pass + 1;
        end
    end
end

