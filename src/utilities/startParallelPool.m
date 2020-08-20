function startParallelPool(nWorkers, displayFigures, keepPool)

licensed = license('test','Distrib_Computing_Toolbox');

if licensed 
    
    try
        
    ps = parallel.Settings;
    ps.Pool.AutoCreate = false;
    distcomp.feature( 'LocalUseMpiexec', false );

    if nWorkers > 1 & ~displayFigures
        maxWorkers = nWorkers;
        ps.Pool.AutoCreate = true;

        try
            parpool(maxWorkers);
        catch
            try
                parpool; % in case there aren't enough cores on the computer
            catch
                % parpool throws an error if there's a pool already running.
            end
        end  
    else
        if ~keepPool
            try %#ok<TRYNC>
                poolobj = gcp('nocreate');
                delete(poolobj);
            end
        end
    end
    catch
        %sometimes this fails if the license is available but the toolbox
        %isn't installed. 
    end
    
end

end