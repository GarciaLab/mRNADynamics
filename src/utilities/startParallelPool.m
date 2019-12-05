function startParallelPool(nWorkers, displayFigures, keepPool)

ps = parallel.Settings;
ps.Pool.AutoCreate = false;
distcomp.feature( 'LocalUseMpiexec', false )

if nWorkers > 1 && ~displayFigures
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

end