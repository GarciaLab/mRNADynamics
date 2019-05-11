function startParallelPool(nWorkers, varargin)

ps = parallel.Settings;
ps.Pool.AutoCreate = false;

if ~isempty(varargin)
    displayFigures = varargin{1};
else
    displayFigures = false;
end

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
    try %#ok<TRYNC>
        poolobj = gcp('nocreate');
        delete(poolobj);
    end
end

end