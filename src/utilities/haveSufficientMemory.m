function haveSufficientMemory = haveSufficientMemory(directory)

try
    [~,sys] = memory;
catch
    %if this fails, assume there's not enough memory
    haveSufficientMemory = false;
    return;
end
movieSize = 0;
for k = 1:length(directory)
    
    movieSize = movieSize + directory(k).bytes;
    
end

%movie shouldn't use 90% of available memory
if .9*sys.PhysicalMemory.Available < movieSize
    haveSufficientMemory = false;
    warning("Don''t have sufficient memory to load entire"+...
        "movie.");
else
    haveSufficientMemory = true;
end

end