function array = gpuArrayWrapper(cpuArray)
 
%wrapper for gpuArray that checks if
%gpu is present and returns cpu array if not
    try
        array = gpuArray(cpuArray);
    catch
        warning('NVIDIA CUDA GPU unavailable. Defaulting to CPU.');
        array = cpuArray;
    end

end