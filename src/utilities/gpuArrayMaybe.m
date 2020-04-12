function array = gpuArrayMaybe(cpuArray)
 
%wrapper for gpuArray that checks if
%gpu is present and returns cpu array if not
    try
        gpuDevice(1);
        array = gpuArray(cpuArray);
    catch
        warning('NVIDIA CUDA GPU unavailable. Defaulting to CPU.');
        array = cpuArray;
    end

end