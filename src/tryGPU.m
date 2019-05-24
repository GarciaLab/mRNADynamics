function argin = tryGPU(~)
    argin = {};

    try
        gpuDevice(1);
        argin = {'gpuArray'};
    catch
        warning('NVIDIA CUDA GPU unavailable. Defaulting to CPU.');
    end

end