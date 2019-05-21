function argin = tryGPU(~)
    argin = {};

    try
        gpuDevice(1);
        argin = {'gpuArray'};
    catch
        warning('nvidia cuda gpu unavailable.. defaulting to cpu');
    end

end