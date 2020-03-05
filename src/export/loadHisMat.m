function hisMat = loadHisMat(hisFile, frameRange)

    hismatfile = matfile(hisFile);
    if ~isempty(frameRange)
        hisMat = hismatfile.hisMat(:, :, frameRange);
    else
        hisMat = hismatfile.hisMat;
    end

end