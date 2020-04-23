function vmStressTester()

%cleanup
poolobj = gcp('nocreate');
delete(poolobj);

saveFile = 'X:\vmStressTest.txt';

parTime_s = [];
loadTime_s = [];
loadFolder = "X:\DorsalSynthetics\Data\PreProcessedData\2019-12-28-2Dgc_EfEfEf_4\2019-12-28-2Dgc_EfEfEf_4_movieMat.Mat";

imreadTime_s = [];
imreadFolder = 'X:\Simon\LivemRNA\Data\PreProcessedData\2018-07-24-AL07R-APX2_7_4';

benchTime_s = [];
currentTime_s = strings;
hpcTime_s = [];

nWorkers = 6;
nLoops = 1;

for i = 1:nLoops
    
    %step 0
    currentTime_s(i, 1) = string(datestr(now));
    
    %step 1
    tic;
    pool = parpool(nWorkers, 'SpmdEnabled',true);
    parTime_s(i) = toc;
    
    t = benchHPC(pool);
    if i == 1
        hpcResults = t;
    else
        hpcResults = vertcat(hpcResults, t);
    end
    
    hpcTimes(i) = sum(t.Time);
    
    %step 2
    tic
    load(loadFolder);
    loadTime_s(i) = toc;
    
    %step 3
    tic
    D = dir([imreadFolder, '*.tif']);
    for k = 3:length(D)
        imread([imreadFolder, D(k).name]);
    end
    imreadTime_s(i) = toc;
    
    
    %step 4
    benchTime_s(i) = sum(sum(bench(10)));
    
    %cleanup
    poolobj = gcp('nocreate');
    delete(poolobj);
    
end

T = table(currentTime_s, parTime_s, loadTime_s,...
    imreadTime_s, benchTime_s, hpcTimes);
cell2csv(saveFile, table2cell(T));
writetable(T, saveFile, 'WriteRowNames',...
true, 'WriteVariableNames', true);

end

function allResults = benchHPC(pool)

gbPerWorker = 1;
dataSizes = hpccDataSizes(pool.NumWorkers,gbPerWorker);
hplResult = hpccHPL(dataSizes.HPL);
dgemmResult = hpccDGEMM(dataSizes.DGEMM);
streamResult = hpccSTREAM(dataSizes.STREAM);
ptransResult = hpccPTRANS(dataSizes.PTRANS);
addAttachedFiles(pool,{'hpccRandomNumberGeneratorKernel.mexa64',...
    'hpccRandomNumberGeneratorKernel.mexw64',...
    'hpccRandomNumberGeneratorKernel.mexmaci64'});
randomAccessResult = hpccRandomAccess(dataSizes.RandomAccess);
fftResult = hpccFFT(dataSizes.FFT);
allResults = [hplResult; dgemmResult; streamResult; ...
    ptransResult; randomAccessResult; fftResult];

end