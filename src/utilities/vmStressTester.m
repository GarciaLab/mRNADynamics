function vmStressTester()

parTimes = [];
loadTimes = [];
imreadTimes = [];
benchTimes = [];
currentTimes = [];

nLoops = 10;

for i = 1:nLoops
    
    %step 0
    currentTimes(i) = string(datestr(now));
    
    %step 1
    tic;
    parpool(20);
    parTimes(i) = toc;
     
    %step 2
    tic
    load("X:\DorsalSynthetics\Data\PreProcessedData\2019-12-28-2Dgc_EfEfEf_4\2019-12-28-2Dgc_EfEfEf_4_movieMat.Mat");
    loadTimes(i) = toc; 
    
    %step 3
    tic
    folder = 'X:\Simon\LivemRNA\Data\PreProcessedData\2018-07-24-AL07R-APX2_7_4';
    D = dir([folder, '*.tif']);
    for k = 3:length(D)
        imread([folder, D(k).name]);
    end
    imreadTimes(i) = toc;
    
    %step 4
    benchTimes(i) = sum(bench(10));
    
    %cleanup
    poolobj = gcp('nocreate');
    delete(poolobj);
end

T = table(currentTimes, parTimes, loadTimes, imreadTimes, benchTimes);

writeTable(T, 'vmStressTester.txt', 'WriteRowNames',...
    true, 'WriteVariablesNames', true, 'WriteMode', 'append');

end