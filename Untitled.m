PreProcPath = 'P:\Simon\LivemRNA\Data\PreProcessedData';

Prefixes = {'2021-02-21-UBC1_cer_60G_3xpower','2021-02-21-UBC1_cer_60G_4xpower',...
    '2021-02-21-UBC1-cer-60G_5xpower','2021-02-21-UBC1-cer-60G_6xpower'};

%%
figure
hold on
ch = 1;
sigma = 5;
Palette = viridis(length(Prefixes));
for p = 1:length(Prefixes)
    Prefix = Prefixes{p};
    PrefixFluos = [];
    PrefixPixelValues = [];
    for frame = 1%:10
        Stack = getMovieFrame(LiveExperiment(Prefix), frame, 1);
        StackDOG = imgaussfilt3(Stack,sigma);
        MaxProjectStackDOG = max(StackDOG,[],3);
        PrefixPixelValues = [PrefixPixelValues MaxProjectStackDOG(:)];
    end
    H = histogram(PrefixPixelValues,'DisplayStyle','stairs','EdgeColor',Palette(p,:),'LineWidth',2,...
        'Normalization','pdf');
    HistVals = H.values;
    [value idx] = find(HistVals==max(HistVals));
    
end
hold off
legend('3','4','5','6')
%     %segmentSpots(Prefix,4500,'nWorkers',8)
%     PrefixPreProcPath = [PreProcPath '/' Prefix];
%     GFPFiles = dir([PrefixPreProcPath '/*ch01.tif']);
%     
%     for s = 1:10
%         seriesFileName = GFPFiles(s).name;
%         
%     end
%     
% end


