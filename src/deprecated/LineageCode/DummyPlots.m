%% Compute Statistics
% Define constant ectogenic boundary
%function StatPlots(lrep,alltracesy,diff,frames,nc12)

%display(['The area being checked is ',num2str(diff),' pixels large']);
means1=zeros(frames,1);
means2=zeros(frames,1);

min1=zeros(frames,1);
min2=zeros(frames,1);

max1=zeros(frames,1);
max2=zeros(frames,1);

std1=zeros(frames,1);
std2=zeros(frames,1);

count1=zeros(frames,1);
count2=zeros(frames,1);

sum1=zeros(frames,1);
sum2=zeros(frames,1);

er1=zeros(frames,1);
er2=zeros(frames,1);

for i=nc12:frames
    %if (inverted=='u')
        traces1=find(alltracesy{i,1}(:)>128 & alltracesy{i,1}(:)<384); % 1 - Mesoderm, 2- Neuro
        traces2=find(alltracesy{i,1}(:)<129 | alltracesy{i,1}(:)>383);
    %elseif (inverted=='d')
%         traces1=find(alltracesy{i,1}(:)<lrep(2) & alltracesy{i,1}(:)>lrep(2)-diff);
%         traces2=find(alltracesy{i,1}(:)>lrep(2) & alltracesy{i,1}(:)<lrep(2)+diff);
    %end
    if ~isempty(traces1)
        min1(i)=min(alltraces{i,1}(traces1));
        if min1(i)<0
            fprintf('Check frame %i with intensity %d\n',i,min1(i));
        end
        max1(i)=max(alltraces{i,1}(traces1));
        std1(i)=std(alltraces{i,1}(traces1));
        count1(i)=length(traces1);
        sum1(i)=sum(alltraces{i,1}(traces1));
        er1(i)=sum(er{i,1}(traces1));
    end
    
    if ~isempty(traces2)
        min2(i)=min(alltraces{i,1}(traces2));
        if min1(i)<0
            fprintf('Check frame %i with intensity %d\n',i,min1(i));
        end
        max2(i)=max(alltraces{i,1}(traces2));
        std2(i)=std(alltraces{i,1}(traces2));
        count2(i)=length(alltraces{i,1}(traces2));
        sum2(i)=sum(alltraces{i,1}(traces2));
        er2(i)=sum(er{i,1}(traces2));
    end
    
    means1(i)=mean(alltraces{i,1}(traces1));
    means2(i)=mean(alltraces{i,1}(traces2));
end

numecto=zeros(frames,1);
nummeso=zeros(frames,1);
for i=nc12:frames
    %if (inverted=='u')
        for ii=1:size(Ellipses{i,1},1)
            if Ellipses{i,1}(ii,2)>383 || Ellipses{i,1}(ii,2)<129
                numecto(i)=numecto(i)+1;
            end
            if Ellipses{i,1}(ii,2)<384 && Ellipses{i,1}(ii,2)>128
                nummeso(i)=nummeso(i)+1;
            end
        end
%     elseif (inverted=='d')
%       for ii=1:size(Ellipses{i,1},1)
%             if Ellipses{i,1}(ii,2)>lrep && Ellipses{i,1}(ii,2)<lrep+diff
%                 nummeso=nummeso+1;
%             end
%             if Ellipses{i,1}(ii,2)<lrep && Ellipses{i,1}(ii,2)>lrep-diff
%                 numecto=numecto+1;
%             end
%       end
%    end
end
    
            

%% Clean up the statistics
for i=nc12:frames
    bd=isnan(max1); % Finds the nans (bad data)
    gd=find(~bd); % Indices of actual points (good data)
    bd([1:(min(gd)-1) (max(gd)+1):end])=0; % Fills extremes with 0
    max1(bd)=interp1(gd,max1(gd),find(bd)); % Interpolates data
    min1(bd)=interp1(gd,min1(gd),find(bd)); % Interpolates data
    std1(bd)=interp1(gd,std1(gd),find(bd)); % Interpolates data
    means1(bd)=0; % Interpolates data
    
    bd=isnan(max2); % Finds the nans (bad data)
    gd=find(~bd); % Indices of actual points (good data)
    bd([1:(min(gd)-1) (max(gd)+1):end])=0; % Fills extremes with 0
    max2(bd)=interp1(gd,max2(gd),find(bd)); % Interpolates data
    min2(bd)=interp1(gd,min2(gd),find(bd)); % Interpolates data
    std2(bd)=interp1(gd,std2(gd),find(bd)); % Interpolates data
    means2(bd)=0; % Interpolates data
end

%% Plot the statistics
figure('Name','Mean Intensity')
plot(nc12:frames, sum1(nc12:frames)./count1(nc12:frames),'g')
hold on
plot(nc12:frames, sum2(nc12:frames)./count2(nc12:frames),'b')
xlabel('Frame (Only nc12, nc13 and nc14 are included)')
ylabel('Mean Intensity (arb. units)')
legend('Mesoderm','Neurogenic Ectoderm');
title('Mean Fluorescence for all Expressing Particles');
saveas(gcf,'Mean Intensity Plots.png');

figure('Name','Minimum Intensity')
plot(nc12:frames, min1(nc12:frames),'g')
hold on
plot(nc12:frames, min2(nc12:frames),'b')
xlabel('Frame (Only nc12, nc13 and nc14 are included)');
ylabel('Minimum Intensity (arb. units)');
legend('Mesoderm','Neurogenic Ectoderm');
title('Minimum Fluorescence for all Expressing Particles');
saveas(gcf, 'Minimum Intensity Plots.png');

figure('Name','Maximum Intensity Plots')
plot(nc12:frames, max1(nc12:frames),'g')
hold on
plot(nc12:frames, max2(nc12:frames),'b')
xlabel('Frame (Only nc12, nc13 and nc14 are included)')
ylabel('Maximum Intensity (arb. units)');
legend('Mesoderm','Neurogenic Ectoderm');
title('Maximum Fluorescence for all Expressing Particles');
saveas(gcf, 'Maximum Intensity Plots.png');

figure('Name','Standard Deviation of Intensity')
plot(nc12:frames, std1(nc12:frames),'g')
hold on
plot(nc12:frames, std2(nc12:frames),'b')
xlabel('Frame (Only nc12, nc13 and nc14 are included)')
ylabel('Standard Deviation of Intensity (arb. units)');
legend('Mesoderm','Neurogenic Ectoderm');
title('Standard Deviation of  Fluorescence for all Expressing Particles');
saveas(gcf, 'Standard Deviation Plots.png');

figure('Name','Density of Particles')
plot(nc12:frames, count1(nc12:frames)./nummeso(nc12:frames),'g')
hold on
plot(nc12:frames, count2(nc12:frames)./numecto(nc12:frames),'b')
xlabel('Frame (Only nc12, nc13 and nc14 are included)')
ylabel('Fraction of Active Nuclei');
legend('Mesoderm','Neurogenic Ectoderm');
title('Fraction of Active Nuclei for Sog');
saveas(gcf, 'Spatial Density Plots.png');

figure('Name','Sum intensity of Particles')
errorbar(nc12:frames, sum1(nc12:frames),zeros(-nc12+frames+1,1),er1(nc12:frames),'g')
hold on
errorbar(nc12:frames, sum2(nc12:frames),zeros(-nc12+frames+1,1),er2(nc12:frames),'b')
xlabel('Frame (Only nc12, nc13 and nc14 are included)')
ylabel('Total Intensity (arb. units)');
legend('Mesoderm','Neurogenic Ectoderm');
title('Total Intensity of all Expressing Particles for Sog');
saveas(gcf, 'Sum Intensity Plots.png');