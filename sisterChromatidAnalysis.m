load 'D:\Data\Armando\livemRNA\Data\Dropbox\DynamicsResults\2015-05-31-89B8-3-P2P\Spots.mat'
load 'D:\Data\Armando\livemRNA\Data\Dropbox\DynamicsResults\2015-05-31-89B8-3-P2P\FrameInfo.mat'
% uiopen
% uiopen
distances = [];
times = [];
for i = 1:length(Spots)
    for j = 1:length(Spots(i).Fits)
        for k = 1:length(Spots(i).Fits(j).z)
            d = Spots(i).Fits(j).sisterSeparation{k};
            times = [times,Spots(i).Fits(j).frame];
            if ~isempty(d)
                distances = [distances, d];
%                 imshow(imresize(Spots(i).Fits(j).rawSpot{k},10),[])
            else 
                distances = [distances, 0];
            end
        end
    end
end

dur = length(FrameInfo);
sister_distance_sum = zeros(1,dur);
frames = 1:dur;
n_sisterpairs = zeros(1,dur);
n_sisterpairs_nozeros = zeros(1,dur);
for i = 1:dur
    for j = 1:length(times)
        if times(j) == i
            sister_distance_sum(i) = sister_distance_sum(i) + distances(j);
            n_sisterpairs(i) = n_sisterpairs(i)+1;
            if distances(j)~=0
                n_sisterpairs_nozeros(i) = n_sisterpairs_nozeros(i) + 1;
            end
        end
    end
end

average_sister_separation = sister_distance_sum./n_sisterpairs;
average_nozeros = sister_distance_sum./n_sisterpairs_nozeros;
framescal = frames*(FrameInfo(2).Time)/60;
figure(1)
plot(framescal, average_sister_separation)
title('Sister chromatid visibility in P2P transcription loci')
xlabel('Time (min)')
ylabel('Mean sister separation (um)')

figure(2)
plot(framescal, average_nozeros)
title('Sister chromatid visibility in P2P transcription loci')
xlabel('Time (min)')
ylabel('Mean sister separation (um)')

figure(3)
plot(framescal, n_sisterpairs_nozeros./n_sisterpairs)
title('Sister chromatid visibility in P2P transcription loci')
xlabel('Time (min)')
ylabel('Mean sister separation (um)')


