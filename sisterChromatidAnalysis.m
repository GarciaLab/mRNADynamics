% load 'D:\Data\Armando\livemRNA\Data\Dropbox\DynamicsResults\2015-05-31-89B8-3-P2P\Spots.mat'
% load 'D:\Data\Armando\livemRNA\Data\Dropbox\DynamicsResults\2015-05-31-89B8-3-P2P\FrameInfo.mat'
% uiopen
% uiopen
a = [];
t = [];
l = length(FrameInfo);
for i = 1:length(Spots)
    for j = 1:length(Spots(i).Fits)
        for k = 1:length(Spots(i).Fits(j).z)
            d = Spots(i).Fits(j).sisterSeparation{k};
            t = [t,Spots(i).Fits(j).frame];
            if ~isempty(d)
                a = [a, d];
                imshow(imresize(Spots(i).Fits(j).rawSpot{k},10),[])
            else 
                a = [a, 0];
            end
        end
    end
end

m = zeros(1,length(FrameInfo));
l = 1:length(FrameInfo);
n = zeros(1,length(FrameInfo));
for i = 1:length(l)
    for j = 1:length(t)
        if t(j) == i
            m(i) = m(i) + a(j);
            n(i) = n(i)+1;
        end
    end
end


plot(l.*20.82/60, m./n)
