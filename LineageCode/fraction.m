%Fraction of Active Nuclei
numect
for i=nc12:numofframes
    %if (inverted=='u')
        for ii=1:size(Ellipses{i,1},1)
            if Ellipses{i,1}(ii,2)>384 || Ellipses{i,1}(ii,2)<128
                numecto=numecto+1;
            end
            if Ellipses{i,1}(ii,2)<384 && Ellipses{i,1}(ii,2)>128
                nummeso=nummeso+1;
            end
        end
end

        
        
        
    