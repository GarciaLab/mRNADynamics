function D = munge(T, k, p, s)

%generate synthetic data
%
%T- training data matrix n x d
%k- size parameter
%p- probability parameter
%s- local variance parameter

numAttributes = size(T, 2);
numObservations = size(T, 1);

D = zeros(k*numObservations, numAttributes);


for i = 1:k
    T0 = T;
    neighbors = knnsearch(T0, T0, 'K', 2, 'Distance', 'euclidean');
    for j = 1:numObservations
        e0 = T0(neighbors(j, 2), :);
        for a = 1 : numAttributes-1 %don't swap labels
            
            r = rand;
           
            if r < p
                
                %swap continous
                e = T0(j, :);
                sd = abs(e(a)-e0(a))/s;   
                      
                T0(j, a) = normrnd(e0(a), sd);
                e0(a) = normrnd(e(a), sd);
                
                %swap nominal
%                 val = T0(j, a);
%                 T0(j, a) = e0(a);
%                 e0(a) = val;
                
                
                
            end
        end
    end
    
    ind0 = ( (i-1)*numObservations) + 1;
    ind1 = i*numObservations;
    D(ind0:ind1, :) = T0;

end


end