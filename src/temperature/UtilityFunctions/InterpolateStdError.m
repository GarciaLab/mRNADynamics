function vqstd = InterpolateStdError(x, vstd, xq)

if (sum(size(vstd) == 1) > 0)
    vqstd = NaN(1, length(xq));
    for i = 1:length(xq)
        lowbound = find(x <= xq(i), 1, 'last');
        highbound = find(x >= xq(i), 1);
        if isempty(lowbound) | isempty(highbound)
            continue
        end
        if lowbound == highbound
            vqstd(i) = vstd(lowbound);
        else
            df1 = 1-(xq(i)-x(lowbound))/(x(highbound)-x(lowbound));
            df2 = (xq(i)-x(lowbound))/(x(highbound)-x(lowbound));
            vqstd(i) = sqrt(df1^2*vstd(lowbound)^2+df2^2*vstd(highbound)^2);
        end
    end
else
    vqstd = NaN(length(xq), size(vstd, 2));
    
    for i = 1:length(xq)
        lowbound = find(x <= xq(i), 1, 'last');
        highbound = find(x >= xq(i), 1);
        if isempty(lowbound) | isempty(highbound)
            continue
        end
        for col = 1:size(vstd, 2)
            if lowbound == highbound
                vqstd(i, col) = vstd(lowbound, col);
            else
                df1 = 1-(xq(i)-x(lowbound))/(x(highbound)-x(lowbound));
                df2 = (xq(i)-x(lowbound))/(x(highbound)-x(lowbound));
                vqstd(i, col) = sqrt(df1^2*vstd(lowbound, col)^2+df2^2*vstd(highbound, col)^2);
            end
        end
        
    end
    
end