function save2(file, var)
%
%save input as either v6 (preferred)
%or v7.3 if necessary. Yes, the eval is 
%necessary. It's so that the variable
%is saved as its original name and not
%as 'var'. 
%AR 4/20


varStr = inputname(2); 
eval([varStr, '=var;']);

if whos(varStr).bytes < 2E9
    save(file, varStr,  '-v6');
else
    save(file,varStr, '-v7.3');
end


end