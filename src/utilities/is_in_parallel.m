function answer = is_in_parallel ()

%from https://www.mathworks.com/matlabcentral...
%/answers/58228-am-i-running-in-parallel-best-way-to-check

try
    answer = ~isempty(getCurrentTask());
catch err
    if ~strcmp(err.identifier, 'MATLAB:UndefinedFunction')
        rethrow(err);
    end
    answer = false;
end


end
