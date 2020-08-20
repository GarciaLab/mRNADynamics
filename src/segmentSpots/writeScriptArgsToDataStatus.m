function writeScriptArgsToDataStatus(DropboxFolder, dataType, Prefix, args, functionString1, functionString2)
    
    %sometimes the structure of the args is different. this hopefully
    %accounts for that so there's a cell array in the end
    args0 = args;
    args1 = [args0{:}];
    if iscell(args1)
        args = args1;
    else
        args = args0;
    end
    
    alphabet = ['ABCDEFGHIJKLMNOPQRSTUVWXYZ'];
    longphabet = {};
    for i = 1:26
        longphabet{i} = ['A', alphabet(i)];
    end

        
    D=dir([DropboxFolder,filesep,'DataStatus.*']);
    [~,StatusTxt]=xlsread([DropboxFolder,filesep,D(1).name],dataType);
    PrefixRow= find(strcmpi(StatusTxt(:,1),'Prefix:'));
    
    
    functionRow=find(contains(StatusTxt(:,1),functionString1, 'IgnoreCase',true));
    PrefixCol=find(contains(StatusTxt(PrefixRow,:),Prefix));
    
    newPrefix = false;
    if isempty(PrefixCol)
        newPrefix = true;
        %we can create the moviedatabase entry.
        PrefixCol = size(StatusTxt, 2); %the column after the latest entry
    end
    
    if PrefixCol <= 26
        charCol = alphabet(PrefixCol);
    elseif PrefixCol > 26 && PrefixCol <= 52
        charCol = longphabet{i};
    else
        warning('not supported above 52 prefixes');
    end
    charRange = [charCol, num2str(functionRow)];
    
    str = [functionString2, '('];
    for i = 1:length(args)
        if i == 1
            str = [str, args{i}];
        else
            str = [str, ',', num2str(args{i})];
        end
    end
    str = [str, ')'];
    if ~iscell(str)
        try
            xlswrite([DropboxFolder,filesep,D(1).name], {str}, dataType, charRange);
            if newPrefix
                prefixRange = [charCol, num2str(PrefixRow)];
                prefixString = ['''prefix=''', Prefix, ''''];
                xlswrite([DropboxFolder,filesep,D(1).name], {prefixString}, dataType, prefixRange);
            end
        catch
            error('is data status open? close it then try again');
        end
    else
        %to prevent rewriting multiple cells on accident
        warning('failed to properly write function arguments');
    end
    
    
end

