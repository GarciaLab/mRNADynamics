function writeScriptArgsToDataStatus(DropboxFolder, DataType, Prefix, args)

    functionString = 'Made filtered spot channel files';
    
    args =  [args{:}];
    D=dir([DropboxFolder,filesep,'DataStatus.*']);
    [~,StatusTxt]=xlsread([DropboxFolder,filesep,D(1).name],DataType);
    PrefixRow= find(strcmpi(StatusTxt(:,1),'Prefix:'));
    functionRow=find(strcmpi(StatusTxt(:,1),functionString));
    PrefixCol=find(contains(StatusTxt(PrefixRow,:),Prefix));
    
    alphabet = ['ABCDEFGHIJKLMNOPQRSTUVWXYZ'];
    longphabet = {};
    for i = 1:26
        longphabet{i} = ['A', alphabet(i)];
    end
    if PrefixCol <= 26
        charCol = alphabet(PrefixCol);
    elseif PrefixCol > 26 && PrefixCol <= 52
        charCol = longphabet{i};
    else
        warning('not supported above 52 prefixes');
    end
    charRange = [charCol, num2str(functionRow)];
    
    str = 'filterMovie(';
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
            xlswrite([DropboxFolder,filesep,D(1).name], {str}, DataType, charRange);
        catch
            error('is data status open? close it then try again');
        end
    else
        %to prevent rewriting multiple cells on accident
        warning('failed to properly write function arguments');
    end
    
    
end

