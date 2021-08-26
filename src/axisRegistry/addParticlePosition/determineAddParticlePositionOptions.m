function [SkipAlignment, ManualAlignment, NoAP, SelectChannel, InvertHis,...
    optionalResults, yToManualAlignmentPrompt, correctDV] = determineAddParticlePositionOptions(varargin);

    varargin = varargin{1};

    SkipAlignment = false;
    ManualAlignment = true;
    NoAP = false;
    SelectChannel = 0;
    InvertHis = false;
    optionalResults = '';
    yToManualAlignmentPrompt = false;
    correctDV = false;

    for i = 1:length(varargin)
        switch varargin{i}
            case {'SkipAlignment'}
                disp('Skipping alignment step')
                SkipAlignment = 1;
            case {'ManualAlignment'}
                ManualAlignment = true;
            case {'NoAP'}
                NoAP = 1;
            case {'SelectChannel'}
                SelectChannel = 1;
            case {'optionalResults'}
                optionalResults = varargin{i + 1};
            case {'yToManualAlignmentPrompt'}
                yToManualAlignmentPrompt = 1;
            case {'correctDV'}
                correctDV = true;
        end
    end
end