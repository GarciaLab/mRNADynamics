% DESCRIPTION
% This function generates a GUI to explore different ways of generating a
% nuclear channel. You can invert a channel and you can combine multiple
% together. You can also change the projectiontype that is used. The final
% combination of channels and projection type used when pressing 'Save
% Channel Selection' will be what will be used when forming the
% nuclear/histone channel. Note that this does not change your
% MovieDatabase entry--just the nuclear/histone images generated. At the
% moment this function should only be called downstream of
% exportDataForLivemRNA and is only usable with Leica data (and the
% 'nuclearGUI' option must be entered in exportDataForLivemRNA)

function [projectionChannels, ProjectionType] =...
    chooseIHProjections(Prefix, varargin)

cleanupObj = onCleanup(@myCleanupFun);
warning('off', 'MATLAB:ui:Slider:fixedHeight')


skip_factor = 1; % Only uses 1/skip_factor frames

% default custom projection parameters
max_custom = 1; % highest histone channel slice used
min_custom = 5; % lowest histone channel slice used


ProjectionType = 'midsumprojection';
load('ReferenceHist.mat', 'ReferenceHist');

%options must be specified as name, value pairs. unpredictable errors will
%occur, otherwise.
for arg = 1:2:(numel(varargin)-1)
    if arg ~= numel(varargin)
        eval([varargin{arg} '=varargin{arg+1};']);
    end
end


liveExperiment = LiveExperiment(Prefix);



FullRepsMovieMat = getMarkAndFindMovieMat(liveExperiment);
movieMat = getFirstRepMat(liveExperiment);

chooseIHHisProjections(Prefix,'movieMat', movieMat,'FullRepsMovieMat', FullRepsMovieMat);
chooseIHMemProjections(Prefix,'movieMat', movieMat,'FullRepsMovieMat', FullRepsMovieMat);

end