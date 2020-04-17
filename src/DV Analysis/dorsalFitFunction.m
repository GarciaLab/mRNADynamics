function model = dorsalFitFunction(modelType)
%% generalized way of constructing anonymous
% function from parts specific to models
arguments
    modelType char
end

%common model string pieces
offTerm = '+ p(4)';
ampCoeff = 'p(1)*';

%%
switch modelType
    
    case 'hill'
        % amp kd hill y-offset
        numeratorTerm = '(d/p(2)).^p(3)';
        partitionTerm = '1 + ((d)/p(2)).^p(3)';
        
    case  'mwcNoPol'    % amp kd delta opening energy y-offset
        numeratorTerm =  'exp(p(3))*(d/p(2))';
        partitionTerm = '1+exp(p(3)) + exp(p(3))*(d/p(2))';
        
    case 'simpleWithPol'    %simple binding with polymerase
        %p(1)=rate p(2)=Kdd p(3)=[P]/Kdp p(4)=y offset p(5)= omegaDP
        numeratorTerm = 'p(3) + (d/p(2))*p(3)*p(5)';
        partitionTerm = '1 + d/p(2) + p(3) + (d/p(2))*p(3)*p(5)';
        
    otherwise
        error('no valid model type')
        
end



modelStr =  ['@(p,d) ', ampCoeff, '(',...
    '(',numeratorTerm,')',...
    './',...
    '(',partitionTerm,')',...
    ')', offTerm];



model = str2func(modelStr);

end