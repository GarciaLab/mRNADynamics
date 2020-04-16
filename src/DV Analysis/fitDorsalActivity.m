function [fit, model] = fitDorsalActivity(dlfluobins, activity, varargin)
%possible model types- {'hill', 'simpleWithPol', 'mcwNoPol'}
arguments
    dlfluobins double
    activity double
end
arguments(Repeating)
    varargin
end

modelType = 'hill';
% xScale =10^-(round(log10(max(x(:)))));
xScale = 1;
% yScale = 10^-(round(log10(max(y(:)))));
yScale = 1;
y = activity(:);
x = dlfluobins(:);
x(isnan(y)) = [];
y(isnan(y)) = [];

%scale the problem for better fitting
x = x.*xScale;
y = y.*yScale;
%p(1)=rate coefficient, p(2)=kd, p(3)=hill coefficient p(4) y offset


%rate, kd, hill, y offset 
p0 = [max(y),max(x)/2, 1, 0];
lb = [0, 1000.*xScale, 2, -max(y)];
ub = [max(y)*2, Inf, 6, max(y)*10];

for i = 1:length(varargin)
   if strcmpi(varargin{i}, 'fix2')
        kd = varargin{i+1};
        p0(2) = kd; lb(2) = kd; ub(2) = kd;
     elseif strcmpi(varargin{i}, 'fix1')
        r = varargin{i+1};
        p0(1) = r; lb(1) = r; ub(1) = r;
     elseif strcmpi(varargin{i}, 'fix4')
        off = varargin{i+1};
        p0(4) = off; lb(4) = off; ub(4) = off;
    elseif strcmpi(varargin{i}, 'fix3')
        n = varargin{i+1};
        p0(3) = n; lb(3) = n; ub(3) = n;
     elseif strcmpi(varargin{i}, 'fix5')
        w = varargin{i+1};
        p0(5) = w; lb(5) = w; ub(5) = w;
   elseif strcmpi(varargin{i}, 'modelType')
       modelType = varargin{i+1};
    end
end

%common model string pieces
offTerm = '+ p(4)';
ampCoeff = 'p(1)*';

%%
if strcmpi(modelType, 'hill')
    % amp kd hill y-offset
    numeratorTerm = '(d/p(2)).^p(3)';
    partitionTerm = '1 + ((d)/p(2)).^p(3)'; 
end
%%
if strcmpi(modelType, 'mwcNoPol')     % amp kd delta opening energy y-offset
    numeratorTerm =  'exp(p(3))*(d/p(2))';
    partitionTerm = '1+exp(p(3)) + exp(p(3))*(d/p(2))';
    
    off = 0;
    p0(4) = off; lb(4) = off; ub(4) = off;

end
%%
if strcmpi(modelType, 'simpleWithPol')     %simple binding with polymerase
    %p(1)=rate p(2)=Kdd p(3)=[P]/Kdp p(4)=y offset p(5)= omegaDP
    numeratorTerm = 'p(3) + (d/p(2))*p(3)*p(5)';
    partitionTerm = '1 + d/p(2) + p(3) + (d/p(2))*p(3)*p(5)';
    
    %add a fifth parameter
    p0 = [p0, 1]; lb = [lb, 0]; ub = [ub, Inf];
    %no y offset for this model
    off = 0;
    p0(4) = off; lb(4) = off; ub(4) = off;
end
%%
modelStr =  ['@(p,d) ', ampCoeff, '(',...
    '(',numeratorTerm,')',...
    './',...
    '(',partitionTerm,')',...
    ')', offTerm];

model = str2func(modelStr);


options = optimoptions(@lsqcurvefit, 'MaxFunctionEvaluations', length(p0)*1000, 'MaxIterations', 4000,...
    'OptimalityTolerance',1E-10,'FunctionTolerance',1E-10, 'Display','none');

fit = lsqcurvefit(model, p0, x, y, lb, ub, options);

%rescale 
fit(2) = fit(2)./xScale; %kd (x when y is half the plateau)
fit(1) = fit(1)./yScale; %rate coefficient (plateau)
fit(4) = fit(4)./yScale; % y offset

% if length(fit)>4
%     fit(5) = fit(5)./xScale; %kd (x when y is half the plateau)
% end



end