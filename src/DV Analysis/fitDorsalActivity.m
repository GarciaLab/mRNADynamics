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
% yScale = 10^-(round(log10(max(y(:)))));
xScale = 1;
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

for k = 1:length(varargin)
    if strcmpi(varargin{k}, 'fix2')
        kd = varargin{k+1};
        p0(2) = kd; lb(2) = kd; ub(2) = kd;
    elseif strcmpi(varargin{k}, 'fix1')
        r = varargin{k+1};
        p0(1) = r; lb(1) = r; ub(1) = r;
    elseif strcmpi(varargin{k}, 'fix4')
        off = varargin{k+1};
        p0(4) = off; lb(4) = off; ub(4) = off;
    elseif strcmpi(varargin{k}, 'fix3')
        n = varargin{k+1};
        p0(3) = n; lb(3) = n; ub(3) = n;
    elseif strcmpi(varargin{k}, 'fix5')
        w = varargin{k+1};
        p0(5) = w; lb(5) = w; ub(5) = w;
    elseif strcmpi(varargin{k}, 'modelType')
        modelType = varargin{k+1};
    end
end

switch modelType
    
    case 'hill'
        %nothing extra to do
    case 'simpleWithPol'
        %add a fifth parameter
        p0 = [p0, 1]; lb = [lb, 0]; ub = [ub, Inf];
        %no y offset for this model
        off = 0;
        p0(4) = off; lb(4) = off; ub(4) = off;
    case 'mwcNoPol'
        off = 0;
        p0(4) = off; lb(4) = off; ub(4) = off;
    otherwise
        error('no valid modeltype')
end

model = dorsalFitFunction(modelType);

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