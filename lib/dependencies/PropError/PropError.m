function [value, sigma] = PropError(f,varlist,vals,errs)
%[VALUE, SIGMA] = PROPERROR(F,VARLIST,VALS,ERRS)
%
%Finds the propagated uncertainty in a function f with estimated variables
%"vals" with corresponding uncertainties "errs".
%
%varlist is a row vector of variable names. Enter in the estimated values
%in "vals" and their associated errors in "errs" at positions corresponding 
%to the order you typed in the variables in varlist.
%
%Example using period of a simple harmonic pendulum:
%
%For this example, lets say the pendulum length is 10m with an uncertainty
%of 1mm, and no error in g.
%syms L g
%T = 2*pi*sqrt(L/g)
%type the function T = 2*pi*sqrt(L/g)
%
%PropError(T,[L g],[10 9.81],[0.001 0])
%ans =
%
%    [       6.3437]    '+/-'    [3.1719e-004]
%    'Percent Error'    '+/-'    [     0.0050]
%
%(c) Brad Ridder 2007. Feel free to use this under the BSD guidelines. If
%you wish to add to this program, just leave my name and add yours to it.

%Modified 4/19/2015 by Armando Reimer -> Now this function returns the
%function value and its error. eg: [3.1719e-004, .05] from the previous
%example
n = numel(varlist);
sig = vpa(ones(1,n));
for i = 1:n
    sig(i) = diff(f,varlist(i),1);
end
value = double(subs(f,varlist,vals));
if any(isnan(errs)) == 1
    sigma = NaN;
    warning('Input errors contained NaNs');
else
    sigma = double(sqrt( sum((subs(sig,varlist,vals).^2).*errs.^2)  ));
end