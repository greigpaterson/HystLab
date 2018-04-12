function p = Get_fcdf(x, v1, v2)
%
% Function to generate the cumulative density function (CDF) of the F-distribution.
%
% Inputs:
%        x - x points to evaluate the CDF
%        v1, v2 - degrees of freedom
%
% Outputs:
%        p - the CDF evaluated at the points in x
%
%
% Written by Greig A. Paterson
%
% Updated: 17 May, 2016
%

%% Check inputs

if nargin < 3
    error('Get_fcdf:TooFewInputs', 'At least 3 inputs are required.');
end

if any(size(x) == 1) ~=1
    error('Get_fcdf:Input', 'Data input should be a vector.');
end

% Transpose to a row vector
if size(x,2) ==1
    x = x';
end

if any([size(v1), size(v2)] ~= 1)
    error('Get_fcdf:DoF', 'Degrees of freedom should be scalars.');
end

% Check for negatives, NaNs, and non-integers
if any([v1, v2] <= 0) || any(isnan([v1, v2]) == 1) || rem(v1,1) || rem(v2,1)
    error('Get_fcdf:DoF', 'Degrees of freedom should be integers greater than zero.');
end

%% Get the CDF

% Initialize the p values
% This is valid for x <= 0
p = zeros(size(x));

% Set NaNs in x to NaNs in p
p(isnan(x)) = NaN;

% Set infs in x to ones in p
p(isinf(x)) = 1;

% Identify the indices where x predefines the p values
Defined_Inds = isnan(x) | isinf(x) | (x <= 0);

% Transform the remaining x
xprime = (v1.*x(~Defined_Inds)) ./ (v1.*x(~Defined_Inds) + v2);

% Get the remaining p values
p(~Defined_Inds) = betainc(xprime, v1/2, v2/2);

