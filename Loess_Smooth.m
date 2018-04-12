function [Sy] = Loess_Smooth(y, span)
%
% Function to loees smooth y using a windong of span points
% This is used to smooth the hysteresis noise curve for drift correction
% The function fits a locally weighted 2nd order polynomial to y using a
% tri-cube weigting function
% Note: this function assumes uniformly spaced points
%
% Inputs:
%         y - vector of data to smooth [n x 1] (n >= 5)
%         n - the number of points in the smoothing window (span >= 3)
%
% Output:
%         Sy - the smoothed version of Y
%

%% Do some input checks

if nargin < 2
    error('Loess_Smooth:Input', 'At least 2 input arguments are requried');
end

[n, l] = size(y);

% Transpose to [n x 1]
if n == 1
    y = y';
    [n, l] = size(y);
end

if sum([n,l] == 1) ~= 1 || n < 5
    error('Loess_Smooth:Data', 'Input data should be a [n x 1] vector where n is >= 5');
end

if max(size(span)) ~= 1 || rem(span,1) || span < 3
    error('Loess_Smooth:Span', 'The span size, n, must be a scalar integer >= 3');
end

% Define the X data
x = (1:1:n)';

% Make the span an odd number
span = 2*floor(span/2) + 1;


%% The main smoothing

half_width = (span-1)/2;

Sy = NaN(n,1);

for ii = 1:n
    
    % Get the current x and distance from all other values
    xi = x(ii);
    d = x - xi;
    
    % get the weights
    if ii <= half_width
        % At the left end
        w = (1 - abs(d./(1/2 + half_width - ii/2 )/2).^3 ).^3;
    elseif n-ii <= half_width
        % At the end
        w = (1 - abs(d./(n/2 - half_width - ii/2)/2).^3 ).^3;
    else
        % in the middle
        w = (1 - abs(d./half_width).^3 ).^3;
    end
    
    
    % find the positive ones
    inds = w > 0;
    
    w = sqrt(w(inds));
    tmp_x = x(inds);
    tmp_y = y(inds);
    
    
    V = [ones(length(tmp_x),1), tmp_x, tmp_x.^2];
    
    V = V .* repmat(w, 1, size(V,2));
    coeffs = V\(w.*tmp_y);
    
    
    Sy(ii) = coeffs(1) + coeffs(2)*x(ii) + coeffs(3)*x(ii).^2;
    
end



