function [R2, M0] = Center_Hyst_R2(Top_Curve, Bot_Curve, Hoff)
%
% Function for determining the
%
% Input:
%       Top_curve - the top curve of the hysteresis loop [Fields, Magnetization]
%       Bot_curve - the bottom curve of the hysteresis loop [Fields, Magnetization]
%       Hoff - the field offset of the hysteresis loop in mT
%
% Output:
%        R2 - the R^2 value between the top curve and the offset corrected
%             inverted bottom curve
%        M0 - the magnetization offset
%

%% Process inputs

if nargin < 3
    error('Center_Hyst_R2:Input', 'Three inputs are required.');
end

%% Main function
% invert bottom curve and offset the fields
Binv = [-Bot_Curve(:,1)+2*Hoff, -Bot_Curve(:,2)];

% Interpolate the inverted bottom curve to the top curve fields
% Keep linear interpolation to avoid smoothing
Binv_int = Interpolate_To_Field(Binv(:,1), Binv(:,2), Top_Curve(:,1), 'linear', 0);


XY = [Top_Curve(:,2), Binv_int];

% Remove NaNs - used if no extrapolation is used for interpolation
XY(sum(isnan(XY),2)>0,:) =[];

% Do the major axis linear fit
[b, a] = MAxisReg(XY(:,1), XY(:,2));

% The magnetization offset
M0 = -a/2;

% Get the correlation and linear fit
R2 = GetR2(XY(:,1), XY(:,2));


function [b, a] = MAxisReg(X, Y)
%
% Function to perform a major axis regression on input data X and Y
%
% Input:
%           X - data vector 1
%           Y - data vector 2
%
% Output:
%           b - slope of relationship
%           a - y-axis intercept
%

n = length(X);

xbar = mean(X);
ybar = mean(Y);


Sxx = sum( (X - xbar).^2 ) / (n-1);
Syy = sum( (Y - ybar).^2 ) / (n-1);
Sxy = sum( (X - xbar) .* (Y - ybar) ) / (n-1);


b = (1/(2*Sxy)) * (Syy - Sxx + sqrt( (Syy-Sxx)^2 + 4*Sxy^2 ) ) ;

a = ybar - b*xbar;
