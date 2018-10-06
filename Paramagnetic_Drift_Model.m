function [MisFit, Model_Err_H, Drift_Results] = Paramagnetic_Drift_Model(Noise, Xhf, Tr, k)
%
% Function that evaluate the fit of a paramagnetic drift model to a hysteresis
% noise curve. Used the apply the paramagnetic deift correction of [1]. The
% model is fitted to the smoothed noise curve to minimize the influence of
% high-frequencey noise.
%
% Input:
%           Noise - The measured noise curve [nx3]. Column 1 is the fields,
%                   column 2 the noise curve, column 3 the smoothed noise.
%
%           Xhf - The estimate of teh high-field paramagentic
%                 susceptibility.
%
%           Tr - The temperature ratio of specimne intial to ambient 
%                temperature (T0/TA).
%
%           k - The rate constant for Newtonian cooling.
%
% Output:
%           MisFit - The norm (RMSE) of the difference between the observed
%                    and modeled noise curves.
%
%           Model_Err_H - The modeled noise curve.
%
%           Drit_Results - The paramagnetic model results [nx3]. Column 1
%                          is the fields, column 2 the non-drifted 
%                          paramagneitc loop, column 3 the drifted 
%                          paramagneitc loop
%
% References:
%
% [1] Paterson, G. A., X. Zhao, M. Jackson, D. Heslop, Measuring, processing,
%     and analyzing hysteresis data, Geochemistry, Geophysics, Geosystems, 19,
%     doi: 10.1029/2018GC007620
%

%% Input processing and checking

if nargin < 4
    error('Paramagnetic_Drift_Model:Input', 'At least 4 input arguments are required.');
end

[n1, n2] = size(Noise);

if n1 < 5
      error('Paramagnetic_Drift_Model:Noise', 'At least 5 noise data are required for fitting.');
end

if n2 ~= 3
        error('Paramagnetic_Drift_Model:Noise', 'Noise data input should have 3 column [field, noise, smoothed noise].');
end

if ~isequal(size(Xhf), size(Tr), size(k))
    error('Paramagnetic_Drift_Model:Input', 'Xhf, Tr, k should all be the same size.');
end

if any(size(Xhf)~= 1)
    error('Paramagnetic_Drift_Model:Input', 'Xhf, Tr, k should all be scalar inputs.');
end


%% Find the model

% Double up the fields for both upper and lower branch sweeps
Fields = [Noise(:,1); -Noise(:,1)];
H = Fields./((4*pi)/1e4); % Convert to A/m

% Get the paramgnetic moments
Moments = H.*Xhf;

% Get the number of points for the loop to use as a "time" index
% Assumes that the measurement steps are approximately equally spaced
npts = length(Fields);
idx = (1:npts)';

% Get the temperature ratios during the loop measurement (T0/Ti)
T_ratio = Tr + exp(-k.*idx) - Tr.*exp(-k.*idx);

% Calculate the momenets after thermal drift
Drift_Moments = Moments .* T_ratio;

% Get the indices for splitting the curves
Ind1 = npts/2;
Ind2 = npts/2;

% Split into upper and lower branched
Top_Curve = [Fields(1:Ind1), Drift_Moments(1:Ind1)];
Bot_Curve = [Fields(Ind2+1:end), Drift_Moments(Ind2+1:end)];

% Calculate the outputs
Drift_Results = [Fields, Moments, Drift_Moments];
Model_Err_H = [Top_Curve(:,1), ( Top_Curve(:,2) + Bot_Curve(:,2) )];
MisFit = norm(Noise(:,3)-Model_Err_H(:,2));

