function [Chi_HF, Ms, Mhat, alpha, beta, err_flag] = Approach_to_Saturation(HF_Data, Beta_Flag, Bvals, Call_type)
%
% Function to dertmine the hysteresis slope correction for unsaturated
% specimen following [1] and [2].
%
% Input:
%        HF_Data - [n x 2] matrix of high-field data taken from the lower
%                  branch of the ascending hysteresis loop. Column 1 is the
%                  field data, column 2 is the magnerization
%
%        Bvals - Scalar or vector of beta approach to saturation coefficients
%
%       Call_type - Flag for how function is call (0 - manual call, 1 - automatic)
%
% Output:
%         Chi_HF - high field susceptibility
%         Ms - saturation magnetization
%         Mhat - best-fit magnetization prediction
%         alpha - approach to saturation coefficient
%         beta - approach to saturation coefficient (bound to be within [-1 -2])
%         err_flag - flag returned to catch errors
%                    0 - No problems
%                    1 - Not enough data to fit
%
%
% References:
% [1] Fabian K. (2006), Approach to saturation analysis of hysteresis
%     measurements in rock magnetism and evidence for stress dominated
%     magnetic anisotropy in young mid-ocean ridge basalt. Phys. Earth
%     Planet. Inter., 154, 299-307, doi: 10.1016/j.pepi.2005.06.016.
%
% [2] Jackson, M., and P. Solheid (2010), On the quantitative analysis and
%     evaluation of magnetic hysteresis data, Geochem. Geophys. Geosyst.,
%     11, Q04Z15, doi:10.1029/2009GC002932.
%

%% Input checking and processing

if nargin < 2
    error('Approach_to_Saturation:Input', 'At least two inputs are required.');
end

if nargin < 4
    % Assume automatic to suppress modal warnings
    Call_type = 1;
end


if Beta_Flag == 1
    if nargin < 3 || isempty(Bvals)
        error('Approach_to_Saturation:FixedBeta', 'No fixed beta value specified.');
    end
    
    if length(Bvals) > 1
        error('Approach_to_Saturation:FixedBeta', 'Only singel fixed values of beta are supported.');
    end
    
end

[n1, n2] = size(HF_Data);

if all(size(HF_Data) ~= 2)
    error('Approach_to_Saturation:Input', 'HF_Data should be a [n x 2] matrix.');
end

% Flip round if the user enters [2 x n] data
if n1 < n2
    HF_Data = HF_Data';
end

if size(HF_Data, 1) < 16
    
    % Set to NaNs and error_flag=1 and return
    Chi_HF = NaN;
    Ms = NaN;
    Mhat = NaN(size(HF_Data));
    alpha = NaN;
    beta = NaN;
    err_flag = 1;
    
    if Call_type == 0
        % Manual call so return warning
        MSG = {'Insufficient data provided for approach to saturation calculation.'; ' At least 4 data points per high-field segement are needed (a total of 16 data).';...
            'The results may be poor, try again with more data points.'};
        warndlg(MSG, 'Approach to Saturation Data');
    end
    
    return;
    
end

% Define the range of Bvals to explore

if nargin < 2 || Beta_Flag ~= 1 || isempty(Bvals)
    npts = 1e2;
    Bvals = linspace(-1, -2, npts);
else
    npts = length(Bvals);
end


%% Setup the problem

% The function to generatre the H matrix (eqn 19) for eqn 18 in Ref [2]
Hmat = @(B) [HF_Data(:,1).^1, HF_Data(:,1).^0, HF_Data(:,1).^B];


%% Loop through the Bvals and solve the least squares problem

% Pre-allocate
Fit_Params = NaN(3,npts);
MisFit = NaN(npts, 1);

for ii = 1:npts
    Fit_Params(:,ii) = Hmat(Bvals(ii))\HF_Data(:,2);
    MisFit(ii,1) = norm(HF_Data(:,2) - Hmat(Bvals(ii))*Fit_Params(:,ii), 'fro');
end

% Find the parameters that give the smallest misfit
ind = find(MisFit == min(MisFit), 1, 'first');

Chi_HF = Fit_Params(1,ind);
Ms = Fit_Params(2,ind);
alpha= Fit_Params(3,ind);
beta = Bvals(ind);

Params = [Chi_HF; Ms; alpha];

% Get the model fit
Mhat = Hmat(beta)*Params;


err_flag = 0;

