function [Param_Int] = Interpolate_To_Field(Fields, Moments, Field_Int, Meth, Extrap_Flag)
%
% Hysteresis field sweeps for the upper and lower branches should be
% monotonic, but sometimes this not true. To overcome this during
% interpolation we progressivley remove bad points until the data fit the
% expected monotonicity.
%
% Input:
%        Fields - the field steps of the moments to be interpolated
%        Moments - the moments to be interpolated
%        Field_Int - the fields to interpolate to
%        Meth - the method of interpolation (linear or shape-preserving piecewise cubic)
%        Extrap_Flag - the extrapolation flag (=0 no extrapolation; =1 extrapolate where needed)
%
% Output:
%        Param_Int - the interpolated parameter
%
% Last Modified 2019/05/07
%

%%  Input processing

if nargin < 3
    error('Interpolate_To_Field:Input', 'Not enough input arguements');
end

if ~isequal(length(Fields), length(Moments))
    error('Interpolate_To_Field:Input', 'Field data and parameter inputs should be of equal size.');
end

if nargin < 4
    % Default to linear interpolation and no extrapolation
    Meth = 'linear';
    Extrap_Flag = 0;
end

if nargin < 5
    % Default to no extrapolation
    Extrap_Flag = 0;
end


% Accepatable interpoaltion methods
Interpolate_Meths = {'linear', 'pchip'};

if ~any(strcmp(Meth, Interpolate_Meths))
    warning('Interpolate_To_Field:IntrepolationMethod', 'Unrecognized interpolation methods. This shoud not happen. Defaulting to linear interpolation.');
    Meth = 'linear';
end


%% Do the main function

% Save the input data for comparison during debugging
Orig_data = [Fields, Moments];

try
    % Try the raw interpolation first
    if Extrap_Flag == 1
        Param_Int = Extrapolate_Function(Fields, Moments, Field_Int, Meth);
    else
        Param_Int = interp1(Fields, Moments, Field_Int, Meth);
    end
    
catch
    % Take the gradients and find the points with zero gradient
    dF = diff(Fields);
    Bad_Inds = find(dF==0);
    nBad = length(Bad_Inds);
    
    
    % Keep looping until we remove the bad points
    while nBad > 0
        
        % Get the average of the parameters
        tmp_Param = mean([Moments(Bad_Inds), Moments(Bad_Inds+1)],2);
        
        % set the new parameter values
        Moments(Bad_Inds) = tmp_Param;
        
        % Remove the second step
        Moments(Bad_Inds+1) = [];
        Fields(Bad_Inds+1) = [];
        
        % Recalculate the bad indexs
        dF = diff(Fields);
        Bad_Inds = find(dF==0);
        nBad = length(Bad_Inds);
        
    end
    
    % Get the interpolated data
    if Extrap_Flag == 1
        Param_Int = Extrapolate_Function(Fields, Moments, Field_Int, Meth);
    else
        Param_Int = interp1(Fields, Moments, Field_Int, Meth);
    end
    
end



function Param_Int = Extrapolate_Function(Fields, Moments, Field_Int, Meth)
%
% Function to perform interpolation with linear smoothed extrapolation
%
% Input:
%        Fields - the field steps of the moments to be interpolated
%        Moments - the moments to be interpolated
%        Field_Int - the fields to interpolate to
%        Meth - the method of interpolation (linear or shape-preserving piecewise cubic)
%
% Output:
%        Param_Int - the interpolated parameter
%


% % First do normal interpolation
Param_Int = interp1(Fields, Moments, Field_Int, Meth);

% To avoid extrapolation artifacts, we use linear smoothing
% Fit straight line to data above and below 80% of peak field and
% use this fit for extrapolation
Threshold = 0.8;
tmp_data = sortrows([Fields, Moments]);
minF = min(Fields);
maxF = max(Fields);

% Fit and extrapolate the negative field data
if sum(Field_Int < minF) > 0
    inds = tmp_data(:,1) <= Threshold*minF;
    
    % Some sparse data are too sparse above 80% peak so take 5% steps until
    % we have enough data
    while sum(inds) < 2
        Threshold = Threshold - 0.05;
        inds = tmp_data(:,1) <= Threshold*minF;
    end
    
    tmp_neg_fit = polyfit(tmp_data(inds,1), tmp_data(inds,2), 1);
    Param_Int(Field_Int < minF) = Field_Int(Field_Int < minF)*tmp_neg_fit(1) + tmp_neg_fit(2);
end

% Fit and extrapolate the positive field data
if sum(Field_Int > maxF) > 0
    inds = tmp_data(:,1) >= Threshold*maxF;
    
    % Some sparse data are too sparse above 80% peak so take 5% steps until
    % we have enough data
    while sum(inds) < 2
        Threshold = Threshold - 0.05;
        inds = tmp_data(:,1) >= Threshold*maxF;
    end
    
    tmp_pos_fit = polyfit(tmp_data(inds,1), tmp_data(inds,2), 1);
    Param_Int(Field_Int > maxF) = Field_Int(Field_Int > maxF)*tmp_pos_fit(1) + tmp_pos_fit(2);
end

