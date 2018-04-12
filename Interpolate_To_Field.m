function [Param_Int] = Interpolate_To_Field(Fields, Moments, Field_Int, Meth, Extrap_Flag)
%
% Hysteresis field sweeps for the upper and lower branches should be
% monotonic, but sometimes this not true. To overcome this during
% interpolation we progressivley remove bad points until the data fit the
% expected monotonicity.
%
% Input:
%        F_Data - the field steps of the parameter to be interpolated
%        Param - the parameter to be interpolated
%        F_Int - the fields to interpolate to
%        Meth - the method of interpolation (linear or shape-preserving piecewise cubic)
%        Int_Dir - the direction of the expected fields (< 0 is decreasing, > 0 is increasing
%
% Output:
%        Param_Int - the interpolated parameter
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
        Param_Int = interp1(Fields, Moments, Field_Int, Meth, 'extrap');
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
        Param_Int = interp1(Fields, Moments, Field_Int, Meth, 'extrap');
    else
        Param_Int = interp1(Fields, Moments, Field_Int, Meth);
    end
    
end

