function [Fields, Moments, Initial_Mag_Data] = Remove_Initial_Loop(Fields, Moments, Type, File_ID)
%
% Function to remove the initial magentization aquisition present in some data files
% This is only used for reading data files where there is no clear measurement script
% information with which to directly remove the initial magentization.
%
% Identification is done by looking at the field gradients to determine when the
% sign changes, indicating the field has reached maximum and is sweeping down for
% the upper branch.
% NOTE: This assumes loop measurements always start from postive saturation.
%
%
% Input:
%       Fields -  Vector of hysteresis fields
%
%       Moments - Vector of hysteresis moments corresponding to the above fields
%
%       Type - String flag for the data file format
%               -1 - Negative to positive fields
%                1 - Positive to negative fields
%
% Output:
%        Fields -  Vector of hysteresis fields
%
%        Moments - Vector of hysteresis moments corresponding to the above fields
%
%        Initial_Mag_Data - nx2 array for the initial moment data
%                           [Fields, Moments]
%
% Last Modified 2021/06/20
%

%% Some input processing and defaults

if nargin < 2
    error('Remove_Initial_Loop:Input', 'At least 2 input arguments are required.')
end

if nargin < 3
    Type = 'Unknown';
end

if nargin < 4
    File_ID = 'Unknown';
end

% Check sizes
Fs = size(Fields);
Ms = size(Moments);

if ~any(Fs == 1)
    error('Remove_Initial_Loop:Input_Fields', 'Input fields must be a vector.')
end

if ~any(Ms == 1)
    error('Remove_Initial_Loop:Input_Moments', 'Input moments must be a vector.')
end

% Convery to column vectors
if Fs(1) == Fs(2)
    error('Remove_Initial_Loop:Input_Fields', 'Input fields have more than 1 value.')
end

if Ms(1) == Ms(2)
    error('Remove_Initial_Loop:Input_Moments', 'Input moments have more than 1 value.')
end


% Convery to column vectors
if Fs(1) < Fs(2)
    Fields = Fields';
end

if Ms(1) < Ms(2)
    Moments = Moments';
end

%% Main body

% The inital moment data
Initial_Mag_Data = [];

% Check intial field is near peak field
if abs(Fields(1)) >= 0.75*max(abs(Fields))
    
    % First field is large, relative to max field
    % Assume no initial loop, but noisey field steps
    % Do nothing
    return;
    
end

% The field gradient and sign of gradient acceleration
dF = diff(Fields);
ddF = diff(sign(dF));

switch sum(ddF~=0)
    
    case 0
        % All fields are the same
        err_ident = ['Read_Hyst_Files:', Type];
        error(err_ident', 'All fields are idetical. Please please check data file %s.', File_ID);
        
    case 1
        % All is good - do nothing
        
    case 2
        % We have 3 segments
        
        % Check not numerical error
        % Not exhaustive checks
        % TODO - Add more checks
        
        % Otherwise assume we have intial moment curve
        
        idx = find(ddF~=0); % indices of sign change
        
        % Get the intial curve
        Initial_Mag_Data = [Fields(1:idx(1)), Moments(1:idx(1))];
        
        % Remove inital curve
        Fields(1:idx(1)) = [];
        Moments(1:idx(1)) =[];
                
    otherwise
        % Duplicate field measurements are present we land here        
        
        if sum(ddF~=0) > 2
            % Multiple sign changes
            % Either noisy or has replicate fields
            % Assume first change at maximum field is the
            % one we want
            
            idx = find(ddF~=0);
            max_idx = 1 + idx(find(Fields(idx+1)==max(Fields), 1, 'first'));
            
            % Get the intial curve
            Initial_Mag_Data = [Fields(1: max_idx), Moments(1:max_idx)];
            
            % Remove inital curve
            % Keeping the maximum field for the loop
            Fields(1:(max_idx-1)) = [];
            Moments(1:(max_idx-1)) =[];
            
            
        else
            err_ident = ['Read_Hyst_Files:', Type];
            error(err_ident, 'Unrecognized measurement sequence - Cannot isolate initial moment data. Please contact the authors.');
        end
        
end

