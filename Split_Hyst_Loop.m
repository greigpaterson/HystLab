function [Top_Curve, Bot_Curve] = Split_Hyst_Loop(Fields, Moments)
%
% Function to split a hysteresis loop into upper and lower branches.
% Field monoticity is enforced for each branch and duplicate field data are averaged
% within each branch.
% The function also returns the minimally processed data that are needed to
% calculate the Mrh and Mih curves
%
%
% Input:
%       Fields -  Vector of hysteresis fields from postive to negative
%                 saturation
%
%       Moments - Vector of hysteresis moments corresponding to the above fields
%
%       Order - Flags for the data measurement order
%               -1 - Negative to positive fields
%                1 - Positive to negative fields
%
% Output:
%        Top_Curve - Upper branch of the hysteresis loop
%                    [n x 2] matrix [Fields, Moments]
%
%        Bot_Curve - Lower branch of the hysteresis loop
%                    [n x 2] matrix [Fields, Moments]
%
%        Minimally_Processed_Data - The minimally processed measurement data with upper and
%                           lower branches on the same field grid (interpolated to either
%                           the upper or lower branch fields depending on which does
%                           not add points to the data
%                           [Fields(U,L), Top_curve, Bottom_Curve, Mih, Mrh, Noise]
%
%
% Last Modified 2021/06/11
%

%% Supress unwanted code analyzer warnings

% Unused variables
%#ok<*ASGLU>

%% Some input processing and defaults

if nargin < 2
    error('Split_Hyst_Loop:Input', 'At least 2 input arguments are required.')
end

% Check sizes
Fs = size(Fields);
Ms = size(Moments);

if ~any(Fs == 1)
    error('Split_Hyst_Loop:Input_Fields', 'Input fields must be a vector.')
end

if ~any(Ms == 1)
    error('Split_Hyst_Loop:Input_Moments', 'Input moments must be a vector.')
end

% Convery to column vectors
if Fs(1) == Fs(2)
    error('Split_Hyst_Loop:Input_Fields', 'Input fields have more than 1 value.')
end

if Ms(1) == Ms(2)
    error('Split_Hyst_Loop:Input_Moments', 'Input moments have more than 1 value.')
end


% Convery to column vectors
if Fs(1) < Fs(2)
    Fields = Fields';
end

if Ms(1) < Ms(2)
    Moments = Moments';
end


%% Split data based on minimum negative saturation field

Min_field = min(Fields);


Min_idx = find(Fields==Min_field);

switch length(Min_idx)
    case 0
        error('Split_Hyst_Loop:Min_Field', 'No minimum field found')
    case 1
        % All good do nothing
    case 2
        Min_idx = Min_idx(1);
    otherwise
        % split based on odd/even points
        if rem(length(Min_idx),2)
            % Odd
            Min_idx = Min_idx( ceil(length(Min_idx)/2) );
        else
            % even
            Min_idx = Min_idx( length(Min_idx)/2 );
        end
end


Top_Curve = [Fields(1:Min_idx), Moments(1:Min_idx)];
Bot_Curve = [Fields(Min_idx+1:end), Moments(Min_idx+1:end)];



%% Sort the curves by field to ensure monotonic field sweeps
Top_Curve = sortrows(Top_Curve, -1);

Bot_Curve = sortrows(Bot_Curve, 1);

%% Average duplicate field data

% Top curve
if length(Top_Curve(:,1)) ~= length(unique(Top_Curve(:,1)))
    
    % Find indices that are duplicate
    [unique_data, row_idx] = unique(Top_Curve(:,1),'first');
    Duplicate_idx = find(not(ismember(1:numel(Top_Curve(:,1)),row_idx)));
    
    % Loop through and average the moments
    for jj = 1: length(Duplicate_idx)
        if Top_Curve(Duplicate_idx(jj)-1,1) == Top_Curve(Duplicate_idx(jj),1)
            Top_Curve(Duplicate_idx(jj)-1,2) = mean([Top_Curve(Duplicate_idx(jj)-1,2),Top_Curve(Duplicate_idx(jj),2)]);
        end
    end
    
    % Remove the excess data
    Top_Curve(Duplicate_idx,:) = [];
    
end

% Bottom curve
if length(Bot_Curve(:,1)) ~= length(unique(Bot_Curve(:,1)))
    
    % Find indices that are duplicate
    [unique_data, row_idx] = unique(Bot_Curve(:,1),'first');
    Duplicate_idx = find(not(ismember(1:numel(Bot_Curve(:,1)),row_idx)));
    
    % Loop through and average the moments
    for jj = 1: length(Duplicate_idx)
        if Bot_Curve(Duplicate_idx(jj)-1,1) == Bot_Curve(Duplicate_idx(jj),1)
            Bot_Curve(Duplicate_idx(jj)-1,2) = mean([Bot_Curve(Duplicate_idx(jj)-1,2), Bot_Curve(Duplicate_idx(jj),2)]);
        end
    end
    
    % Remove the excess data
    Bot_Curve(Duplicate_idx,:) = [];
    
end

