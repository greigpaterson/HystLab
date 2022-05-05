function [Field_Grid, Moment_Grid] = Regularize_Fields(Upper_Branch, Lower_Branch, Regularize_Flag)
%
% Function put the upper and lower branches onto the same field steps.
% Either interpolate to grid of regularly spaced points or interpolate to the fields
% of the upper branch.
% In both cases we interpolate such that the fields of a single branch
% are symmetric about zero.
%
%
% Input:
%
%        Upper_Branch - Upper branch of the hysteresis loop
%                    [n x 2] matrix [Fields, Moments]
%
%        Lower_Branch - Lower branch of the hysteresis loop
%                    [n x 2] matrix [Fields, Moments]
%
%       Regularize_Flag - Flag for putting the data on regular grid of fields
%                         1 - Grid of even spaced fields
%                         0 - Upper branch fields [DEFAULT]
%
% Output:
%
%        Field_Grid - The hysteresis fields on a regular field grid
%                     [nData x 2]. Column 1 is the upper branch sweeping
%                     positive to negative fields. Column 2 is the lower branch
%                     sweeping negative to postive fields
%
%        Moment_Grid - The hysteresis moments on a regular field grid
%                      [nData x 2]. Column 1 is the upper branch sweeping
%                      positive to negative fields. Column 2 is the lower branch
%                      sweeping negative to postive fields

%        Minimally_Processed_Data - The minimally processed measurement data with upper and
%                           lower branches on the same field grid (interpolated to either
%                           the upper or lower branch fields depending on which does
%                           not add points to the data
%                           [Fields(U,L), Top_curve, Bottom_Curve, Mih, Mrh, Noise]
%
%
% Last Modified 2021/06/11
%

%% Some input processing and defaults

if nargin < 2
    error('Regularize_Fields:Input', 'At least 2 input arguments are required.')
end

if nargin < 3
    Regularize_Flag = 0;
end

% TODO: Exapnd input checks

%% Do the regularization

if Regularize_Flag == 1
    
    % Interpolate onto a regular grid
    % The value to round the field to
    Fstep = mean([abs(diff(Upper_Branch(:,1))); abs(diff(Lower_Branch(:,1)))]);
    Field_Round = RoundField(Fstep);
    
    % Get the minimum field at the 4 field maxima
    % This approach linear avoids extrapolation in the high field region
    Max_Field = min([max(Upper_Branch(:,1)), abs(min(Upper_Branch(:,1))), max(Lower_Branch(:,1)), abs(min(Lower_Branch(:,1)))]);
    Max_Field = floor(Max_Field/Field_Round)*Field_Round; % round down
    
    
    % Down size the number of points to be restricted to the number of
    % real data points measured in the field range of the grid
    % We take this to be the minimum number from the upper or lower
    % branches
    npts_interp = min([sum(abs(Upper_Branch(:,1)<=Max_Field)), sum(abs(Lower_Branch(:,1)<=Max_Field))]);
    
    %     if rem(npts_interp,2) == 0
    %         % Subtract 1 point to make odd
    %         % This ensures a zero field point
    %         npts_interp = npts_interp -1;
    %     end
    
    % Make the field grid
    % Field_Grid = [Upper branch fields (+ to -), Lower branch fields (- to +)];
    Field_Grid = [linspace(Max_Field, -Max_Field, npts_interp)', linspace(-Max_Field, Max_Field, npts_interp)'];
    
    % Interpolate the moments
    % Moment_Grid = [Upper branch moments (+ to -), Lower branch moments (- to +)];
    Moment_Grid = [ Interpolate_To_Field(Upper_Branch(:,1), Upper_Branch(:,2), Field_Grid(:,1), 'linear', 0), ...
        Interpolate_To_Field(Lower_Branch(:,1), Lower_Branch(:,2), Field_Grid(:,2), 'linear', 0)];
    
else
    
    % Interpolate to the upper branch fields of the measured data
    
    % Split the upper branch into positive and neagive fields, invert
    % the negative fields, and average with the positive fields
    % These mean fields are used for interpolation
    tmp_F1 = Upper_Branch(Upper_Branch(:,1)>0,1);
    tmp_F2 = flipud(-Upper_Branch(Upper_Branch(:,1)<0,1));
    nMin = min([length(tmp_F1), length(tmp_F2)]);
%     mean_Fields = mean([tmp_F1(1:nMin), tmp_F2(1:nMin)],2);
    mean_Fields = mean([tmp_F1((length(tmp_F1)-nMin+1):end), tmp_F2((length(tmp_F2)-nMin+1):end)],2);
    
    % Add in the zero field step or lowest absolute field if present, not used,
    % and it doesn't add more data than was measured.
    if any(mean_Fields == min(abs(Upper_Branch(:,1)))) == 0
        % Min field not used
        if min(abs(Upper_Branch(:,1))) == 0
            % Min field is zero
            if 4*length(mean_Fields)+2 <= length(Upper_Branch(:,1))
                % Not adding in more data
                mean_Fields = [mean_Fields; min(abs(Upper_Branch(:,1)))];
            end
        else
            if 4*(length(mean_Fields)+1) <= length(Upper_Branch(:,1))
                % Not adding in more data
                mean_Fields = [mean_Fields; min(abs(Upper_Branch(:,1)))];
            end
        end
    end
    
    % Make the field grid
    % Field_Grid = [Upper branch fields (+ to -), Lower branch fields (- to +)];
    if min(abs(mean_Fields)) == 0
        % Don't over replicate the zero field step if present
        Field_Grid = [[mean_Fields(1:end-1); flipud(-mean_Fields)], [-mean_Fields(1:end-1); flipud(mean_Fields)]];
    else
        Field_Grid = [[mean_Fields; flipud(-mean_Fields)], [-mean_Fields; flipud(mean_Fields)]];
    end
    
    
    % Interpolate the moments
    % Moment_Grid = [Upper branch moments (+ to -), Lower branch moments (- to +)];
    Moment_Grid = [ Interpolate_To_Field(Upper_Branch(:,1), Upper_Branch(:,2), Field_Grid(:,1), 'linear', 1), ...
        Interpolate_To_Field(Lower_Branch(:,1), Lower_Branch(:,2), Field_Grid(:,2), 'linear', 1)];
        
end