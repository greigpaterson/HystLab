function [Processed_Data, Uncorrected_Data, Noise_Data, Fitted_Data, Data_Parameters, Processing_Parameters, Basis_Coeffs, Error_Flag] = Process_Hyst_Data(Data, Order, varargin)
%
% Function that processes the hysteresis data
%
% Input:
%       Data - cell containing the data to porcess [nData x 1] {nFields x 2}
%
%       Order - vector of flags for the data measurement order
%               -1 - Negative to positive fields
%                1 - Positive to negative fields
%
%       varargin - name/value paired data input
%
%       Drift - the flag for the drift correction
%                    0 - No correction
%                    1 - Automatic (DEFAULT)
%                    2 - Positive field correction [1]
%                    3 - Upper branch only
%                    4 - Symmetric averaging and tip-to-tip closure [3]
%                    5 - Paramagnetic correction [4]
%
%       Saturation - the flag for the high-field saturation correction
%                         0 - No correction
%                         1 - Automatic correction (DEFAULT) [1,2]
%                         2 - Linear correction
%                         3 - Approach to saturation [2]
%
%       SaturationField - the field in mT, above which saturation
%                         correction is applied
%
%       FixedBeta - Value for constant approach to saturation coefficient
%                   Set to empty [] for non-fixed values.
%
%       Trim - the flag for trimming off high field data
%                   0 - Do not trim the data
%                   1 - Trim the data
%
%       TrimField - the field (in mT) above which data are trimmed
%
%       PoleSaturation - flag for pole saturation correction - NOT FULLY SUPPORTED!!
%                    0 - No correction (DEFAULT)
%                    1 - Truncate data
%                    2 - Full correction
%
%      PoleData - the data for pole saturation.
%                 PoleFlag = 0, then [];
%                 PoleFlag = 1, then field value in mT
%                 PoleFlag = 2, then
%
%
% Output:
%        Processed_Data - the measurement data adjusted for offset, placed
%                         on a regular field grid, with all applied
%                         corrections [Fields(U,L), Upper, Lower, Mih, Mrh]
%
%        Uncorrected_Data - The minimally processed measurement data with upper and 
%                           lower branches on the same field grid (interpolated to either
%                           the upper or lower branch fields depending on which does 
%                           not add points to the data
%                           [Fields(U,L), Top_curve, Bottom_Curve, Mih, Mrh]
%
%        Fitted_Data - cell containing the fitted hysteresis loop
%                         [Fields(U,L),Upper, Lower, Mih, Mrh]
%
%        Noise_Data - cell with the noise data [Fields, Noise]
%
%        Data_Parameters - array of the hysteresis loop and field parameters
%
%        Processing_Parameters - vector of the parameters used to process the data
%
%        Basis_Coeffs - cell array of the coefficients for the basis
%                       functions used to fit the data
%
%        Error_Flag - error flag for checking the behavior of the automated
%                     saturation correction
%
%
% References:
%
% [1] Jackson, M., and P. Solheid (2010), On the quantitative analysis and
%     evaluation of magnetic hysteresis data, Geochem. Geophys. Geosyst.,
%     11, Q04Z15, doi:10.1029/2009GC002932.
%
% [2] Fabian K. (2006), Approach to saturation analysis of hysteresis
%     measurements in rock magnetism and evidence for stress dominated
%     magnetic anisotropy in young mid-ocean ridge basalt. Phys. Earth
%     Planet. Inter., 154, 299-307, doi: 10.1016/j.pepi.2005.06.016.
%
% [3] von Dobeneck, T. (1996), A systematic analysis of natural magnetic mineral
%     assemblages based on modelling hysteresis loops with coercivity-related
%     hyperbolic basis functions, Geophys. J. Int., 124, 675?694,
%     doi:10.1111/j.1365-246X.1996.tb05632.x
%
% [4] Paterson, G. A., X. Zhao, M. Jackson, D. Heslop, Measuring, processing,
%     and analyzing hysteresis data, Geochemistry, Geophysics, Geosystems, 19,
%     doi: 10.1029/2018GC007620
%
%

% TODO - Full pole saturation. Create GUI for pole saturation data input,
% which runs on first input, but is stored for any other reanalysis
%

% Last Modified 2021/06/011
%

%% Supress unwanted code analyzer warnings

% Unused variables
%#ok<*ASGLU>

%% Process the inputs

if nargin < 2
    error('Process_Hyst_Data:Input', 'At least two input are required');
end

% Set the defaults for the processing

Offset_Flag = 1;
Drift_Flag = 1; %
Saturation_Flag = 1; % Automatic correction
Saturation_Field = [];
FixedBeta_Flag = 0;
FixedBeta_Val = -1.5;

Trim_Flag = 0;
Trim_Field = NaN;

% Internally defined flag to regularize the fields to a regular grid (= 1),
% or to interpolate lower branch field to upper branch fileds (~= 1)
Regularize_Fields_Flag = 0;

Pole_Flag = 0; % No correction
Pole_Data = [];

Fit_Param_Flag = 0;

if nargin > 2
    prescription = varargin;
    
    nt=length(prescription);
    n=nt/2;
    
    if floor(n) ~= n
        error('Process_Hyst_Data:Input', 'Input arguments must have a corresponding value.')
    end
    
    for ii=1:n
        idx=(ii-1)*2+1;
        prop=lower(prescription{idx});
        val=lower(prescription{idx+1});
        
        if strcmpi(prop, 'Offset')
            Offset_Flag = val;
        elseif strcmpi(prop, 'Drift')
            Drift_Flag = val;
        elseif strcmpi(prop, 'Saturation')
            Saturation_Flag = val;
        elseif strcmpi(prop, 'SaturationField')
            Saturation_Field = val;
        elseif strcmpi(prop, 'PoleSaturation')
            Pole_Flag = val;
        elseif strcmpi(prop, 'PoleData')
            Pole_Data = val;
        elseif strcmpi(prop, 'Trim')
            Trim_Flag = val;
        elseif strcmpi(prop, 'TrimField')
            Trim_Field = val;
            if isempty(Trim_Field)
                Trim_Field = NaN;
            end
        elseif strcmpi(prop, 'FitEst')
            Fit_Param_Flag = val;
        elseif strcmpi(prop, 'FixedBeta_Flag')
            FixedBeta_Flag = val;
        elseif strcmpi(prop, 'FixedBeta_Val')
            FixedBeta_Val = val;
        else
            error('Process_Hyst_Data:Input', 'Unrecognized input argument.')
        end
        
    end
    
end

% Override pole saturation input to prevet its use - not fully implemented
% Pole_Flag = 0;
% Catch if no pole saturation data are provided
if Pole_Flag ~= 0 && isempty(Pole_Data)
    %     warning('Process_Hyst_Data:PoleData', 'No pole saturation correction data provided, no correction will be applied.')
    %     Pole_Flag = 0;
    % TODO - handle pole saturation data
end


%% Pre-allocation

Error_Flag = 0;

nData = length(Data);

% Preallocate some data
B0 = NaN(nData, 1);
M0 = NaN(nData, 1);
Q = NaN(nData, 1);
Qf = NaN(nData, 1);
Qrh = NaN(nData, 1);
Qih = NaN(nData, 1);
Noise_RMS = NaN(nData, 1);

Linear_F = NaN(nData,1);
Linear_p = NaN(nData,1);

Nfit = NaN(nData, 1);
Fit_F = NaN(nData, 1);
Fit_p = NaN(nData, 1);
Fit_RMS = NaN(nData, 1);


Ms = NaN(nData,1);
Mrs = NaN(nData,1);
Sat_pct = NaN(nData,1);
Bc = NaN(nData,1);
S_star = NaN(nData,1);
Brh = NaN(nData,1);
Bih = NaN(nData,1);
Bsat = NaN(nData,1);
Xhf = NaN(nData,1);
X0 = NaN(nData,1);
Shape = NaN(nData,1);
C_err = NaN(nData,1);

AS_alpha = NaN(nData,1);
AS_beta = NaN(nData,1);
AS_pVals1 = NaN(nData, 4);
AS_FVals1 = NaN(nData, 4);
AS_pVals2 = NaN(nData, 4);
AS_FVals2 = NaN(nData, 4);

RH_SNR_HN = NaN(nData, 4); % High-field to noise SNR
RH_SNR_HL = NaN(nData, 4); % High-field to low-field SNR

Fitted_Data = cell(nData, 1);
Processed_Data = cell(nData, 1);
Uncorrected_Data = cell(nData, 1);
Noise_Data = cell(nData, 1);
Processing_Parameters = NaN(nData, 14);

Basis_Coeffs = cell(nData, 2);


%% Loop through the files

for ii = 1:1:nData
    
    
    %% Get the relevant specimen data
    Spec_Data = cell2mat(Data(ii));
    Fields = Spec_Data(:,1);
    Moments = Spec_Data(:,2);
    
    if Order(ii) == -1
        % Flip the direction of the fields and moments
        Fields = -Fields;
        Moments = -Moments;
    end
    
    
    %% Trim the fields
    
    if Trim_Flag == 1
        
        if ~isnan(Trim_Field)
            Moments(abs(Fields) > Trim_Field) = [];
            Fields(abs(Fields) > Trim_Field) = [];
        else
            warning('Process_Hyst_Data:TrimField', 'Field trimming has be requested, but no valid field supplied');
        end
        
    else
        Trim_Field = max(abs(Fields));
    end
    
    % The number of data points
    %npts = length(Fields);
    
    
    %% Do the pole saturation correction
    % This is not currently accessible in the user interface
    % The only method that is currently implemented is data truncation
    % This is a hard coded version of the correction of [1], which may not
    % be suitable for all machines
    
    a = [0.7957; 0.6955; -0.8908; 0.5091; -0.1097];
    PS_Field = 1000;
    
    switch Pole_Flag
        case 0 % no correction
            % Do nothing
            
        case 1 % full correction
            
            % Do the posistive fields
            x = Fields(Fields >= PS_Field)./1e3; %Convert back to T for nicely valued coefficents
            Moments(Fields >= PS_Field) = Moments(Fields >= PS_Field) ./ ( a(1) + a(2).*x + a(3).*x.^2 + a(4).*x.^3 + a(5).*x.^4  );
            
            % Do the negative fields
            x = Fields(Fields <= -PS_Field)./1e3; %Convert back to T for nicely valued coefficents
            Moments(Fields <= -PS_Field) = Moments(Fields <= -PS_Field) ./ ( a(1) - a(2).*x + a(3).*x.^2 - a(4).*x.^3 + a(5).*x.^4  );
            
        otherwise
            error('Process_Hyst_Data:PoleSaturationCorrection', 'Unrecognized pole saturation correction type.')
    end
    
    
    %% Split the loop into upper and lower branches
    
    [Top_Curve, Bot_Curve] = Split_Hyst_Loop(Fields, Moments);
    

    %% Get the minimally processed data
    % Used for the "raw" Mih and Mrh curves
    % Known as "Uncorrected_Data" for legacy reasons
    % Uncorrected_Data - The minimally processed measurement data with upper and
    %                    lower branches on the same field grid (interpolated to either
    %                    the upper or lower branch fields depending on which does
    %                    not add points to the data
    %                    [Fields(U,L), Moments(U,L), Mih, Mrh, Noise]
    %
    
    % Regularize the fields
    [Raw_Field_Grid, Raw_Moment_Grid] = Regularize_Fields(Top_Curve, Bot_Curve, 0);
    
    Uncorrected_Data(ii) = {[Raw_Field_Grid, Raw_Moment_Grid, (Raw_Moment_Grid(:,1) + flipud(Raw_Moment_Grid(:,2)))/2,...
        (Raw_Moment_Grid(:,1) - flipud(Raw_Moment_Grid(:,2)))/2, Raw_Moment_Grid(:,1) + Raw_Moment_Grid(:,2)]};
    

    %% Correct for offset
    
    % Get the top curve moment on unique fields
    Top_interp = []; % Clear from previous iteration
    Top_interp(:,1) = sort(unique(Top_Curve(:,1)), 'descend');
    
    if length(Top_Curve(:,1)) ~= length(unique(Top_Curve(:,1)))
        % Top curve has one or more duplicate fields
        % This takes the average moment of the duplicate fields
        Top_interp(:,2) = Interpolate_To_Field(Top_Curve(:,1), Top_Curve(:,2), Top_interp(:,1), 'linear', 0);
    else
        Top_interp(:,2) = Top_Curve(:,2);
    end
    
    % Get the center of symmetry
    % fminsearchbnd uses Nelder-Mead Simplex algorithm
    B0(ii) = fminsearchbnd(@(x) -Center_Hyst_R2(Top_interp, Bot_Curve, x), 0, [], []);
    [R2_B0, M0(ii)] = Center_Hyst_R2(Top_interp, Bot_Curve, B0(ii));
    
    % fminbnd is based on Brent's method, but requires the (not free) Optimization Toolbox
    %     B0(ii) = fminbnd(@(x) -Center_Hyst_R2(Top_interp, Bot_Curve, x), -Top_interp(1,1)/2, Top_interp(1,1)/2);
    %     [R2_H0, M0(ii)] = Center_Hyst_R2(Top_interp, Bot_Curve, B0(ii));
    
    
    % Get the quality of the loop - Here for testing only
    %     R2_Q = GetR2(Top_interp(:,2), Binv_interp);
    %     Q(ii) = log10(1/ sqrt(1-R2_H0) );
    
    % Check if we want to correct it or not
    switch Offset_Flag
        
        case 0
            % No offset correction
            OFC_Top = [Top_Curve(:,1), Top_Curve(:,2)];
            OFC_Bot = [Bot_Curve(:,1) , Bot_Curve(:,2)];
            
        case 1
            % Correct the offset
            OFC_Top = [Top_Curve(:,1) - B0(ii), Top_Curve(:,2) - M0(ii)];
            OFC_Bot = [Bot_Curve(:,1) - B0(ii), Bot_Curve(:,2) - M0(ii)];
            
        otherwise
            error('Process_Hyst_Data:OffsetCorrection', 'Unrecognized offset correction flag.')
    end
    
    % OFC_Top and OFC_Bot are the offset corrected data. This includes zero
    % offset correction results.
    
    %% Test for whole loop linearity
    
    % invert the bottom curve
    Binv = [-OFC_Bot(:,1), -OFC_Bot(:,2)];
    
    % Get the top curve moment on unique fields
    Top_interp = []; % Clear from previous iteration
    Top_interp(:,1) = sort(unique(OFC_Top(:,1)), 'descend');
    
    if length(OFC_Top(:,1)) ~= length(unique(OFC_Top(:,1)))
        % Top curve has one or more duplicate fields
        % This takes the average moment of the duplicate fields
        Top_interp(:,2) = Interpolate_To_Field(OFC_Top(:,1), OFC_Top(:,2), Top_interp(:,1), 'linear', 0);
    else
        Top_interp(:,2) = OFC_Top(:,2);
    end
    
    % Interpolate the inverted bottom curve to the top curve fields
    Binv_interp = Interpolate_To_Field(Binv(:,1), Binv(:,2), Top_interp(:,1), 'linear', 0);
    
    tmp_Fields = [Top_interp(:,1); -Top_interp(:,1)];
    tmp_Moments = [Top_interp(:,2); -Binv_interp];
    
    Bad_Inds = isnan(tmp_Moments); % values where there is no interpolated value
    Bad_Fields = abs(tmp_Fields(Bad_Inds)); % Get the inverse field to remove (this ensures equal numbers of data in the top and bottom curves)
    Bad_Inds = isnan(tmp_Moments) | ismember(abs(tmp_Fields), Bad_Fields);
    
    if rem(sum(Bad_Inds),2)
        % Odd number - should not be here
        error('Process_Hyst_Data:Interpolation', 'An unexpected interpolation error occurred. Please contact Greig Paterson.')
    end
    
    tmp_Fields(Bad_Inds) = [];
    tmp_Moments(Bad_Inds) = [];
    
    % Invert the lower branch
    npts = length(tmp_Fields);
    
    tmp_Fields((npts/2+1):end) = -tmp_Fields((npts/2+1):end);
    tmp_Moments((npts/2+1):end) = -tmp_Moments((npts/2+1):end);
    
    % Get the linear fit
    LinFit = polyfit(tmp_Fields,  tmp_Moments, 1);
    Mhat = LinFit(1).*tmp_Fields + LinFit(2);
    
    % Do the F-Lack-of-Fit test
    [Linear_p(ii), Linear_F(ii)] = Get_F_Test(tmp_Moments, Mhat, 2, 1);
    
    %% Get the Q value for the loop
    
    % invert the bottom curve
    Binv = [-OFC_Bot(:,1), -OFC_Bot(:,2)];
    
    % Get the top curve moment on unique fields
    Top_interp(:,1) = sort(unique(OFC_Top(:,1)), 'descend');
    if length(OFC_Top(:,1)) ~= length(unique(OFC_Top(:,1)))
        % Top curve has one or more duplicate fields
        % This takes the average moment of the duplicate fields
        Top_interp(:,2) = Interpolate_To_Field(OFC_Top(:,1), OFC_Top(:,2), Top_interp(:,1), 'linear', 0);
    else
        Top_interp(:,2) = OFC_Top(:,2);
    end
    
    % Interpolate the inverted bottom curve to the top curve fields
    Binv_interp = Interpolate_To_Field(Binv(:,1), Binv(:,2), Top_interp(:,1), 'linear', 0);
    
    % Remove NaNs - used if no extrapolation is used for interpolation
    XY = [Top_interp(:,2), Binv_interp];
    XY(sum(isnan(XY),2)>0,:) =[];
    
    R2_Q = GetR2(XY(:,1), XY(:,2));
    Q(ii) = log10(1/ sqrt(1-R2_Q) );
    
    
    %% Regularize the data
    
    [Field_Grid, Moment_Grid] = Regularize_Fields(OFC_Top, OFC_Bot, Regularize_Fields_Flag);
        
    Max_Field = max(abs(Field_Grid(:,1)));

    %% Do the drift correction
    
    % Check for exponetial drift correction
    % Need to have an estiamte for paramagnetic susceptibility
    
    if Drift_Flag == 5
        
        % Get an estimate of the high-field (paramagnetic) susceptibility
        % Get this from the high-field slope correction after applying an
        % automated correction
        
        % The main input parameters
        if isempty(Saturation_Field)
            Dum_Input_Params(1) = NaN;
        else
            Dum_Input_Params(1) = Saturation_Field;
        end
        
        if Saturation_Flag == 0
            Dum_Input_Params(2) = 1;
        else
            Dum_Input_Params(2) = Saturation_Flag;
        end
        
        Dum_Input_Params(3) = Max_Field;
        Dum_Input_Params(4) = FixedBeta_Flag;
        Dum_Input_Params(5) = FixedBeta_Val;
        Dum_Input_Params(6) = Linear_p(ii);
        
        [Dum_Moment_Grid, Dum_SC_Output1] = Slope_Correction(Field_Grid, Moment_Grid, Dum_Input_Params);
        
        Chi_Hat = Dum_SC_Output1(3);
        
        clear vars Dum_*;
        
        
    else
        Chi_Hat = [];
    end
    
    % TODO - Check for return correct pulldown menu in main window
    % Do the dirft correction
    if Drift_Flag == 5 && isnan(Chi_Hat)
        warndlg([{'Paramagnetic correction cannot be applied.'},...
            {'Applying automatic correction.'}], 'Paramagnetic Correction', 'modal');
        [Moment_Grid, Drift_Type, Drift_Ratio, DummyVar, DummyVar, Temp_Ratio] = Hyst_Drift_Correction(Field_Grid, Moment_Grid, 10, 1, Chi_Hat);
    else
        [Moment_Grid, Drift_Type, Drift_Ratio, DummyVar, DummyVar, Temp_Ratio] = Hyst_Drift_Correction(Field_Grid, Moment_Grid, 10, Drift_Flag, Chi_Hat);
    end
    
    % Drift_Type
    Drift_Type = Drift_Type +1; % Add one for indexing purposes
    
    
    %% High-Field Slope Correction
    
    % The main input parameters
    if isempty(Saturation_Field)
        Input_Params(1) = NaN;
    else
        Input_Params(1) = Saturation_Field;
    end
    Input_Params(2) = Saturation_Flag;
    Input_Params(3) = Max_Field;
    Input_Params(4) = FixedBeta_Flag;
    Input_Params(5) = FixedBeta_Val;
    Input_Params(6) = Linear_p(ii);
    
    [Moment_Grid, SC_Output1, SC_Output2, SC_Error_Flag] = Slope_Correction(Field_Grid, Moment_Grid, Input_Params);
    
    % Check for error output
    switch SC_Error_Flag
        
        case 1
            
            MSG = {'Insufficient high-field data for slope correction above user selected field value.'; 'At least 3 data points per high-field segement are needed (a total of 12 data).';...
                'The results may be poor, try again with more data points.'};
            warndlg(MSG, 'Insufficient Data');
            
            Error_Flag = 1;
            
            Processed_Data = [];
            Uncorrected_Data = [];
            Noise_Data = [];
            Fitted_Data = [];
            Data_Parameters = [];
            Processing_Parameters = [];
            Basis_Coeffs = [];
            return;
            
        case 2
            MSG = {'Approach to saturation calculation cannot be applied.'; 'The correction is poorly conditioned.';...
                'Please select more high-field data points or apply another correction.'};
            warndlg(MSG, 'Approach to Saturation Correction');
            
            Error_Flag = 1;
            
            Processed_Data = [];
            Uncorrected_Data = [];
            Noise_Data = [];
            Fitted_Data = [];
            Data_Parameters = [];
            Processing_Parameters = [];
            Basis_Coeffs = [];
            return;
            
    end
    
    
    % Get the output data
    Saturation_Method = SC_Output1(1);
    Saturation_Field = SC_Output1(2);
    Xhf(ii) = SC_Output1(3);
    tmp_Ms = SC_Output1(4);
    AS_alpha(ii) = SC_Output1(5);
    AS_beta(ii) = SC_Output1(6);
    Bsat = SC_Output1(7);
    
    AS_pVals1(ii,:) = SC_Output2{1};
    AS_FVals1(ii,:) = SC_Output2{2};
    AS_pVals2(ii,:) = SC_Output2{3};
    AS_FVals2(ii,:) = SC_Output2{4};
    RH_SNR_HN(ii,:) = SC_Output2{5};
    RH_SNR_HL(ii,:) = SC_Output2{6};
    
    
    %% Fit/Filter the data
    
    % Get the fit to the data
    [Fitting_Results, Basis_Coeffs(ii,:), Nfit(ii), Fit_F(ii), Fit_p(ii), Fit_RMS(ii)] = Fit_Hyst_Data(Moment_Grid, Field_Grid(:,1));
    
    % redo the slope correction if needed
    if Fit_Param_Flag == 1
        
        % The main input parameters
        if isempty(Saturation_Field)
            Input_Params(1) = NaN;
        else
            Input_Params(1) = Saturation_Field;
        end
        Input_Params(2) = Saturation_Flag;
        Input_Params(3) = Max_Field;
        Input_Params(4) = FixedBeta_Flag(ii);
        Input_Params(5) = FixedBeta_Val(ii);
        Input_Params(6) = Linear_p(ii);
        
        tmp_Noise = [Field_Grid(:,1), (Moment_Grid(:,1) + Moment_Grid(:,2)) ];
        
        [Fitting_Results(:,3:4), SC_Output1, SC_Output2, SC_Error_Flag] = Slope_Correction(Fitting_Results(:,1:2), Fitting_Results(:,3:4), Input_Params, 2, tmp_Noise);
        
        % Check for error output
        if SC_Error_Flag == 1
            
            MSG = {'Insufficient data provided for approach to saturation calculation.'; ' At least 4 data points per high-field segement are needed (a total of 16 data).';...
                'The results may be poor, try again with more data points.'};
            warndlg(MSG, 'Approach to Saturation Data');
            
            Error_Flag = 1;
            
            Processed_Data = [];
            Uncorrected_Data = [];
            Noise_Data = [];
            Fitted_Data = [];
            Data_Parameters = [];
            Processing_Parameters = [];
            Basis_Coeffs = [];
            return;
            
        end
        
        
        % Get the output data
        Saturation_Method = SC_Output1(1);
        Saturation_Field = SC_Output1(2);
        Xhf(ii) = SC_Output1(3);
        tmp_Ms = SC_Output1(4);
        AS_alpha(ii) = SC_Output1(5);
        AS_beta(ii) = SC_Output1(6);
        Bsat = SC_Output1(7);
        
        AS_pVals1(ii,:) = SC_Output2{1};
        AS_FVals1(ii,:) = SC_Output2{2};
        AS_pVals2(ii,:) = SC_Output2{3};
        AS_FVals2(ii,:) = SC_Output2{4};
        RH_SNR_HN(ii,:) = SC_Output2{5};
        RH_SNR_HL(ii,:) = SC_Output2{6};
        
        % Recalculate the Mrh and Mih curves
        Fitting_Results(:,5) = 0.5 .* (Fitting_Results(:,3) + flipud(Fitting_Results(:,4)));
        Fitting_Results(:,6) = 0.5 .* (Fitting_Results(:,3) - flipud(Fitting_Results(:,4)));
        
        
        % Recalculate the corrected data loop
        if Saturation_Method ~= 1
            Moment_Grid = Moment_Grid - Xhf.*Field_Grid;
        end
        
    end
    
    %% Get the hysteresis parameters
    
    if Fit_Param_Flag == 1
        % Get parameters from the fitted data
        % Fitting_Results(:,1) - Fields (+ve to -ve)
        % Fitting_Results(:,2) - Fields (-ve to +ve)
        % Fitting_Results(:,3) - Upper branch
        % Fitting_Results(:,4) - Lower branch
        % Fitting_Results(:,5) - Mih
        % Fitting_Results(:,6) - Mrh
        
        % redo the slope correction
        
        if Saturation_Method ~= 3
            % We haven't got Ms from approach to saturation
            % Take means of the 4 ends
            
            Ms(ii) = mean(max(abs(Fitting_Results(:,3:4))));
            AS_alpha(ii) = NaN;
            AS_beta(ii) = NaN;
            
        else
            % Use approach to saturation Ms
            
            Ms(ii) = tmp_Ms;
            clear var tmp_Ms;
            
        end
        
        % Get Mrs
        Mrs(ii) = mean( [abs(Interpolate_To_Field(Fitting_Results(:,1), Fitting_Results(:,3), 0, 'linear')), abs(Interpolate_To_Field(Fitting_Results(:,2), Fitting_Results(:,4), 0, 'linear'))] );
        
        % Get Bc
        Bc(ii) = mean( [abs(Interpolate_To_Field(Fitting_Results(:,3), Fitting_Results(:,1), 0, 'linear')), abs(Interpolate_To_Field(Fitting_Results(:,4), Fitting_Results(:,2), 0, 'linear')) ]);
        
        % Get S*
        % TODO - check S* calculations
        S_star(ii) = NaN;
        
        
        % Get shape
        Ehys = polyarea([Fitting_Results(:,1); Fitting_Results(:,2)], [Fitting_Results(:,3); Fitting_Results(:,4)]);
        Shape(ii) = log( Ehys/(4*Ms(ii)*Bc(ii)) );
        
        
        % Convert Xhf to m^3
        Xhf(ii) = Xhf(ii) / (1e4/(4*pi));
        
        
        % Get the ih and rh curves
        Mih = Fitting_Results(:,5);
        Mrh = Fitting_Results(:,6);
        
        
        % Get Bih
        try
            Bih(ii) = mean([abs(interp1(Mih./Ms(ii), Fitting_Results(:,1), 0.5)), abs(interp1(Mih./Ms(ii), Fitting_Results(:,1), -0.5))]);
        catch
            % Flat curve can cause issues with interpolation
            % Trim off the extreme tails
            idx1 = find(Mih./Ms(ii) <= 0.9, 1, 'first');
            idx2 = find(Mih <= 0, 1, 'first');
            idx3 = find(Mih./Ms(ii) <= -0.9, 1, 'first');
            
            Bih(ii) = mean([interp1(Mih(idx1:idx2)./Ms(ii), Fitting_Results(idx1:idx2,1), 0.5), abs(interp1(Mih(idx2:idx3)./Ms(ii), Fitting_Results(idx2:idx3,1), -0.5))]);
            
        end
        
        
        % Get Brh
        % Get the zero field point or the one closest to it
        Mid = find(Fitting_Results(:,1)==0);
        if isempty(Mid)
            Mid = floor(length(Mrh)/2);
        end
        
        % Check we are not at a point the is below the half mark
        % This happens if the number of data near zero field is low
        if abs(Mrh(Mid)/Mrs) < 0.5
            if sign(Fitting_Results(Mid,1)) >0
                Mid = Mid+1;
            else
                Mid = Mid-1;
            end
        end
        
        
        try
            Brh(ii) = mean( abs( [interp1(Mrh(1:Mid)./Mrs(ii), Fitting_Results(1:Mid,1), 0.5), interp1(Mrh(Mid:end)./Mrs(ii), Fitting_Results(Mid:end,1), 0.5)] ) );
        catch
            % Excess zeros cause problems with the interpolation
            % Trim off excess zeros
            idx1 = find(Mrh./Mrs(ii) > 0.01, 1, 'first');
            idx2 = find(Mrh./Mrs(ii) > 0.01, 1, 'last');
            
            try
                Brh(ii) = mean([interp1(Mrh(idx1:Mid)./Mrs(ii), Fitting_Results(idx1:Mid,1), 0.5), abs(interp1(Mrh(Mid+1:idx2)./Mrs(ii),Fitting_Results(Mid+1:idx2,1), 0.5))] );
            catch
                % TODO - Better catch return/warning
                Brh(ii) = NaN;
            end
            
        end
        
        % Saturation percentage
        Sat_pct(ii) = 100 .* abs(Fitting_Results(1,3)) / Ms(ii);
        
    else
        % Get parameters from the data
        
        
        if Saturation_Method ~= 3
            % We haven't got Ms from approach to saturation
            
            if Saturation_Method == 1
                % No slope correction
                % take the mean of the peak moments
                Ms(ii) = mean(abs([Moment_Grid(abs(Field_Grid(:,1)) == max(abs(Field_Grid(:,1))),1); Moment_Grid(abs(Field_Grid(:,2)) == max(abs(Field_Grid(:,2))),2)]));
            else
                % Take the average moment from fields >= the saturation field
                Ms(ii) = mean(abs([Moment_Grid(abs(Field_Grid(:,1)) >= Saturation_Field,1); Moment_Grid(abs(Field_Grid(:,2)) >= Saturation_Field,2)]));
            end
            
            AS_alpha(ii) = NaN;
            AS_beta(ii) = NaN;
        else
            % Use approach to saturation Ms
            Ms(ii) = tmp_Ms;
            clear var tmp_Ms;
        end
        
        % Get Mrs
        Mrs(ii) = mean( [abs(Interpolate_To_Field(Field_Grid(:,1), Moment_Grid(:,1), 0, 'linear')), abs(Interpolate_To_Field(Field_Grid(:,2), Moment_Grid(:,2), 0, 'linear'))] );
        
        % Get Bc
        try
            Bc(ii) = mean( [abs(Interpolate_To_Field(Moment_Grid(:,1), Field_Grid(:,1), 0, 'linear')), abs(Interpolate_To_Field(Moment_Grid(:,2), Field_Grid(:,2), 0, 'linear')) ]);
        catch
            
            % Get a series of indices to interpolate from
            idx1a = find(Moment_Grid(:,1) <= Moment_Grid(1,1)./2, 1, 'first');
            idx1b = find(Moment_Grid(:,1) <= -Moment_Grid(1,1)./2, 1, 'first');
                       
            idx2a = find(Moment_Grid(:,2) >= Moment_Grid(1,2)./2, 1, 'first');
            idx2b = find(Moment_Grid(:,2) >= -Moment_Grid(1,2)./2, 1, 'first');
            
            % Catch cases where this fails and just take 2 points
            if idx1a == idx1b
                idx1a = find(Moment_Grid(:,1) > 0, 1, 'last');
                idx1b = find(Moment_Grid(:,1) < 0, 1, 'first');
            end
            
            if idx2a == idx2b
                idx2a = find(Moment_Grid(:,2) < 0, 1, 'last');
                idx2b = find(Moment_Grid(:,2) > 0, 1, 'first');
            end
            
            Bc(ii) = mean( [abs(Interpolate_To_Field(Moment_Grid(idx1a:idx1b,1), Field_Grid(idx1a:idx1b,1), 0, 'linear')), abs(Interpolate_To_Field(Moment_Grid(idx2a:idx2b,2), Field_Grid(idx2a:idx2b,2), 0, 'linear')) ]);
            
        end
        
        
        
        % Get S*
        % TODO - check S* calculations
        S_star(ii) = NaN;
        % Get the top curve and isolate the upper left quadrant
        %     tmp_top = [Field_Grid(:,1), Moment_Grid(:,1)];
        %     tmp_top(tmp_top(:,1) > 0,:) = []; % remove postive fields
        %     tmp_top(tmp_top(:,2) < 0,:) = []; % remove negative moments
        %
        %     B2 = interp1(tmp_top(:,2), tmp_top(:,1), Mrs(ii)/2);
        %     S1 = 2*B2/-Bc(ii) - 1;
        %
        %     keyboard
        %
        %     tmp_bot = [Field_Grid(:,2), Moment_Grid(:,2)];
        %         tmp_bot(tmp_bot(:,1) < 0,:) = []; % remove negative fields
        %     tmp_bot(tmp_bot(:,2) > 0,:) = []; % remove positive moments
        %     B2 = interp1(tmp_bot(:,2), tmp_bot(:,1), -Mrs(ii)/2);
        %     S2 = 2*B2/Bc(ii) -1 ;
        %
        %
        %      S_star(ii) = mean(abs([S1, S2]));
        %
        %
        %     dF = mean(abs(diff(Field_Grid(:,1))));
        %     dM = gradient(Moment_Grid(:,1), dF);
        %     dM_Bc = abs(interp1(Moment_Grid(:,1), dM, 0));
        %     S_star(ii) = 1-(Mrs(ii)./Bc(ii))*(1/dM_Bc);
        
        %     if (S_star(ii) < 0)
        %         keyboard
        %     end
        
        % Get shape
        Ehys = polyarea([Field_Grid(:,1); (Field_Grid(:,2))], [Moment_Grid(:,1); (Moment_Grid(:,2))]);
        Shape(ii) = log( Ehys/(4*Ms(ii)*Bc(ii)) );
        
        
        % Convert Xhf to m^3
        Xhf(ii) = Xhf(ii) / (1e4/(4*pi));
        
        % Get the ih and rh curves
        Mih = 0.5.*(Moment_Grid(:,1) + flipud(Moment_Grid(:,2)));
        Mrh = 0.5.*(Moment_Grid(:,1) - flipud(Moment_Grid(:,2)));
        
        
        % Get Bih and Brh
        try
            Bih(ii) = mean([interp1(Mih./Ms(ii), Field_Grid(:,1), 0.5), abs(interp1(Mih./Ms(ii), Field_Grid(:,1), -0.5))]);
            
        catch
            idx1 = find(Mih./Ms(ii) <= 0.9, 1, 'first');
            idx2 = find(Mih <= 0, 1, 'first');
            idx3 = find(Mih./Ms(ii) <= -0.9, 1, 'first');
            
            Bih(ii) = mean([interp1(Mih(idx1:idx2)./Ms(ii), Field_Grid(idx1:idx2,1), 0.5), abs(interp1(Mih(idx2:idx3)./Ms(ii), Field_Grid(idx2:idx3,1), -0.5))]);
            
        end
        
        
        % Get the zero field point or the one closest to it
        Mid = find(Field_Grid(:,1)==0);
        if isempty(Mid)
            Mid = floor(length(Mrh)/2);
        end
        
        % Check we are not at a point the is below the half mark
        % This happens if the number of data near zero field is low
        if abs(Mrh(Mid)/Mrs) < 0.5
            if sign(Field_Grid(Mid,1)) >0
                Mid = Mid+1;
            else
                Mid = Mid-1;
            end
        end
        
        
        
        try
            Brh(ii) = mean([interp1(Mrh(1:Mid)./Mrs(ii), Field_Grid(1:Mid,1), 0.5), abs(interp1(Mrh(Mid:end)./Mrs(ii), Field_Grid(Mid:end,1), 0.5))] );
        catch
            
            idx1 = find(Mrh./Mrs(ii) > 0.01, 1, 'first');
            idx2 = find(Mrh./Mrs(ii) > 0.01, 1, 'last');
            
            try
                Brh(ii) = mean([interp1(Mrh(idx1:Mid)./Mrs(ii), Fitting_Results(idx1:Mid,1), 0.5), abs(interp1(Mrh(Mid+1:idx2)./Mrs(ii),Fitting_Results(Mid+1:idx2,1), 0.5))] );
            catch
                % TODO - Better catch return/warning
                Brh(ii) = NaN;
            end
            
        end
        
        
        % Saturation percentage
        Sat_pct(ii) = 100 .* mean([mean(abs(Moment_Grid(1,:))) , mean(abs(Moment_Grid(end,:)))]) / Ms(ii);
        
    end
    
    
    %% Do some final processing
    
    % Call drift correction to get the error/noise curve
    [DummyVar, DummyVar, DummyVar, Err_H, Smooth_Err_H] = Hyst_Drift_Correction(Field_Grid, Moment_Grid, 10, 0);
    
    
    % Get the closure error
    %     C_err(ii) = Err_H(Field_Grid(:,1)==max(Field_Grid(:,1))) - Err_H(Field_Grid(:,1)==min(Field_Grid(:,1)));
    C_err(ii) = Moment_Grid(1,1) - Moment_Grid(end,2);
    
    % Get the Q value
    Binv = (-Moment_Grid(:,2)); % Inverted lower branch
    R2 = GetR2(Moment_Grid(:,1), Binv);
    Qf(ii) = log10( 1 / sqrt((1-R2)) );
    
    
    % Get the Q values for the ih and rh curves
    
    % Get the ih and rh curves
    Mih = 0.5.*(Moment_Grid(:,1) + flipud(Moment_Grid(:,2)));
    Mrh = 0.5.*(Moment_Grid(:,1) - flipud(Moment_Grid(:,2)));
    
    tmp_Mih1 = Mih(Field_Grid(:,1) > 0);
    tmp_Mih2 = flipud(-Mih(Field_Grid(:,1) < 0));
    R2 = GetR2(tmp_Mih1, tmp_Mih2);
    Qih(ii) = log10( 1 / sqrt((1-R2)) );
    
    tmp_Mrh1 = Mrh(Field_Grid(:,1) > 0);
    tmp_Mrh2 = flipud(Mrh(Field_Grid(:,1) < 0));
    R2 = GetR2(tmp_Mrh1, tmp_Mrh2);
    Qrh(ii) = log10( 1 / sqrt((1-R2)) );
    
    % Collate the data to return, reordering the data if needed
    Fitted_Data(ii) = {Fitting_Results};
    Noise_Data(ii) = {[Field_Grid(:,1), Err_H, Smooth_Err_H]};
    Processed_Data(ii) = {[Field_Grid, Moment_Grid, Mih, Mrh]};
    
    Noise_RMS(ii) = sqrt(mean(Err_H.^2));
    
    
    Processing_Parameters(ii,:) = [Drift_Flag, Drift_Type, Drift_Ratio, Temp_Ratio, Saturation_Flag, Saturation_Field, Saturation_Method, Pole_Flag, Offset_Flag,...
        Trim_Flag, Trim_Field, Fit_Param_Flag, FixedBeta_Flag, FixedBeta_Val];
    
    
    % Set to empty for default batch processing
    Saturation_Field =[];
    
    
end

% TODO - Add fit RMS
% The hysteresis statistics
%                   1                       7             10   11  12     13     14
Data_Parameters = [B0, M0, Q, Qf, Ms, Mrs, Bc, Brh, Bih, Xhf, X0, Shape, C_err, Nfit,...
    AS_alpha, AS_beta, AS_pVals1, AS_FVals1, S_star, Fit_F, Fit_p, Linear_F, Linear_p, AS_pVals2, AS_FVals2, Qrh, Qih, Bsat, Sat_pct, Noise_RMS, Fit_RMS, RH_SNR_HN, RH_SNR_HL];
%    15                  17-20    21-24	    25      26    27      28        29w        30-33       34-37    38   39   40      41


