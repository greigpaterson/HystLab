function [Moment_Grid, SC_Output1, SC_Output2, Error_Flag] = Slope_Correction(Field_Grid, Moment_Grid, Input_Params, Type_Flag, Noise_Curve)
%
% Function to determine the various high-field slope corrections following [1-3].
%
% Input:
%        Field_Grid - The hysteresis fields on a regular field grid
%                     [nData x 2]. Column 1 is the upper branch sweeping
%                     positive to negative fields. Column 2 is the lower branch
%                     sweeping negative to postive fields
%
%        Moment_Grid - The hysteresis moments on a regular field grid
%                      [nData x 2]. Column 1 is the upper branch sweeping
%                      positive to negative fields. Column 2 is the lower branch
%                      sweeping negative to postive fields
%
%        Input_Params - Matrix of correction parameters
%                      (1) - Saturation Field
%                      (2) - Saturation Method Flag
%                      (3) - Max Field of data
%                      (4) - FixedBeta Application Flag
%                      (5) - FixedBeta Valuse
%                      (6) - Linear p - Pvalue of wole loop linearity
%
%        Type_Flag - The type of data being fitted
%                    1 - Noisy measurement data
%                    ~1 - Model fit to a loop
%
%        Noise_Curve - [nDatax2] matrix of the noise curve. This is used if the correction
%                      is being aplied to the fitted data since, by defintiion, the model
%                      noise curve is zero. If data
%
% Output:
%           Moment_Grid - The slope corrected hysteresis moments
%
%           SC_Ouput1 - vector of scalar outputs
%                       (1) -  Saturation method to apply
%                       (2) -  Field above which correction is applied
%                       (3) -  High-field susceptibility
%                       (4) -  Ms esimtate (not necessarily final value)
%                       (5) -  Approach to saturation alpha value
%                       (6) -  Approach to saturation beta value
%                       (7) -  Saturation field estimate
%
%           SC_Ouput2 - cell of vector outputs
%                       (1) -  p-values for F-Lack-of-Fit test
%                       (2) -  F-values for F-Lack-of-Fit test
%                       (3) -  p-values for F-Model Comparison test
%                       (5) -  F-values for F-Model Comparison test
%                       (6) -  p-values for the sign test for zero median high-field Mrh
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
% [3] Paterson, G. A., X. Zhao, M. Jackson, D. Heslop, Measuring, processing,
%     and analyzing hysteresis data, in prep.
%
%% Process the inputs

if nargin < 4
    Type_Flag = 1; % Defalt to data fitting
    Noise_Curve = [Field_Grid(:,1), (Moment_Grid(:,1) + Moment_Grid(:,2)) ];
end

if nargin < 5
    if Type_Flag ==1
        % Data input so calcualte the noise curve
        Noise_Curve = [Field_Grid(:,1), (Moment_Grid(:,1) + Moment_Grid(:,2)) ];
    else
        error('Slope_Correction:Noise', 'Noise curve must input if data are not being corrected.');
    end
end


Saturation_Field = Input_Params(1);
if isnan(Saturation_Field)
    Saturation_Field = [];
end

Saturation_Flag = Input_Params(2);
Max_Field = Input_Params(3);
FixedBeta_Flag = Input_Params(4);
FixedBeta_Val = Input_Params(5);
Linear_p = Input_Params(6);


% Set default values for approach to saturation parameters
tmp_Ms = NaN;
AS_alpha = NaN;
AS_beta = NaN;

% Default error flag
Error_Flag = 0;

% Tolerance for zero moment
% Mtol = sqrt(eps);
% Mtol = eps;

% TODO - Tidy up slope correction code

%% Test the various possible corrections

% Define a % Field Upper Limit (FUL) to exclude some very high fields that
% may still be influenced by pole saturation
FUL = 0.99 * Max_Field;

if isempty(Saturation_Field) || Saturation_Field > Max_Field
    % We are applying saturation correction, but have no field value
    % Default to a correction at 70% of the maximum field
    Saturation_Field = 0.7 * Max_Field;
end


% Get values for the automatic correction
% The percentage values to test the correction
Percentages = [0.7, 0.8, 0.9];
nPercent = length(Percentages);

% Preallocate some variable
% Linear model F-test
FpVals1 = NaN(nPercent, 1); % P vals for F-test
FFVals1 = NaN(nPercent, 1); % F vals for F-test

% Model comparison F-test
FpVals2 = NaN(nPercent, 1); % P vals for F-test
FFVals2 = NaN(nPercent, 1); % F vals for F-test

% Loop closure test
SNR_HF_Noise = NaN(nPercent, 1); % P vals for sign-test
SNR_HF_Area = NaN(nPercent, 1); % S vals for sign-test

% Get the Mrh fields
% Get the averaged Mrh curve by flipping the negative field half to
% positive fields
Mrh = [Field_Grid(:,1), 0.5 .* (Moment_Grid(:,1) - flipud(Moment_Grid(:,2)))];
tmp_var = [Mrh(Mrh(:,1) > 0,:), flipud(Mrh(Mrh(:,1) < 0,2))];
Mrh_Mean = [tmp_var(:,1), mean(tmp_var(:,2:3),2)];
Mrh_Mean(Mrh_Mean(:,2)<0,2) = 0; % Set negatives to zero for postive areas only (negative bits are noise and should be excluded from SNR)
Area_Mrh = abs(trapz(Mrh_Mean(:,1), Mrh_Mean(:,2))); % The area of the total Mrh curve


% Loop through the field percentages
for jj = 1:nPercent
    
    Saturation_Field_AC = Percentages(jj) * Max_Field;
    
    % Get the high field data
    % Take the 4 individual branches
    % Odd numbers are decreasing from saturation
    % Even numbers are approaching saturation
    HF_Data_AC_1 = [Field_Grid(Field_Grid(:,1) >= Saturation_Field_AC, 1), Moment_Grid(Field_Grid(:,1) >= Saturation_Field_AC, 1)];
    HF_Data_AC_2 = [-Field_Grid(Field_Grid(:,1) <= -Saturation_Field_AC, 1), -Moment_Grid(Field_Grid(:,1) <= -Saturation_Field_AC, 1)];
    HF_Data_AC_3 = [-Field_Grid(Field_Grid(:,2) <= -Saturation_Field_AC, 2), -Moment_Grid(Field_Grid(:,2) <= -Saturation_Field_AC, 2)];
    HF_Data_AC_4 = [Field_Grid(Field_Grid(:,2) >= Saturation_Field_AC, 2), Moment_Grid(Field_Grid(:,2) >= Saturation_Field_AC, 2)];
    
    % Compile the data and trim the very high fields
    HF_Data_AC = [HF_Data_AC_1;HF_Data_AC_2;HF_Data_AC_3;HF_Data_AC_4];
    HF_Data_AC(abs(HF_Data_AC(:,1)) > FUL, :) = [];
    
    
    % Test the linear and approach to saturation models
    
    % Get the best linear model
    LinFit = polyfit(HF_Data_AC(:,1), HF_Data_AC(:,2), 1);
    Mhat_lin = LinFit(1).*HF_Data_AC(:,1) + LinFit(2);
    
    % Get the best fit Approach_to_Saturation
    [DummyVar, DummyVar, Mhat_nl, DummyVar, DummyVar, AS_err_flag] = Approach_to_Saturation(HF_Data_AC, FixedBeta_Flag, FixedBeta_Val, 1); %#ok<ASGLU>
    
    
    % Linear model lack-of-fit F-test
    if Type_Flag ==1
        % Real data
        [FpVals1(jj), FFVals1(jj)] =  Get_F_Test(HF_Data_AC(:,2), Mhat_lin, 2, 1);
    end
    
    
    % Check the approach to saturation returned valid data and get the
    % F-test model comparison results
    if AS_err_flag ~= 1
        [FpVals2(jj), FFVals2(jj)] =  Get_F_Test([HF_Data_AC(:,2), HF_Data_AC(:,2)], [Mhat_lin, Mhat_nl], [2,4], 2);
    end
    
    
    
    % Do test for loop closure
    
    % Get the averaged Mrh curve by flipping the negative field half to
    % positive fields
    %     tmp_var = [Mrh(Mrh(:,1) > 0,:), flipud(Mrh(Mrh(:,1) < 0,2))];
    %     Mrh_Mean = [tmp_var(:,1), mean(tmp_var(:,2:3),2)];
    %     Mrh_Mean(Mrh_Mean(:,2)<0,2) = 0; % Set neagtives to zero for postive areas only (negative bits are noise and should be excluded from SNR)
    %     Area_Mrh = abs(trapz(Mrh_Mean(:,1), Mrh_Mean(:,2))); % The area of the total Mrh curve
    
    
    % Get the high field average Mrh section
    HF_Mrh = Mrh_Mean(Mrh_Mean(:,1) >= Saturation_Field_AC, :);
    HF_Mrh(abs(HF_Mrh(:,1)) > FUL, :) = []; % trim the extreme fields
    Area_HF_Mrh = abs(trapz(HF_Mrh(:,1), HF_Mrh(:,2))); % The area of the high field Mrh curve
    RMS_HF_Mrh = sqrt(mean(HF_Mrh(:,2).^2)); % RMS
    
    % Get the high field noise
    % Here we take both branches for calculating the average power
    HF_Noise = Noise_Curve(abs(Noise_Curve(:,1)) >= Saturation_Field_AC, :);
    HF_Noise(abs(HF_Noise(:,1)) > FUL, :) = [];
    RMS_HF_Noise = sqrt(mean(HF_Noise(:,2).^2));
    
    
    % Get the closure stats
    SNR_HF_Area(jj) = 20*log10(Area_HF_Mrh / Area_Mrh );
    SNR_HF_Noise(jj) = 20*log10(RMS_HF_Mrh / RMS_HF_Noise );
    
    
end


% Get the slope correction stats for the input saturation field
Saturation_Field = abs(Saturation_Field);

% Define the high field data for the non-automatic corrections
% Take the 4 individual branches
% Odd numbers are decrease from saturation
% Even numbers are approach to saturation
HF_Data_1 = [Field_Grid(Field_Grid(:,1) >= Saturation_Field, 1), Moment_Grid(Field_Grid(:,1) >= Saturation_Field, 1)];
HF_Data_2 = [-Field_Grid(Field_Grid(:,1) <= -Saturation_Field, 1), -Moment_Grid(Field_Grid(:,1) <= -Saturation_Field, 1)];
HF_Data_3 = [-Field_Grid(Field_Grid(:,2) <= -Saturation_Field, 2), -Moment_Grid(Field_Grid(:,2) <= -Saturation_Field, 2)];
HF_Data_4 = [Field_Grid(Field_Grid(:,2) >= Saturation_Field, 2), Moment_Grid(Field_Grid(:,2) >= Saturation_Field, 2)];

% Compile the data and trim the very high fields
HF_Data = [HF_Data_1;HF_Data_2;HF_Data_3;HF_Data_4];
HF_Data(abs(HF_Data(:,1)) > FUL, :) = [];


% Do test for loop closure

% Get the high field average Mrh section
HF_Mrh = Mrh_Mean(Mrh_Mean(:,1) >= Saturation_Field, :);
HF_Mrh(abs(HF_Mrh(:,1)) > FUL, :) = []; % trim the extreme fields
Area_HF_Mrh = abs(trapz(HF_Mrh(:,1), HF_Mrh(:,2))); % The area of the high field Mrh curve
RMS_HF_Mrh = sqrt(mean(HF_Mrh(:,2).^2)); % RMS

% Get the high field noise
% Here we take both branches for calculating the average power
HF_Noise = Noise_Curve(abs(Noise_Curve(:,1)) >= Saturation_Field, :);
HF_Noise(abs(HF_Noise(:,1)) > FUL, :) = [];
RMS_HF_Noise = sqrt(mean(HF_Noise(:,2).^2));


% Get the closure stats
Current_HF_Area = 20*log10(Area_HF_Mrh / Area_Mrh );
Current_SNR_HF_Noise = 20*log10(RMS_HF_Mrh / RMS_HF_Noise );

RH_SNR_HF_Area = [SNR_HF_Area', Current_HF_Area];
RH_SNR_HF_Noise = [SNR_HF_Noise', Current_SNR_HF_Noise];



% Return a warning if too few data are available
if size(HF_Data, 1) < 16
    
    % Don't need to return a warning since there is no correction or
    % automatic
    if Saturation_Flag == 0 || Saturation_Flag == 1
        
        % Return NaN for currently selected field value
        AS_pVals1 = [FpVals1', NaN];
        AS_FVals1 = [FFVals1', NaN];
        
        AS_pVals2 = [FpVals2', NaN];
        AS_FVals2 = [FFVals2', NaN];
        
    else
        % Return warning
        
        %         MSG = {'Insufficient data provided for approach to saturation calculation.'; ' At least 4 data points per high-field segement are needed (a total of 16 data).';...
        %             'The results may be poor, try again with more data points.'};
        %         warndlg(MSG, 'Approach to Saturation Data');
        
        Error_Flag = 1;
        
        %         Processed_Data = [];
        %         Uncorrected_Data = [];
        %         Noise_Data = [];
        %         Fitted_Data = [];
        %         Data_Parameters = [];
        %         Processing_Parameters = [];
        %         Basis_Coeffs = [];
        return;
        
    end
    
else
    % We have enough data so to the correction
    
    % Get the best linear model
    LinFit = polyfit(HF_Data(:,1), HF_Data(:,2), 1);
    Mhat_lin = LinFit(1).*HF_Data(:,1) + LinFit(2);
    
    % Get the best fit Approach_to_Saturation model
    [DummyVar, DummyVar, Mhat_nl] = Approach_to_Saturation(HF_Data, FixedBeta_Flag, FixedBeta_Val, 1); %#ok<ASGLU>
    
    
    % Linear model lack-of-fit F-tet
    if Type_Flag == 1
        [Current_Field_FpVals, Current_Field_FFVals] =  Get_F_Test(HF_Data(:,2), Mhat_lin, 2, 1);
    else
        Current_Field_FpVals = NaN;
        Current_Field_FFVals = NaN;
    end
    
    % F-test model comparison results
    [Current_Field_FpVals2, Current_Field_FFVals2] =  Get_F_Test([HF_Data(:,2), HF_Data(:,2)], [Mhat_lin, Mhat_nl], [2,4], 2);
    
    
    AS_pVals1 = [FpVals1', Current_Field_FpVals];
    AS_FVals1 = [FFVals1', Current_Field_FFVals];
    
    AS_pVals2 = [FpVals2', Current_Field_FpVals2];
    AS_FVals2 = [FFVals2', Current_Field_FFVals2];
    
end


%% Find the correction to apply
%
% All the above is testing or finding the appropriate correction
% Now we can determine the correction paramaters
%

% Switch through the options
switch Saturation_Flag
    case 0 % no correction
        % Do nothing
        
        Xhf = NaN;
        Saturation_Method = 1;
        
    case 1 % Automatic correction
        
        
        % Matrix for the slope saturation logic tests
        Logic_Inds = NaN(3,3);
        
        % Closed loop?
        % We have a low Mrh SNR or the high field Mrh portion is very small
        Logic_Inds(:,1) = SNR_HF_Noise < 8 | SNR_HF_Area < -48;
        
        % Linear high field? (lack-of-fit result)
        % Lack of fit p-values > 5% indicate that linear fit cannot be rejected
        Logic_Inds(:,2) = FpVals1 > 0.05;
        
        % Linear high field? (model comparison result)
        % Model comparion p-values > 5% indicate that the simpler linear model is better
        Logic_Inds(:,3) = FpVals2 > 0.05;
        
        
        if Linear_p > 0.05
            % F-Test cannot reject the hypothesis that the loop is just
            % a straight line, so do no correction
            
            % Get sucsptibility from the whole loop
            LinFit = polyfit(Field_Grid(:), Moment_Grid(:), 1);
            Xhf = LinFit(1);
            
            Saturation_Method = 1;
            
        elseif sum(Logic_Inds(:,1)) == 0
            % Loop is open at all tested values
            % Do not apply any correction
            
            Xhf = NaN;
            Saturation_Method = 1;
            
        elseif sum(sum(Logic_Inds(:,2:3))) ~= 0
            % Loop is closed, but we cannot reject linearity for at least one part
            % Apply correction to the lowest field with closure where we cannot reject linearity
            
            Min_Ind = find(Logic_Inds(:,1) == 1 & (Logic_Inds(:,2) == 1 | Logic_Inds(:,3) == 1), 1, 'first');
            
            % For noisy data the indices of the logic array (closed,
            % linear, linear) don't match up, yield an empty Min_Ind
            % Simply take the lowest field that supports closure
            if isempty(Min_Ind)
                Min_Ind = find(Logic_Inds(:,1) == 1, 1, 'first');
            end
            
            Saturation_Field = Percentages(Min_Ind) * Max_Field;
            
            % Get the high field data
            HF_Data_1 = [Field_Grid(Field_Grid(:,1) >= Saturation_Field, 1), Moment_Grid(Field_Grid(:,1) >= Saturation_Field, 1)];
            HF_Data_2 = [-Field_Grid(Field_Grid(:,1) <= -Saturation_Field, 1), -Moment_Grid(Field_Grid(:,1) <= -Saturation_Field, 1)];
            HF_Data_3 = [-Field_Grid(Field_Grid(:,2) <= -Saturation_Field, 2), -Moment_Grid(Field_Grid(:,2) <= -Saturation_Field, 2)];
            HF_Data_4 = [Field_Grid(Field_Grid(:,2) >= Saturation_Field, 2), Moment_Grid(Field_Grid(:,2) >= Saturation_Field, 2)];
            
            HF_Data = [HF_Data_1;HF_Data_2;HF_Data_3;HF_Data_4];
            
            % Trim the very high fields
            HF_Data(abs(HF_Data(:,1)) > FUL, :) = [];
            
            
            LinFit = polyfit(HF_Data(:,1), HF_Data(:,2), 1);
            
            Xhf = LinFit(1);
            
            Saturation_Method = 2;
            
            % Update the Fval/Sval table
            AS_pVals1(end) = FpVals1(Min_Ind);
            AS_FVals1(end) = FFVals1(Min_Ind);
            
            AS_pVals2(end) = FpVals2(Min_Ind);
            AS_FVals2(end) = FFVals2(Min_Ind);
            
            
            RH_SNR_HF_Noise(end) = SNR_HF_Noise(Min_Ind);
            RH_SNR_HF_Area(end) = SNR_HF_Area(Min_Ind);
            
            
        elseif sum(sum(Logic_Inds(:,2:3))) == 0
            % Loop is closed, but we reject linearity based on both tests
            % Apply correction to the lowest field with closure where we cannot reject linearity
            
            Min_Ind = find(Logic_Inds(:,1) == 1, 1, 'first');
            Saturation_Field = Percentages(Min_Ind) * Max_Field;
            
            % Get the high field data
            HF_Data_1 = [Field_Grid(Field_Grid(:,1) >= Saturation_Field, 1), Moment_Grid(Field_Grid(:,1) >= Saturation_Field, 1)];
            HF_Data_2 = [-Field_Grid(Field_Grid(:,1) <= -Saturation_Field, 1), -Moment_Grid(Field_Grid(:,1) <= -Saturation_Field, 1)];
            HF_Data_3 = [-Field_Grid(Field_Grid(:,2) <= -Saturation_Field, 2), -Moment_Grid(Field_Grid(:,2) <= -Saturation_Field, 2)];
            HF_Data_4 = [Field_Grid(Field_Grid(:,2) >= Saturation_Field, 2), Moment_Grid(Field_Grid(:,2) >= Saturation_Field, 2)];
            
            HF_Data = [HF_Data_1;HF_Data_2;HF_Data_3;HF_Data_4];
            
            
            % Trim the very high fields
            HF_Data(abs(HF_Data(:,1)) > FUL, :) = [];
            
            
            % Get the approach to saturation fit and parameters
            [Xhf, tmp_Ms, Mhat_nl, AS_alpha, AS_beta] = Approach_to_Saturation(HF_Data, FixedBeta_Flag, FixedBeta_Val, 1); %#ok<ASGLU>
            
            % Set the used method
            Saturation_Method = 3;
            
            % Update the Fval/Sval table
            AS_pVals1(end) = FpVals1(1);
            AS_FVals1(end) = FFVals1(1);
            
            AS_pVals2(end) = FpVals2(1);
            AS_FVals2(end) = FFVals2(1);
            
            RH_SNR_HF_Noise(end) = SNR_HF_Noise(1);
            RH_SNR_HF_Area(end) = SNR_HF_Area(1);
            
        else
            % Shouldn't be here - Logic error... how illogical
            MSG = [{'Automatic correction logic test has failed.'};...
                {'We really should not be here... how illogical.'};...
                {'No slope correction has been applied'}];
            
            warndlg(MSG, 'Logic Failure')
            Xhf = NaN;
            Saturation_Method = 1;
        end
        
    case 2 % Linear correction
        
        Saturation_Method = 2;
        
        % Get the high field data
        HF_Data_1 = [Field_Grid(Field_Grid(:,1) >= Saturation_Field, 1), Moment_Grid(Field_Grid(:,1) >= Saturation_Field, 1)];
        HF_Data_2 = [-Field_Grid(Field_Grid(:,1) <= -Saturation_Field, 1), -Moment_Grid(Field_Grid(:,1) <= -Saturation_Field, 1)];
        HF_Data_3 = [-Field_Grid(Field_Grid(:,2) <= -Saturation_Field, 2), -Moment_Grid(Field_Grid(:,2) <= -Saturation_Field, 2)];
        HF_Data_4 = [Field_Grid(Field_Grid(:,2) >= Saturation_Field, 2), Moment_Grid(Field_Grid(:,2) >= Saturation_Field, 2)];
        
        HF_Data = [HF_Data_1;HF_Data_2;HF_Data_3;HF_Data_4];
        
        % Trim the very high fields
        HF_Data(abs(HF_Data(:,1)) > FUL, :) = [];
        
        LinFit = polyfit(HF_Data(:,1), HF_Data(:,2), 1);
        
        Xhf = LinFit(1);
        
        
    case 3 % Approach to saturation
        
        Saturation_Method = 3;
        
        % Get the high field data
        HF_Data_1 = [Field_Grid(Field_Grid(:,1) >= Saturation_Field, 1), Moment_Grid(Field_Grid(:,1) >= Saturation_Field, 1)];
        HF_Data_2 = [-Field_Grid(Field_Grid(:,1) <= -Saturation_Field, 1), -Moment_Grid(Field_Grid(:,1) <= -Saturation_Field, 1)];
        HF_Data_3 = [-Field_Grid(Field_Grid(:,2) <= -Saturation_Field, 2), -Moment_Grid(Field_Grid(:,2) <= -Saturation_Field, 2)];
        HF_Data_4 = [Field_Grid(Field_Grid(:,2) >= Saturation_Field, 2), Moment_Grid(Field_Grid(:,2) >= Saturation_Field, 2)];
        
        HF_Data = [HF_Data_1;HF_Data_2;HF_Data_3;HF_Data_4];
        
        % Trim the very high fields
        HF_Data(abs(HF_Data(:,1)) > FUL, :) = [];
        
        [Xhf, tmp_Ms, Mhat_nl, AS_alpha, AS_beta, AS_err_flag] = Approach_to_Saturation(HF_Data, FixedBeta_Flag, FixedBeta_Val, 0);
        
        
        if AS_err_flag == 1
            
            % Not enough data
            % The user has been warned
            Error_Flag = 1;
            return;
            
        end
        
    otherwise
        error('Process_Hyst_Data:HighFieldSaturationCorrection', 'Unrecognized high field saturation correction type.')
end


%% Apply saturation correction

% Get the saturation field
switch Saturation_Method
    case 1
        % No correction
        Bsat = NaN;
        
    case 2
        % Linear slope correction
        Bsat = Saturation_Field;
        
        Moment_Grid = Moment_Grid - Xhf.*Field_Grid;
        
    case 3
        % Approach to saturation
        Sat_lim = 0.99; % Take 99% saturation
        
        func = @(x) ((1-Sat_lim) + (AS_alpha .* x.^AS_beta)./tmp_Ms ).^2;
        options = optimset('MaxFun', 1e4, 'MaxIter', 1e4, 'TolX', 1e-6, 'TolFun', 1e-6, 'Display', 'off');
        Bsat = fminsearchbnd(func, Max_Field, min(HF_Data_AC(:,1)), [], options);
        
        Moment_Grid = Moment_Grid - Xhf.*Field_Grid;
        
    otherwise
        error('Process_Hyst_Data:Bsat', 'Unrecognized high field saturation correction method.');
end



%% Collate the output

SC_Output1(1) = Saturation_Method;
SC_Output1(2) =  Saturation_Field;
SC_Output1(3) =  Xhf;
SC_Output1(4) = tmp_Ms;
SC_Output1(5) =  AS_alpha;
SC_Output1(6) = AS_beta;
SC_Output1(7) = Bsat;


SC_Output2{1} = AS_pVals1;
SC_Output2{2} = AS_FVals1;
SC_Output2{3} = AS_pVals2;
SC_Output2{4} = AS_FVals2;
SC_Output2{5} = RH_SNR_HF_Area;
SC_Output2{6} = RH_SNR_HF_Noise;

