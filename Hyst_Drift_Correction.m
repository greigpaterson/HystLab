function [Corrected_Moments, Type, Drift_Ratio, Err_H, Smooth_Err_H, Temp_Ratio] = Hyst_Drift_Correction(Field_Grid, Moment_Grid, Field_Span, Type, Chi_hat)
%
% Function that corrects for hysteresis drift using multiple routines. All
% routines use the smooth error curve, which is LOESS smoothed.
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
%        Field_Span - The fraction denominator of the maximum field over
%                     which to smooth (e.g., a value of 2 smooths over a
%                     field range equivalent to half the maximum field).
%                     This is converted into an equivalent number of points
%                     using the average field spacing. Default is 10 (i.e.,
%                     1/10 of the maximum field).
%
%        Type - Flag for the drift correction type
%               0 - No correction
%               1 - Automatic decision based on comapring low and high
%                   field noise (DEFAULT)
%               2 - Positive field correction [1]
%               3 - Upper branch only correction
%               4 - Symmetric averaging and tip-to-tip closure [2]
%               5 - Paramagnetic thermal dirft correction [3]
%
%        Chi_hat - An estimate for the high-field susceptibility. This is
%                  used for the exponential correction
%
% Output:
%         Corrected_moments - The drift corrected hysteresis moments [nData x 2]
%
%         Type - The drift correction type applied
%
%         Drift_Ratio - The ratio of median high-field drift to low-field
%                       drift. Used to decide between postive field (type 2)
%                       or upper branch correction (type 3) when automatic
%                       correction is selected.
%
%         Err_H - The hysteresis error/noise curve before drift correction.
%                 This corresponds to the upper branch fields (postive to
%                 negative sweep).
%
%         Smooth_Err_H - The smoothed hysteresis error/noise curve before drift correction.
%
%         Temp_Ratio - The ratio of specimen to ambient temperature for the
%                      paramagnetic thermal drift correction
%
% References:
%
% [1] Jackson, M., and P. Solheid (2010), On the quantitative analysis and
%     evaluation of magnetic hysteresis data, Geochem. Geophys. Geosyst.,
%     11, Q04Z15, doi:10.1029/2009GC002932.
%
% [2] von Dobeneck, T. (1996), A systematic analysis of natural magnetic mineral
%     assemblages based on modelling hysteresis loops with coercivity-related
%     hyperbolic basis functions, Geophys. J. Int., 124, 675-694,
%     doi:10.1111/j.1365-246X.1996.tb05632.x
%
% [3] Paterson, G. A., X. Zhao, M. Jackson, D. Heslop, Measuring, processing,
%     and analyzing hysteresis data, Geochemistry, Geophysics, Geosystems, 19,
%     doi: 10.1029/2018GC007620
%
% Last Modified 2019/05/07
%

%% Input processing and checking

if nargin < 2
    error('Hyst_Drift_Correction:Input', 'At least 2 input arguments are required.');
end

if nargin < 3
    Field_Span = 10;
    Type = 1; % Automatic correction
    Chi_hat = [];
end

if nargin < 4
    Type = 1; % Automatic correction
    Chi_hat =[];
end

if nargin < 5
    Chi_hat =[];
end

if isempty(Chi_hat) && Type ==5
    error('Hyst_Drift_Correction:Input', 'An estimate of paramagnetic susceptibility is need for the paramagnetic correction.');
end


%% Set some default parameters

% The percentage for defining "high fields"
HF_Lim = 0.75;

% The ratio threshold for applying positive field
Drift_Lim = 0.75;

% Limits on the number of points for smoothing the noise curve
nLim1 = 15; % The maximum limit - should never exceed this
nLim2 = 7; % The minimim limit - should never be below this

%% Get the error curve and smooth it

Err_H = Moment_Grid(:,1) + Moment_Grid(:,2);

% Define the number of points to smooth based on a field span
Max_Field = max(abs(Field_Grid(:,1)));
Field_Step = mean(abs(diff(Field_Grid(:,1))));
npts_smooth = ceil((Max_Field/Field_Span)/Field_Step);

if rem(npts_smooth,2) == 0
    % Even number so add 1
    npts_smooth = npts_smooth + 1;
end

% If too few points are used enforce a manual limit
npts_smooth = min([npts_smooth, nLim1]);
npts_smooth = max([npts_smooth, nLim2]);

Smooth_Err_H = Loess_Smooth(Err_H, npts_smooth);


%% Define correction type for automatic option

% Get the high field threshold
High_Field = Max_Field * HF_Lim;

% Get the high and low field absolute drifts
HFD = abs( Smooth_Err_H( abs(Field_Grid(:,1)) >= High_Field ) );
LFD = abs( Smooth_Err_H( abs(Field_Grid(:,1)) <  High_Field ) );

Drift_Ratio = median(HFD) / median(LFD);

% TEST CODE ONLY
% KS-test with the null hypothesis that the distribution are equal, against
% the alternative hypothesis that low-field drift ecdf is larger (plots
% above) than high-field (i.e., low-field drift tends to have lower values).
% if p <= 0.05 then high-field drift tends to be higher
% [h,p]=kstest2(LFD, HFD, 'tail', 'larger');

if Type == 1
    if Drift_Ratio >= Drift_Lim
        % Most drift is at high fields
        % Apply postive field correction
        Type = 2;
    else
        % Most drift is at low fields
        % Upper branch correction
        Type = 3;
    end
end


%% Do the correction

Corrected_Moments = Moment_Grid;

% Set to NaN for all other methods except paramagnetic
Temp_Ratio = NaN;

switch Type
    
    case 0
        % No correction, so do nothing
        
    case 2
        % Postive field correction
        
        % Correct the upper branch
        Pos_Fields = Field_Grid(Field_Grid(:,1)>=0,1);
        npts_PF = size(Pos_Fields,1);
        
        for jj = 1:npts_PF
            tmp_Field = Pos_Fields(jj);
            Err_corr = Smooth_Err_H(jj);
            Corrected_Moments(Field_Grid(:,1) == tmp_Field,1) = Moment_Grid(Field_Grid(:,1) == tmp_Field, 1) - Err_corr;
        end
        
        
        % Correct the lower branch
        Pos_Fields = Field_Grid(Field_Grid(:,2)>=0,2);
        npts_PF = size(Pos_Fields,1);
        
        for jj = 1:npts_PF
            tmp_Field = Pos_Fields(jj);
            Err_corr = Smooth_Err_H(Field_Grid(:,1) == -tmp_Field);
            Corrected_Moments(Field_Grid(:,2) == tmp_Field,2) = Moment_Grid(Field_Grid(:,2) == tmp_Field, 2) - Err_corr;
        end
        
    case 3
        % Upper branch correction
        
        Corrected_Moments(:,1) = Moment_Grid(:,1) - Smooth_Err_H;
        
        
    case 4
        % Symmetric averaging and tip-to-tip closure
        
        % Invert bottom curve
        Binv = [-Field_Grid(:,2), -Moment_Grid(:,2)];
        
        % Interpolate the inverted bottom curve to the top curve fields
        Binv_interp = Interpolate_To_Field(Binv(:,1), Binv(:,2), Field_Grid(:,1), 'linear');
        
        Loop_bar = mean([Moment_Grid(:,1), Binv_interp], 2);
        
        % Get the max/min values for the upper/lower branches
        Max_Upper_Moment = Moment_Grid(Field_Grid(:,1)==max(Field_Grid(:,1)), 1);
        Min_Upper_Moment = Moment_Grid(Field_Grid(:,1)==min(Field_Grid(:,1)), 1);
        Max_Lower_Moment = Moment_Grid(Field_Grid(:,2)==max(Field_Grid(:,2)), 2);
        Min_Lower_Moment = Moment_Grid(Field_Grid(:,2)==min(Field_Grid(:,2)), 2);
        
        dTip = ( Max_Upper_Moment - Max_Lower_Moment + Min_Upper_Moment - Min_Lower_Moment )/4;
        
        Corrected_Moments = [Loop_bar - dTip, dTip - Loop_bar];
        
    case 5
        % Paramagnetic correction based on a paramagnetic thermal drift model
        
        
        % The closure error and noise
        Mce = Err_H(1) - Err_H(end);
        Noise = [Field_Grid(:,1), Err_H, Smooth_Err_H];
        
        % The number of points for the whole loop
        npts = 2*size(Field_Grid(:,1),1);
        
        % Get some intial guess for the temperatures
        if Mce < 0
            % Upper branch is below the lower branch
            % Specimen temeprature decrease
            % Ambient is lower than the specimen
            % Tr > 1
            Tr_0 = 1.02;
            Tr_min = 1;
            Tr_max = 1.2;
        else
            % Upper branch is above the lower branch
            % Specimen temeprature decrease
            % Ambient is higher than the specimen
            % Tr < 1
            Tr_0 = 1/1.02;
            Tr_min = 1/1.2;
            Tr_max = 1;
        end
        
        % Do the optimization
        options = optimset('TolX', 1e-8, 'Tolfun', 1e-8, 'MaxIter', 1e4, 'MaxFunEvals', 1e4, 'Display', 'off');
        
        Opt_vals = fminsearchbnd(@(x) Paramagnetic_Drift_Model(Noise, Chi_hat, x(1), x(2)),...
            [Tr_0, 1/(0.2*npts)], [Tr_min, 0], [Tr_max, 1], options );
        
        
        % Output for checking
%         disp('  ')
%         disp('%%%%%%%%%%%%%%%%%%%%%%')
%         disp(['Xp = ', sprintf('%1.6e ', Chi_hat)]);
%         disp(['Tr = ', sprintf('%1.6f', Opt_vals(1))]);
%         disp(['c = ', sprintf('%1.6e', Opt_vals(2))]);
%         disp('%%%%%%%%%%%%%%%%%%%%%%')
%         disp('  ')
        
        [MisFit, Para_Err_H, Para_Drift] = Paramagnetic_Drift_Model(Noise, Chi_hat, Opt_vals(1), Opt_vals(2));
        
        Temp_Ratio = Opt_vals(1);
        
        % Get the deviation of the drifted paramagnetic signal from the
        % true paramagnetic signal
        Paramagnetic_Dev = [Para_Drift(:,1), Para_Drift(:,2) - Para_Drift(:,3)];
        
        % Subtract the paramagnetic deviation from the measured moment to
        % correct for paramagnetic drift
        Corrected_Moments = [Moment_Grid(:,1)+Paramagnetic_Dev(1:npts/2,2), Moment_Grid(:,2) + Paramagnetic_Dev(1+npts/2:end,2)];
        
        % Apply a second drift correction
        % 0 - None
        % 1 - Automatic
        % 2 - Positive field
        % 3 - Upper branch
        
        [Corrected_Moments, Type, Drift_Ratio, Err_H, Smooth_Err_H] = Hyst_Drift_Correction(Field_Grid, Corrected_Moments, Field_Span, 1, []);
        %          [Corrected_Moments, Type, Drift_Ratio, Err_H, Smooth_Err_H] = Hyst_Drift_Correction(Field_Grid, Corrected_Moments, Field_Span, 2, []);
%                  [Corrected_Moments, Type, Drift_Ratio, Err_H, Smooth_Err_H] = Hyst_Drift_Correction(Field_Grid, Corrected_Moments, Field_Span, 0, []);
        
        % Get the drift correction type flag - Add 5 for paramagnetic plus
        % other method
        Type = Type + 5;
        
    otherwise
        error('Hyst_Drift_Correction:Type', 'Unrecognized drift correction type. Please check.');
end


