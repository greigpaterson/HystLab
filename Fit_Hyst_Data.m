function [Fitted_Data, Basis_Coeffs, P, Fit_F, Fit_p, Fit_RMS] = Fit_Hyst_Data(Moments, Fields)
%
% Function to fit Mih and Mrh curves using hyperbolic and logistic
% functions following the procedures of [1-3]
%
% Input:
%       Moments - The hysteresis moments on a regular field grid
%                 [nData x 2]. Column 1 is the upper branch sweeping
%                 positive to negative fields. Column 2 is the lower branch
%                 sweeping negative to postive fields
%
%       Fields - The gridded fields (only one sweep, positive to negative)
%
% Output:
%        Fitted_Data - Fitted hysteresis loop using the gridded data [nField x 5]
%                     [Fields, Top_curve, Bottom_Curve, Mih, Mrh]
%        Basis_Coeff - The basis function coefficients used for fitting
%        P - The number of basis functions fitted
%        Fit_F - The Lack of Fit F-test F-value
%        Fit_p - The Lack of Fit F-test p-value
%
% References:
%
% [1] von Dobeneck, T. (1996), A systematic analysis of natural magnetic mineral
%     assemblages based on modelling hysteresis loops with coercivity-related
%     hyperbolic basis functions, Geophys. J. Int., 124, 675?694,
%     doi:10.1111/j.1365-246X.1996.tb05632.x
%
% [2] Jackson, M., and P. Solheid (2010), On the quantitative analysis and
%     evaluation of magnetic hysteresis data, Geochem. Geophys. Geosyst.,
%     11, Q04Z15, doi:10.1029/2009GC002932.
%
% [3] Paterson, G. A., X. Zhao, M. Jackson, D. Heslop, Measuring, processing,
%     and analyzing hysteresis data, Geochemistry, Geophysics, Geosystems, 19,
%     doi: 10.1029/2018GC007620
%
%
% Last Modified 2021/05/27
%

%% Some input processing and defaults

if nargin < 2
    error('Fit_Hyst_Data:Input', 'At least 2 input arguments are required.')
end

% Reinterpolate the moments to ensure there is a zero field step
% This ensures that better fitting as zero field IH curves are zero, and RH
% curves are Mrs.
npts = length(Fields);
if ~any(Fields == 0)
    
    Old_Fields = Fields;
    Old_Moments = Moments;
    
    % Redefine the fields using 1 less point
    % Each branch should now have an odd number of points and a zero field
    Fields = linspace(Old_Fields(1), Old_Fields(end), npts-1)';
    npts = npts-1;
    Fields(ceil(npts/2)) = 0; % Ensure the mid field is zero
    
    % Reinterpolate the moments
    Moments = [interp1(Old_Fields, Old_Moments(:,1), Fields), interp1(-Old_Fields, Old_Moments(:,2), -Fields)];
    Moments(:,1) = interp1(Old_Fields, Old_Moments(:,1), Fields);
    Moments(:,2) = interp1(-Old_Fields, Old_Moments(:,2), -Fields);
    
end

% Number of basis functions for each type (total is 4*nBasis + 4)
% Impose a minimum based on the number of data points
% This ensures we have no problems when doing the F-test for lack-of-fit

nBasis_Limit = (npts/2-4)/4; % The limit on the number of bases as imposed by the F- lack-of-fit test
nBasis = min(10, floor(nBasis_Limit));


% The tolerance on zero abundances (100*Atol)%
Atol = 1e-4;

%% Process the input

% Get the curves
Mih = 0.5.*(Moments(:,1) + flipud(Moments(:,2)));
Mrh = 0.5.*(Moments(:,1) - flipud(Moments(:,2)));

% Flip and average the data
Pos_Fields = sort(Fields(Fields>=0));
IH = [abs(Mih(Fields==0)); mean([-Mih(Fields<0), flipud(Mih(Fields>0))], 2) ]';
RH = [abs(Mrh(Fields==0)); mean([ Mrh(Fields<0), flipud( Mrh(Fields>0))], 2)]';

% Get the Qih and Qrh values
% These are used for deciding on the fit
tmp_Mih1 = Mih(Fields > 0);
tmp_Mih2 = flipud(-Mih(Fields < 0));
R2 = GetR2(tmp_Mih1, tmp_Mih2);
Qih = log10( 1 / sqrt((1-R2)) );

tmp_Mrh1 = Mrh(Fields > 0);
tmp_Mrh2 = flipud(Mrh(Fields < 0));
R2 = GetR2(tmp_Mrh1, tmp_Mrh2);
Qrh = log10( 1 / sqrt((1-R2)) );



%% Define some paramters for fitting, normalization etc

% The limit (L) of the sigmoidal functions near zero field
% Used to define the b parameter
% Generally, if the specimen is noisy a larger value is better as it
% supress the manifestion of "steps" in the fitted loop
% This essentially reduces the number of sigmoid funtions being fitted
% This is empirically derived

if Qih < 1.2
    Lih = .2 / 1e-1;
else
    Lih = 1 / 1e-4;
end

if Qrh < 1
    Lrh = 1 / 1e-1;
else
    Lrh = 1 / 1e-4;
end


% Field range
MaxField = max(Pos_Fields);
MinField = min(Pos_Fields(Pos_Fields>0)/1.2); % mT  - the minimum median field for the basis functions

% Fields for interpolating mangetizations
% Setting first field to a small non-zero value avoids some errors
tmp_Fields = Pos_Fields;
tmp_Fields(1) = MinField;



%% Do the RH curve

Unique_RH = unique(RH);

if length(Unique_RH) == 1
    % Deal with the constant Mrh case
    
    if Unique_RH ==0
        % Create a matrix of the basis coefficients
        % [Flag, Abund, Coeff1, Coeff2]
        Basis_Coeffs_rh = [4, 1, NaN, NaN];
        
        rh_Scale = 0;
        
    else
        error('Fit_Hyst_Data:Mrh', 'A non-zero constant Mrh curve is not supported by the fitting routine.')
    end
    
else
    % Get smoothed monotonic curves to find fields of evenly spaced moments
    
    % Normalize the curve, get an isotonic fit, and renormalize
    RH_norm = RH./max(RH);
    RH_norm = IsoReg(fliplr(RH_norm)');
    RH_norm = fliplr(RH_norm'./max(RH_norm));
    
    % Find the unique values for the rh curve
    [vals, idx] = unique(RH_norm);
    idx(vals==max(vals)) = 1;
    idx(vals==min(vals)) = length(RH_norm);
    
    % Interpret between the unique values to define a full smoothed rh curve
    if length(vals) == 1
        % Really non-monotonic behavior
        % This indicates very noise data
        % Define decreasing linear trend
        RH_norm = linspace(1, 0, length(RH_norm));
    else
        RH_norm = interp1(Pos_Fields(idx), vals, Pos_Fields, 'linear');
    end
    
    
    % Define the basis function coefficients
    
    % Fields of evenly spaced moments
    if RH_norm(1) ~= 1
        RH_norm(1) = interp1(Pos_Fields, RH_norm, tmp_Fields(1), 'pchip');
    end
    
    pcts = linspace(RH_norm(1), RH_norm(end), nBasis);
    
    % Logistic functions
    
    if length(RH_norm) ~= length(unique(RH_norm))
        % Catch noisy data where some values in the RH curve are non-unique
        [dum_var, ia]=unique(RH_norm);
        interp_array = sortrows([tmp_Fields(ia), RH_norm(ia)],1);
        
        Brh_m =  interp1(interp_array(:,2), interp_array(:,1), pcts, 'pchip');
        
    else
    
        Brh_m =  interp1(RH_norm, tmp_Fields, pcts, 'pchip');
    
    end
    
    b_rh = log( Lrh - 1 ) ./ Brh_m;
    
    
    % the sech coeffs
    s = log(2+sqrt(3))./Brh_m;
    
    
    % The basis function for the Mrh curves
    sig2_Basis = NaN(nBasis, length(Pos_Fields));
    sech_Basis = NaN(nBasis, length(Pos_Fields));    
    
    for ii = 1:nBasis
        sig2_Basis(ii,:) = 1 + -1./ (1 + exp(-b_rh(ii).*(Pos_Fields - Brh_m(ii)) ) );  % A = 1, C = 0
        sech_Basis(ii,:) = (sech(s(ii).*Pos_Fields)-sech(s(ii).*MaxField))./(1-sech(s(ii).*MaxField));
    end   
    
    
    % Combine the basis functions
    Basis_rh = [ sig2_Basis; sech_Basis; linspace(1, 0, length(Pos_Fields)); linspace(0, 1, length(Pos_Fields))];
    
    
    % Debug
    %
    % % Uncomment this section to examine the basis functions
    %
    % figure()
    % plot(Pos_Fields, Basis_rh)
    % keyboard
    
    
    % Find the RH basis contributions
    
    % Normalize to sum to one to get the abundances
    RH2 = RH./sum(RH);
    rh_Basis2 = bsxfun(@rdivide, Basis_rh, sum(Basis_rh,2));
    Abunds_rh = SUNSAL(RH2, rh_Basis2);
    Abunds_rh = Abunds_rh./sum(Abunds_rh);
    
    
    % Create a matrix of the basis coefficients
    Basis_Coeffs_rh = NaN(size(Basis_rh,1), 4);
    Basis_Coeffs_rh(:,1) = [ones(nBasis,1); 2.*ones(nBasis,1); 3; -3]; % Indices for the type of basis
    Basis_Coeffs_rh(:,2) = Abunds_rh';
    Basis_Coeffs_rh(:,3:4) = [ [b_rh', Brh_m']; [s', MaxField.*ones(nBasis,1)]; [1,0]; [1,0] ];
    
    
    % Loop through to remove low abundance basis functions
    nLow = sum(Abunds_rh < Atol);
    iCount_rh = 0;
    
    while nLow > 0
        
        % Find "zero" inds
        Zero_inds = find(Abunds_rh < Atol);
        
        % Check that we are not rejecting all basis functions
        % This should not be the case due to sum-to-one constraint
        if length(Zero_inds) >= size(rh_Basis2, 1)
            warning('Fit_Hyst_Data:Zero_Abunds', 'At least 1 basis should be fitted to the remanent hysteresis curve');
            break;
        end
        
                
        % Recalculate the abundances
        rh_Basis2(Zero_inds,:) = [];
        Abunds_rh = SUNSAL(RH2, rh_Basis2);
        Abunds_rh = Abunds_rh./sum(Abunds_rh); % Renormalize
        
        % Remove the zero abundace data and update with new abunds
        Basis_Coeffs_rh(Zero_inds,:) = [];
        Basis_Coeffs_rh(:,2) = Abunds_rh';
        
        nLow = sum(Abunds_rh < Atol);
        
        iCount_rh = iCount_rh +1;
        
    end
    
    
    % Make the model
    Model = Abunds_rh*rh_Basis2;
    
    % Mrs estimates can be noisy, so optimize scaling to the whole Mrh curve
    if any(Fields == 0)
        % Stop the zero field step being doubled
        Model_full = fliplr([fliplr(Model), Model(2:end)]);
    else
        Model_full = fliplr([fliplr(Model), Model(1:end)]);
    end
    
    scale_func = @(x) norm(Model_full.*x - Mrh');
    
    % Get the mean of the low-field branches
    Data_LF = mean(abs(Mrh(abs(Fields)<=0.1*MaxField)));
    Model_LF = mean(abs(Model_full(abs(Fields)<=0.1*MaxField)));
    
    X0 = Data_LF/Model_LF;
    
    rh_Scale = fminsearchbnd(scale_func, X0, 0, max(Mrh)./max(Model_full));
    
end


%% Do the IH curve

% Get smoothed monotonic curves to find fields of evenly spaced moments

% Normalize the curve, get an isotonic fit
IH_norm = IH./max(IH);
IH_norm = IsoReg(IH_norm');


if length(unique(IH_norm)) < 2
    % Non-monotonic data beyond simple noise
    % This usually happens when the specimen is strongly diamagnetic and
    % no slope correction is applied
    
    % Renormalize
    IH_norm = IH./max(IH);
    
    % Get the cumualtive sum of the absolute differences
    dIH = diff(IH_norm);
    IH_norm = [IH_norm(1), cumsum(abs(dIH))];
    IH_norm = IH_norm./max(IH_norm);
    IH_norm = IsoReg(IH_norm');
end

% Renomalize the curve
IH_norm = IH_norm'./max(IH_norm);

% Find the unique values for the ih curve
[vals, idx] = unique(IH_norm);
idx(vals==max(vals)) = length(IH_norm);
idx(vals==min(vals)) = 1;

% Interpret between the unique values to define a full smoothed ih curve
IH_norm = interp1(Pos_Fields(idx), vals, Pos_Fields, 'linear');


% Define the basis function coefficients

% Fields of evenly spaced moments
IH_norm(1) = interp1(Pos_Fields, IH_norm, tmp_Fields(1), 'pchip');
pcts = linspace(IH_norm(1), IH_norm(end), nBasis);

% Logistic functions
% catch and remove duplicates of normalized moments
if length(unique(IH_norm)) < length(IH_norm)
    tmp_var = unique(IH_norm);
    tmp_n = length(tmp_var);
    tmp_IH = NaN(tmp_n, 2);
    
    for ii = 1:tmp_n
        tmp_IH(ii,:) = [tmp_var(ii), mean(tmp_Fields(IH_norm==tmp_var(ii)))];
    end
    
    IH_norm = tmp_IH(:,1);
    tmp_Fields = tmp_IH(:,2);
    
    clear vars tmp_var tmp_n;
end

Bih_m =  interp1(IH_norm, tmp_Fields, pcts, 'pchip');
b_ih = log( Lih - 1 ) ./ (Bih_m);

% the tanh coeffs
t = log(3)./(2*Bih_m);


% The basis function for the Mih curves
sig1_Basis = NaN(nBasis, length(Pos_Fields));
tanh_Basis = NaN(nBasis, length(Pos_Fields));

for ii = 1:nBasis
    sig1_Basis(ii,:) = 1 ./ (1 + exp(-b_ih(ii).*(Pos_Fields - Bih_m(ii)) ) ); % A = 0, C = 1
    tanh_Basis(ii,:) = tanh(t(ii).*Pos_Fields) ./ tanh(t(ii).*MaxField);
end


% Combine the basis functions
Basis_ih = [ sig1_Basis; tanh_Basis; linspace(0, 1, length(Pos_Fields)); linspace(1, 0, length(Pos_Fields)) ];


% Debug
%
% % Uncomment this section to examine the basis functions
%
% figure()
% plot(Pos_Fields, Basis_ih)
% keyboard


    % Find the IH basis contributions

% Normalize to sum to one to get the abundances
IH2 = IH./sum(IH);
ih_Basis2 = bsxfun(@rdivide, Basis_ih, sum(Basis_ih,2));
Abunds_ih = SUNSAL(IH2, ih_Basis2);
Abunds_ih = Abunds_ih./sum(Abunds_ih); % Renormalize


% Create a matrix of the basis coefficients
Basis_Coeffs_ih = NaN(size(Basis_ih,1), 4);
Basis_Coeffs_ih(:,1) = [ones(nBasis,1); 2.*ones(nBasis,1); 3; -3]; % Indices for the type of basis
Basis_Coeffs_ih(:,2) = Abunds_ih';
Basis_Coeffs_ih(:,3:4) = [ [b_ih', Bih_m']; [t', MaxField.*ones(nBasis,1)]; [1,0]; [0,1] ];


% Loop through to remove low abundance basis functions
nLow = sum(Abunds_ih < Atol);
iCount_ih = 0;

while nLow > 0
    
    % Find "zero" inds
    Zero_inds = find(Abunds_ih < Atol);
    
    % Check that we are not rejecting all basis functions
    % This should not be the case due to sum-to-one constraint
    if length(Zero_inds) >= size(ih_Basis2, 1)
        warning('Fit_Hyst_Data:Zero_Abunds', 'At least 1 basis should be fitted to the induced hysteresis curve');
        break;
    end
    
    
    % Recalculate the abundances
    ih_Basis2(Zero_inds,:) = [];
    Abunds_ih = SUNSAL(IH2, ih_Basis2);
    Abunds_ih = Abunds_ih./sum(Abunds_ih); % Renormalize
    
    % Remove the zero abundace data and update with new abunds
    Basis_Coeffs_ih(Zero_inds,:) = [];
    Basis_Coeffs_ih(:,2) = Abunds_ih';
    
    nLow = sum(Abunds_ih < Atol);
    
    iCount_ih = iCount_ih +1;
    
end

% Make the model
Model = Abunds_ih*ih_Basis2;

% Ms estimates can be noisy, so optimize scaling to the whole Mih curve
if any(Fields == 0)
    % Stop the zero field step being doubled
    Model_full = fliplr([fliplr(-Model), Model(2:end)]);
else
    Model_full = fliplr([fliplr(-Model), Model(1:end)]);
end

scale_func = @(x) norm(Model_full.*x - Mih');

% Get the mean of the high-field branches
Data_HF = mean(abs(Mih(abs(Fields)>=0.7*MaxField)));
Model_HF = mean(abs(Model_full(abs(Fields)>=0.7*MaxField)));

X0 = Data_HF/Model_HF;
ih_Scale = fminsearchbnd(scale_func, X0, X0/3, 3*X0);


%% Get the model data at the hysteresis loop fields

% The total number of basis functions used
P = size(Basis_Coeffs_ih,1) + size(Basis_Coeffs_rh,1);

Basis_Coeffs = [{Basis_Coeffs_ih}, {Basis_Coeffs_rh}];

% Get the fit and the measured moment data

% This direct model curve works best for densly measured data sets
% if exist('Old_Fields', 'var')
%     Fitted_Data = Get_Fitted_Data(Old_Fields, Basis_Coeffs, 'Hys', [ih_Scale, rh_Scale]);
%     M = [Old_Moments(:,1); -Old_Moments(:,2)];
% else
%     Fitted_Data = Get_Fitted_Data(Fields, Basis_Coeffs, 'Hys', [ih_Scale, rh_Scale]);
%     M = [Moments(:,1); -Moments(:,2)];
% end

% This interpolation approach works better for sparse data
Fitted_Data = Get_Fitted_Data(Fields, Basis_Coeffs, 'Hys', [ih_Scale, rh_Scale]);
M = [Moments(:,1); -Moments(:,2)];

if exist('Old_Fields', 'var')
    % Reinterpolate to the input fields if needed
    tmp_Fit = NaN(length(Old_Fields), 6);
    tmp_Fit(:,1:2) = [Old_Fields, -Old_Fields];
    tmp_Fit(:,[3,5,6]) = interp1(Fitted_Data(:,1), Fitted_Data(:,[3,5,6]), tmp_Fit(:,1));
    tmp_Fit(:,4) = interp1(Fitted_Data(:,2), Fitted_Data(:,4), tmp_Fit(:,2));
    
    Fitted_Data = tmp_Fit;
    M = [Old_Moments(:,1); -Old_Moments(:,2)];
end


% Get the fit quality to the whole loop

% The modeled moments
Mhat = [Fitted_Data(:,3); -Fitted_Data(:,4)];

% Fit RMS
Fit_RMS = sqrt( mean( (M(:) - Mhat(:)).^2 ) );

[Fit_p, Fit_F] = Get_F_Test(M, Mhat, P, 1);

