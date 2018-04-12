function [p, F] = Get_F_Test(M, Mhat, P, type)
%
% Function to perform F test for lack of fit following [1] and the F-test
% for model variance comaprison
%
%
% Input:
%       M - the measurement data (size depends on test type)
%       Mhat - the model estimate (size depends on test type)
%       P - the number of model parameters
%       type - flag to denote if the call is for:
%              1) Lack of fit test test (e.g., whole loop linear model or whole loop basis fitting)
%              2) F-test model comparison (e.g., linear or nonlinear high field slope correction)
%
% Output:
%       p - the p-value of the test
%       F - the F-value of the test
%
%
% References:
% [1] Jackson, M., and P. Solheid (2010), On the quantitative analysis and
%     evaluation of magnetic hysteresis data, Geochem. Geophys. Geosyst.,
%     11, Q04Z15, doi:10.1029/2009GC002932.
%

%% Check the inputs

if nargin < 4
    error('Get_F_Test:Input', 'Four inputs are required.')
end

if ~isequal(size(M), size(Mhat))
    error('Get_F_Test:Input', 'M and Mhat should be equal in size.')
end


%% Begin

switch type
    
    case 1 % Lack of Fit F-test
        
        % Here we wish to test if we can reject a particular model fit to
        % the data based on comparing the misfit of the model to the noise
        % estimate from repeated measurements (i.e., a lack of fit)
        % Specifically, we are testing wheter deviations from the model are
        % statistically significant or indistinguishable from noise
        
        % If the calculated p value is less than, say 0.05, then we can
        % reject the null hypothesis that there is no lack of fit.
        % That is, the misfit of the model cannot be explained by the noise
        % of the data (i.e., the model not does explains the data)
        
        % M and Mhat should be [npts x 1]
        % M is the observed data
        % Mhat is the model fit to the data
        
        % P is scalar for the number of model parameters
        % P = 2 for a linear model
        % P = 4 for the approach to saturation model
        % P = number of basis funtions for loop fitting
        
        if length(P) ~=1
            warning('Get_F_Test:P', 'For type 1 tests P should be a scalar. Taking first value.');
            P = P(1);
        end
        
        npts = length(M);
        
        % Check the minimum nuber of poits for the F-test
        nMin = 2*(1+P);
        
        if npts < nMin
            F = NaN;
            p = NaN;
            return;
        end
        
        
        if mod(npts, 2) ~= 0
            error('Get_F_Test:npts', 'The number of data points should be even. Should not be here');
        end
        
        % Get the degrees of freedom
        df1 = (npts/2 - P);
        df2 = npts/2;
        
        if df1 <= 0 || df2 <= 0
            error('Get_F_Test:DoF', 'Degrees of freedom must be greater than zero.');
        end
        
        
        SSD = sum(sum( (Mhat - M).^2 ));
        SSPE = sum( 0.5 * (M(1:npts/2) - M(npts/2+1:end)).^2);
        
        SSLF = SSD - SSPE;
        
        MSLF = SSLF / df1;
        MSPE = SSPE / df2;
        
        F = MSLF/MSPE;
        
        p = 1 - Get_fcdf(F, df1, df2);
        
        
    case 2 % Slope correction comparison
        
        % Here we wish to test if a simple linear high-field fit (model 1)
        % yields a statistically signifcant better fit than a non-linear
        % approach to saturation fit (model 2)
        
        % If the calculated p value is less than, say 0.05, then we can
        % reject the null hypothesis that model 1 fits the data best, in
        % favour of model 2 (specimen has not reached magnetic saturation)
        
        % M and Mhat should be [npts x 2]
        % M is lower branch of the hysteresis loop sweeping to high
        % positive fields, replicated in each column
        % Mhat(:,1) is the linear model fit, Mhat(:,2) is the appraoch to
        % saturation fit
        
        % P is [1 x 2] vector of model parameters [model 1, model 2] = [2, 4]
        
        % Do some input checks
        if size(M,2) ~= 2 || length(P) ~=2
            error('Get_F_Test:SlopeCorrection', 'M or P have the wrong size.');
        end
        
        npts = size(M,1);
        SSD = sum( (M - Mhat).^2 );
        
        df = npts - P;
        
        Ft = (SSD(1) - SSD(2)) / (df(1)-df(2));
        Fb = SSD(2) / (df(2));
        
        F = Ft/Fb;
        p = 1 - Get_fcdf(F, df(1)-df(2), df(2));
        
    otherwise
        
        error('Get_F_Test:Type', 'Undefined fit type %d', type);
        
end


