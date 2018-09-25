function Save_HystLab_Stats(Names, Masses, Stats, Params, file, path)
%
% Function to save the hysteresis statistic and processing parameters to a
% tab delimited text file
%
% Input:
%        Names - Specimen names ([nx1] cell)
%        Masses - Vector of specimen masses
%        Stats - Matrix of hysteresis statistics
%        Params - Matrix of processing parameters
%        file - File to save the data to (string)
%        path - Directory path the save the file to (string)
%

%% Input check and definitions


if ~isequal(size(Names, 1), size(Masses,1), size(Stats,1), size(Params,1))
    error('Save_HystLab_Stats:Input', 'Inconsistent number of specimens in the input data.');
end


if sum(cellfun(@iscell, Names)) == length(Names)
%     keyboard
    Names = [Names{:}]';
end

% Define the stat names
Stat_Names = [{'H0'}; {'M0'}; {'Q'}; {'Qf'}; {'Ms'}; {'Mrs'}; {'Bc'}; {'Brh'}; {'Bih'};...
    {'Xhf'}; {'X0'}; {'Shape'}; {'C_err'}; {'Nfit'}; {'AS_alpha'}; {'AS_beta'};...
    {'AS_Lin_pVals1'}; {'AS_Lin_pVals2'}; {'AS_Lin_pVals3'}; {'AS_Lin_pVals4'}; {'AS_Lin_FVals1'}; {'AS_Lin_FVals2'}; {'AS_Lin_FVals3'}; {'AS_Lin_FVals4'};...
    {'S_star'}; {'Fit_F'}; {'Fit_p'}; {'Linear_F'}; {'Linear_p'};...
    {'AS_Model_pVals1'}; {'AS_Model_pVals2'}; {'AS_Model_pVals3'}; {'AS_Model_pVals4'}; {'AS_Model_FVals1'}; {'AS_Model_FVals2'}; {'AS_Model_FVals3'}; {'AS_Model_FVals4'};...
    {'Qrh'}; {'Qih'}; {'Bsat'}; {'Sat_pct'}; {'Noise_RMS'}; {'Fit_RMS'};...
    {'HAR_70'}; {'HAR_80'}; {'HAR_90'}; {'HAR_Selected'}; {'SNR_70'}; {'SNR_80'}; {'SNR_90'}; {'SNR_Selected'}];

if size(Stats, 2) ~= length(Stat_Names)
    error('Save_HystLab_Stats:Stats', 'Number of data statistics are incorrect.')
end

% Define the processing names
Param_Names = [{'Drift_Flag'}; {'Drift_Type'}; {'Drift_Ratio'}; {'Temp_Ratio'}; {'Saturation_Flag'}; {'Saturation_Field'}; {'Method'}; {'Pole_Flag'}; {'Offset_Flag'};...
    {'Trim_Flag'}; {'Trim_Field'}; {'Fit_Param_Flag'}; {'FixedBeta_Flag'}; {'FixedBeta_Val'}];

if size(Params, 2) ~= length(Param_Names)
    error('Save_HystLab_Stats:Parameters', 'Number of data processing parameters are incorrect.')
end


%% Defined the table for output

% Names to print
Print_Stat_Names = [{'Bc [mT]'}, {'Brh [mT]'}, {'Bih [mT]'}, {'Ms [Am^2]'}, {'Mrs [Am^2]'},  {'Xhf [m^3]'},  {'Shape'},...
    {'alpha'}, {'beta'}, {'Bsat [mT]'}, {'Sat %'},...
    {'B0 [mT]'}, {'M0 [Am^2]'}, {'Q'}, {'Qf'}, {'Qrh'}, {'Qih'}, {'Closure err. [Am^2]'}, {'Noise RMS [Am^2]'},...
    {'N fit'}, {'Fit F'}, {'Fit p'}, {'Fit RMS [Am^2]'}];

% Order to print
Print_Stat_Order = [{'Bc'}, {'Brh'}, {'Bih'}, {'Ms'}, {'Mrs'},  {'Xhf'},  {'Shape'},...
    {'AS_alpha'}, {'AS_beta'}, {'Bsat'}, {'Sat_pct'},...
    {'H0'}, {'M0'}, {'Q'}, {'Qf'}, {'Qrh'}, {'Qih'}, {'C_err'}, {'Noise_RMS'},...
    {'Nfit'}, {'Fit_F'}, {'Fit_p'},  {'Fit_RMS'}];

% Format to print
Print_Stat_Format = [{'%3.1f'}, {'%3.1f'}, {'%3.1f'}, {'%1.2e'}, {'%1.2e'},  {'%1.2e'},  {'%1.3f'},...
    {'%1.3f'}, {'%1.3f'}, {'%3.1f'}, {'%3.1f'},...
    {'%1.2f'}, {'%1.2e'}, {'%1.2f'}, {'%1.2f'}, {'%1.2f'}, {'%1.2f'}, {'%1.2e'}, {'%1.2e'},...
    {'%d'}, {'%5.1f'}, {'%3.2f'}, {'%1.2e'}];

nSpec = length(Names);
nPrint = length(Print_Stat_Order);
Print_Stats = cell(nSpec, nPrint);

for ii = 1:nPrint
    Print_Stats(:,ii) = num2cell(Stats(:,strcmp(Stat_Names, Print_Stat_Order{ii})));
end


% Names to print
Print_Param_Names = [{'Correct for Offset'}, {'Drift Method'}, {'Drift Correction Applied'}, {'Drift Ratio'}, {'Temperature Ratio'}, {'Trim Fields?'}, {'Trim Field >= [mT]'},...
    {'Saturation Method'}, {'Applied Saturation'}, {'Saturation Field [mT]'}, {'Used Fixed Beta'}, {'Fixed Beta'}, ...
    {'Parameters Estimated From Fit?'}];

% Order to print
Print_Param_Order = [{'Offset_Flag'}, {'Drift_Flag'}, {'Drift_Type'}, {'Drift_Ratio'}, {'Temp_Ratio'}, {'Trim_Flag'}, {'Trim_Field'},...
    {'Saturation_Flag'}, {'Method'}, {'Saturation_Field'}, {'FixedBeta_Flag'}, {'FixedBeta_Val'}, ...
    {'Fit_Param_Flag'}];


% Format to print
Print_Param_Format = [{'%s'}, {'%s'}, {'%s'}, {'%3.2f'}, {'%1.3f'}, {'%s'}, {'%4.1f'},...
    {'%s'}, {'%s'}, {'%4.1f'}, {'%s'}, {'%2.2f'}, ...
    {'%s'}];


nSpec = length(Names);
nPrint = length(Print_Param_Order);
Print_Params = cell(nSpec, nPrint);

% Param_Names = [{'Drift_Flag'}; {'Drift_Type'}; {'Drift_Ratio'}; {'Temp_Ratio'}; {'Saturation_Flag'}; {'Saturation_Field'}; {'Method'}; {'Pole_Flag'}; {'Offset_Flag'};...
%     {'Trim_Flag'}; {'Trim_Field'}; {'Fit_Param_Flag'}; {'FixedBeta_Flag'}; {'FixedBeta_Val'}];



for ii = 1:nPrint
    
    tmp_Params = num2cell(Params(:,strcmp(Param_Names, Print_Param_Order{ii})));
    
    switch ii
        
        case {1, 6, 11, 13}
            % Offset, trim field, fixed beta, params from fit
            % All yes/no
            tmp_Params(cellfun(@(x) isequal(x,1), tmp_Params)) = {'Yes'};
            tmp_Params(cellfun(@(x) isequal(x,0), tmp_Params)) = {'No'};
            
            
        case 2
            % Drift flag
            tmp_Params(cellfun(@(x) isequal(x,0), tmp_Params)) = {'None'};
            tmp_Params(cellfun(@(x) isequal(x,1), tmp_Params)) = {'Automatic'};
            tmp_Params(cellfun(@(x) isequal(x,2), tmp_Params)) = {'Positive Field'};
            tmp_Params(cellfun(@(x) isequal(x,3), tmp_Params)) = {'Upper Branch'};
            tmp_Params(cellfun(@(x) isequal(x,4), tmp_Params)) = {'Symmetric Averaging'};
            tmp_Params(cellfun(@(x) isequal(x,5), tmp_Params)) = {'Paramagnetic'};
            
        case 3
            % Applied drift correction
            
            tmp_Params(cellfun(@(x) isequal(x-1,0), tmp_Params)) = {'No Correction'};
            tmp_Params(cellfun(@(x) isequal(x-1,1), tmp_Params)) = {'Automatic'};
            tmp_Params(cellfun(@(x) isequal(x-1,2), tmp_Params)) = {'Positive Field'};
            tmp_Params(cellfun(@(x) isequal(x-1,3), tmp_Params)) = {'Upper Branch'};
            tmp_Params(cellfun(@(x) isequal(x-1,4), tmp_Params)) = {'Symmetric Averaging'};
            tmp_Params(cellfun(@(x) isequal(x-1,5), tmp_Params)) = {'Paramagnetic'};
            tmp_Params(cellfun(@(x) isequal(x-1,5), tmp_Params)) = {''};
            tmp_Params(cellfun(@(x) isequal(x-1,7), tmp_Params)) = {'Paramagnetic/Postive'};
            tmp_Params(cellfun(@(x) isequal(x-1,8), tmp_Params)) = {'Paramagnetic/Upper'};
            tmp_Params(cellfun(@(x) isequal(x-1,9), tmp_Params)) = {'Paramagnetic/Symmetric'};
            
        case 8
            % Saturastion flag
            tmp_Params(cellfun(@(x) isequal(x,0), tmp_Params)) = {'None'};
            tmp_Params(cellfun(@(x) isequal(x,1), tmp_Params)) = {'Automatic'};
            tmp_Params(cellfun(@(x) isequal(x,2), tmp_Params)) = {'Linear High-Field'};
            tmp_Params(cellfun(@(x) isequal(x,3), tmp_Params)) = {'Approach to Saturation'};
            
        case 9
            % Applied Saturation
            
            tmp_Params(cellfun(@(x) isequal(x,1), tmp_Params)) = {'No Correction'};
            tmp_Params(cellfun(@(x) isequal(x,2), tmp_Params)) = {'Linear High-Field'};
            tmp_Params(cellfun(@(x) isequal(x,3), tmp_Params)) = {'Approach to Saturation'};
            
    end
    
    Print_Params(:,ii) = tmp_Params;
end



% Create the header line and format
Final_Header = [{'Specimen ID'}, {'Mass [mg]'}, Print_Stat_Names, Print_Param_Names];
nHead = length(Final_Header);
Header_Format = repmat('%s\t', 1, nHead);
Header_Format(end-1:end) = []; % remove the trailing '\t'
Header_Format = strcat(Header_Format, '\n');

% Combine the formats
Final_Format = strjoin([{'%s'}, {'%3.3f'}, Print_Stat_Format, Print_Param_Format], '\\t');
Final_Format = strcat(Final_Format, '\n');

% Combine the final data
Final_Data = [Names, num2cell(Masses), Print_Stats, Print_Params]';



%% Save

FID = fopen(strcat(path, file), 'wt');

fprintf(FID, Header_Format, Final_Header{:});
fprintf(FID, Final_Format, Final_Data{:});

fclose(FID);

