function Save_HystLab_Session(handles, Type)
%
% Function to save HystLab sessions
%

%% Check inputs

if nargin < 2
    error('Save_HystLab_Session:Input', 'Two input arguments are required.');
end

if ~strcmpi(Type, 'Hys')
    error('Save_HystLab_Session:Input', 'The current session type is not supported.');
end

%% The main function

default_name = 'HystLab_Saved_Session.mat';


% The field names to save
fNames = [{'Data_Loaded'}, {'Nspec'}, {'spec_ind'}, {'Norm_Type'}, ...
    {'All_Names'}, {'All_Data'}, {'All_Masses'}, {'All_Data_Order'}, {'All_Processing_Parameters'}, {'All_Processed_Data'},...
    {'All_Uncorrected_Data'}, {'All_Noise_Data'}, {'All_Fitted_Data'}, {'All_Data_Parameters'}, {'All_Basis_Coeffs'},...
    {'Current_Raw_Loop'}, {'Current_Uncorrected_Data'}, {'Current_Processed_Data'}, {'Current_Fitted_Data'},...
    {'Current_Noise_Data'}, {'Current_Data_Parameters'}, {'Current_Spec_Mass'}...
    {'HystLab_Version'}];


% Get the output file
[file,path] = uiputfile(default_name,'Save session...');


if ~ischar(file) && file==0
    % User has cancelled
    % Do nothing and...
    return;
end


nFields = length(fNames);
Session_handles = struct();

for ii = 1:nFields
    Session_handles.(fNames{ii}) = handles.(fNames{ii});
end

save(strcat(path, file), 'Session_handles', 'Type');

