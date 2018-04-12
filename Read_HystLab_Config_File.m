function Defaults = Read_HystLab_Config_File(MyPath)
%
% Function to read in a HystLab configuration file, which contains
% various default settings
%

%%

% Get the basic setting names and types
All_Settings = [{'PlotColorFile'},...
    {'BiPlotFaceColor'}, {'BiPlotSymbol'}, {'BiPlotSymbolSize'},...
    {'DataFormat'}, {'DataFormatOther'}];


Setting_Types = [{'String'},...
    {'String'}, {'String'}, {'Number'},...
    {'String'}, {'String'}];



nSettings = length(All_Settings);

% Find the config file
if exist(strcat(MyPath, 'HystLabUserDefaults.cfg'), 'file') == 2
    % The file exists so load user defaults
    Config_File = strcat(MyPath, 'HystLabUserDefaults.cfg');
elseif exist(strcat(MyPath, 'HystLabDefaults.cfg'), 'file') == 2
    % Use the basic defaults
    Config_File = strcat(MyPath, 'HystLabDefaults.cfg');
else
    error('HystLab:Read_Config',...
        'No HystLab configuration file is present. Please check the HystLab path for the configuration file.');
end

% Open and read the file
FID=fopen(Config_File, 'r');
input=textscan(FID, '%s = %s\n');
fclose(FID);
nDefault = length(input{1});


if nDefault ~= nSettings
    % if the number of settings doesn't match try different line endings
    FID=fopen(Config_File, 'r');
    input=textscan(FID, '%s = %s\r\n');
    fclose(FID);
    nDefault = length(input{1});
end

if nDefault ~= nSettings
    % if the number of settings still doesn't match try the default config file
    warning('HystLab:ConfigFile', 'Configuration file inconsistent, attempting load default configuration file.')
    Config_File = strcat(MyPath, 'HystLabDefaults.cfg');
    FID=fopen(Config_File, 'r');
    input=textscan(FID, '%s = %s\n');
    fclose(FID);
    nDefault = length(input{1});
end


if nDefault ~= nSettings
    % if the number of settings STILL doesn't match - ERROR
    error('HystLab:ConfigFile', 'Default configuration file is corrupted or settings are missing.');
end

% The return structure
Defaults = struct();


for ii=1:nSettings
    
    % Check the input type
    if strcmpi(Setting_Types{ii}, 'String')
        if ~ischar(input{2}{ii})
            error('HystLab:ConfigFile', '%s should be a text input.', All_Settings{ii});
        end
    else
        if ~isnumeric(str2double(input{2}{ii}))
            error('HystLab:ConfigFile', '%s should be a number input.', All_Settings{ii});
        end
    end
    
    if strcmpi(All_Settings{ii}, 'PlotColorFile')
        
        try
            tmp_input = load(strcat(MyPath, input{2}{ii}) );
            Defaults.Hyst_Plot_Color = tmp_input.Hyst_Plot;
            Defaults.Noise_Plot_Color = tmp_input.Noise_Plot;
            Defaults.IR_Plot_Color = tmp_input.IR_Plot;
            Defaults.BiPlot_Color = tmp_input.BiPlot;
            
        catch
            error('HystLab:ConfigFile', 'Default color file (%s) not found.', input{2}{ii})
        end
        
    elseif strcmpi(Setting_Types{ii}, 'String')
        % It is a text input
        Defaults.(All_Settings{ii}) = input{2}{ii};
    else %number
        Defaults.(All_Settings{ii}) = str2double(input{1,2}{ii});
    end
end
