function varargout = HystLab(varargin)
% HYSTLAB MATLAB code for HystLab.fig
%      HYSTLAB, by itself, creates a new HYSTLAB or raises the existing
%      singleton*.
%
%      H = HYSTLAB returns the handle to a new HYSTLAB or the handle to
%      the existing singleton*.
%
%      HYSTLAB('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in HYSTLAB.M with the given input arguments.
%
%      HYSTLAB('Property','Value',...) creates a new HYSTLAB or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before HystLab_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to HystLab_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help HystLab

% Last Modified 2021/06/18
%

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @HystLab_OpeningFcn, ...
    'gui_OutputFcn',  @HystLab_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before HystLab is made visible.
function HystLab_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to HystLab (see VARARGIN)

% Choose default command line output for Select_EndMembers
handles.output = hObject;

% Get the MATLAB and HystLab versions
Ver = ver('MATLAB');
handles.Version = str2double(Ver.Version);

% TODO - Update version number and date to read from a file
% Stops having to update this main function for fixes elsewhere
handles.HystLab_Version = '1.1.0';
handles.HystLab_Date = 'June, 2021';

% Get the screen dpi
set(0, 'Units', 'Pixels');
Sp = get(0, 'ScreenSize');
handles.Screen_Res = Sp(3:4);

try
    handles.Screen_DPI = get(0, 'ScreenPixelsPerInch');
catch %#ok<CTCH>
    set(0, 'Units', 'Inches');
    Si = get(0, 'ScreenSize');
    set(0, 'Units', 'Pixels');
    dpi = Sp./Si;
    handles.Screen_DPI = mean(dpi(3:4));
end

% Set the postion for if we are running standalone
parentPosition = Sp;

% TODO - add fucntionality for changing default processing parameters

% Get default settings
% Get the current path of the main m-file
S = mfilename('fullpath');
name_len = length(mfilename());
MyPath = S(1:end-name_len);

Defaults = Read_HystLab_Config_File(MyPath);

handles.Default_Hyst_Color = Defaults.Hyst_Plot_Color;
handles.Default_IR_Color = Defaults.IR_Plot_Color;
handles.Default_Noise_Color = Defaults.Noise_Plot_Color;


% Get the OS and line endings for outputting files
if ispc
    handles.OS = 'Win';
    handles.Line_End = '\r\n';
elseif ismac
    handles.OS = 'Mac';
    handles.Line_End = '\n';
elseif isunix
    handles.OS = 'Unix';
    handles.Line_End = '\n';
else
    handles.OS = 'Unknown';
    handles.Line_End = '\n';
end


% Check for input and if we are call as a stand alone
if isempty(varargin)
    % We have called it as a standalone GUI
    
        % Position to be relative to the screen:
    currentPosition = get(hObject, 'Position');
    newX = parentPosition(1) + (parentPosition(3)/2 - currentPosition(3)/2);
    newY = parentPosition(2) + (parentPosition(4)/2 - currentPosition(4)/2);
    newW = currentPosition(3);
    newH = currentPosition(4);
    set(hObject, 'Position', [newX, newY, newW, newH]);
    
    % Get the defaults
    
    % Set data loaded status
    handles.Data_Loaded = 0;
    handles.Cancel_Processing = 0;
    
    handles.Stand_Alone = 1;
    guidata(hObject, handles);
    
else
    
    
    % dontOpen = false;
    DataTransfer = find(strcmp(varargin, 'DataTransfer'));
    if (isempty(DataTransfer)) ...
            || (length(varargin) <= DataTransfer)
        clear var DataTransfer;
    else
        DataTransfer = varargin{DataTransfer+1};
    end
    
    
    % Get the main window
    MainWindow = findall(0,'type','figure', 'name', 'MEMA');
    if isempty(MainWindow)
        % Can't find the main window - shouldn't be here
        error('HystLab:NoMainWindow', 'A standalone call to HystLab has been made with input arguments. This is not supported.')
    else
        parentPosition = get(MainWindow, 'Position');
        handles.Version = getappdata(MainWindow, 'Version');
        DPI = getappdata(MainWindow, 'DPI');
    end
    
    % Position to be relative to parent:
    currentPosition = get(hObject, 'Position');
    newX = parentPosition(1) + (parentPosition(3)/2 - currentPosition(3)/2);
    newY = parentPosition(2) + (parentPosition(4)/2 - currentPosition(4)/2);
    newW = currentPosition(3);
    newH = currentPosition(4);
    set(hObject, 'Position', [newX, newY, newW, newH]);
        
    % Set data loaded and cancel flag status
    handles.Data_Loaded = 0;
    handles.Cancel_Processing = 0;
    
    % Check if data has been sent
    if exist('DataTransfer', 'var')
        handles.All_Names = DataTransfer.Names;
        handles.All_Data = DataTransfer.Data;
        handles.All_Masses = DataTransfer.Masses;
        handles.All_Processing_Parameters = DataTransfer.All_Processing_Parameters;
        handles.Nspec = length(handles.All_Names);
        
        % Do the inital default data processing
        handles = Set_Initial_Data(handles, 1);
        
        handles.Data_Loaded = 1;
        handles.Stand_Alone = 0;
        
                guidata(hObject, handles);
        uiwait(handles.HystLab_Fig);
    else        
        error('HystLab:DataTransfer', 'HystLab was called with incorrect input arguements');
    end
    
end

% Set key data to appdata
setappdata(handles.HystLab_Fig, 'LineEnd', handles.Line_End);
setappdata(handles.HystLab_Fig, 'OS', handles.OS);
% setappdata(handles.HystLab_Fig, 'Resolution', handles.Screen_Res);
% setappdata(handles.HystLab_Fig, 'DPI', handles.Screen_DPI);
setappdata(handles.HystLab_Fig, 'Defaults', Defaults );
% setappdata(handles.HystLab_Fig, 'Version', handles.Version);


% --- Outputs from this function are returned to the command line.
function varargout = HystLab_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
% varargout{1} = handles.output;

% delete(hObject)


% --- Executes on button press in PB_OK.
function PB_OK_Callback(hObject, eventdata, handles)
% hObject    handle to PB_OK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in PB_Cancel.
function PB_Cancel_Callback(hObject, eventdata, handles)
% hObject    handle to PB_Cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Call the close function
HystLab_Fig_CloseRequestFcn(handles.HystLab_Fig, [], handles)


%--- Executes when user attempts to close HystLab_Fig.
function HystLab_Fig_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to HystLab_Fig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure

if handles.Data_Loaded ~= 0
    
    choice = questdlg('Do you want to save before quitting?', 'Save Session?', 'Yes', 'No', 'Cancel', 'Cancel');
    
    if strcmpi(choice, 'Yes')
        MB_Save_Session_Callback([], [], handles);
    end
    
    if strcmpi(choice, 'Cancel')
        return;
    else
        % Check if we need to resume
        if handles.Stand_Alone == 0
            uiresume(handles.HystLab_Fig);
        end
        delete(hObject);
    end
    
    
else
    
    % No data loaded, just ask to quit or not
    choice = questdlg('Do you want to quit?', 'Quit HystLab?', 'Yes', 'No', 'No');
    
    if strcmpi(choice, 'Yes')
        
        
        % Check if we need to resume
        if handles.Stand_Alone == 0
            uiresume(handles.HystLab_Fig);
        end
        delete(hObject);
        
    else
        return;
    end
    
end


function handles = Set_Initial_Data(handles, Type)
% Get the data from previous window

% Error flag to catch processing issues
Error_Flag = 0;
Error_MSG = '';

% Do a waitbar
Cancel_Flag = 0;
h = waitbar(0,'Initializing....', 'Name', 'Processing hysteresis data...',...
    'CreateCancelBtn', 'setappdata(gcbf,''Cancelled'',1)');
setappdata(h,'Cancelled',0)


if Type == 0
    % Do the inital default processing
    
    nSpec = handles.Nspec;
    
    % pre-allocate
    Processed_Data = cell(nSpec,1);
    Uncorrected_Data = cell(nSpec,1);
    Noise_Data = cell(nSpec,1);
    Fitted_Data = cell(nSpec,1);
    Data_Parameters = NaN(nSpec, 51);
    Processing_Parameters = NaN(nSpec, 14);
    Basis_Coeffs = cell(nSpec,2);
    
    Bad_Spec = []; % Cell array to catch list of bad specimen name
    
    for ii = 1:nSpec
        
        % Check for Cancel button press
        if Cancel_Flag == 1 || getappdata(h,'Cancelled')
            Cancel_Flag = 1;
            break;
        end
        
        % Update the waitbar and continue
        waitbar((ii-1)/(nSpec), h, strcat(sprintf(' %d', ii-1), ' specimens processed'))
        
        try
        % No processing
%                 [Processed_Data(ii), Uncorrected_Data(ii), Noise_Data(ii), Fitted_Data(ii), Data_Parameters(ii,:), Processing_Parameters(ii,:), Basis_Coeffs(ii,:)]...
%                     = Process_Hyst_Data(handles.All_Data(ii), handles.All_Data_Order(ii), 'Offset', 0, 'Drift', 0, 'Saturation', 0,...
%                     'SaturationField', [], 'FixedBeta_Flag', 0, 'FixedBeta_Val', -1.5,...
%                     'PoleSaturation', 0, 'PoleData', [], 'Trim', 0, 'TrimField', [], 'FitEst',0);
        
        % Default processing
                [Processed_Data(ii), Uncorrected_Data(ii), Noise_Data(ii), Fitted_Data(ii), Data_Parameters(ii,:), Processing_Parameters(ii,:), Basis_Coeffs(ii,:)]...
                    = Process_Hyst_Data(handles.All_Data(ii), handles.All_Data_Order(ii), 'Offset', 1, 'Drift', 1, 'Saturation', 1,...
                    'SaturationField', [], 'FixedBeta_Flag', 0, 'FixedBeta_Val', -1.5,...
                    'PoleSaturation', 0, 'PoleData', [], 'Trim', 0, 'TrimField', [], 'FitEst',0);
        %
        % My custom processing for specific data sets
%         [Processed_Data(ii), Uncorrected_Data(ii), Noise_Data(ii), Fitted_Data(ii), Data_Parameters(ii,:), Processing_Parameters(ii,:), Basis_Coeffs(ii,:)]...
%             = Process_Hyst_Data(handles.All_Data(ii), handles.All_Data_Order(ii), 'Offset', 1, 'Drift', 0, 'Saturation', 0,...
%             'SaturationField', 550, 'FixedBeta_Flag', 0, 'FixedBeta_Val', -1.5,...
%             'PoleSaturation', 0, 'PoleData', [], 'Trim', 1, 'TrimField', 990, 'FitEst',0);
        
% TODO - Tidy up this section


        catch
            % TODO - add more sophisticated error catch to allow other data
            % to be used even though only some cause problems
            
            % Delete the processing waitbar if in debug mode
            %             delete(h);
            
            Bad_Spec = [Bad_Spec; handles.All_Names(ii) ]; %#ok<AGROW>
            
            a = size(handles.All_Data{ii},1);

            Processed_Data(ii) = {NaN(a, 6)};
            Uncorrected_Data(ii) = {NaN(floor(a/2), 7)};
            Noise_Data(ii) = {NaN(floor(a/2), 3)};
            Fitted_Data(ii) = {NaN(floor(a/2), 6)};
            Data_Parameters(ii,:) = NaN(1,51);
            Processing_Parameters(ii,:) = zeros(1,14);
            Basis_Coeffs(ii,:) = cell(1,2);
            
            % Adjust some processing parameters to show no processing
            Processing_Parameters(ii,2) = 1; % Drift correction
            Processing_Parameters(ii,7) = 1; % Drift correction
            
        end % end try
        
        
    end
    
    % Warn of bad specimens
    if ~isempty(Bad_Spec)
        
        Error_MSG = {'The following specimens could not be processed:'};
        
        for ii = 1:length(Bad_Spec)
            Error_MSG(ii+1,1) = Bad_Spec(ii);
        end
        
        Error_MSG(ii+2,1) = {'Please check the raw data plots, remove any bad data, and try reprocessing.'};
        
        warndlg(Error_MSG, 'Processing Failed', 'modal');
        
    end
    
elseif Type == 1
    % Use the given processing parameters
    
    Processing_Parameters = handles.All_Processing_Parameters;
    
    [nSpec, nParams] = size(Processing_Parameters);
    
    % Check the size is correct and pad by NaNs if needed
    % Used for development and debugging
%     if nParams < 14
%         Processing_Parameters = [Processing_Parameters, NaN(nSpec, 13-nParams)];
%     end
    
    
    Processed_Data = cell(nSpec,1);
    Uncorrected_Data = cell(nSpec,1);
    Noise_Data = cell(nSpec,1);
    Fitted_Data = cell(nSpec,1);
    Data_Parameters = NaN(size(handles.All_Data_Parameters));
    Basis_Coeffs = cell(nSpec,2);
    
    
    for ii = 1:nSpec
        
        % Check for Cancel button press
        if Cancel_Flag == 1 || getappdata(h,'Cancelled')
            Cancel_Flag = 1;
            break;
        end
        
        % Update the waitbar and continue
        waitbar((ii-1)/(nSpec), h, strcat(sprintf(' %d', ii-1), ' specimens processed'))
        
        % Get the processing parameters
        Drift_Flag = Processing_Parameters(ii,1);
        %             Drift_Type = Processing_Parameters(ii,2);
        %             Drift_Ratio = Processing_Parameters(ii,3);
        % Temp_Ratio = Processing_Parameters(ii,4);
        Saturation_Flag = Processing_Parameters(ii,5);
        Saturation_Field = Processing_Parameters(ii,6);
        %             Method = Processing_Parameters(ii,7);
        Pole_Flag = Processing_Parameters(ii,8);
        Offset_Flag = Processing_Parameters(ii,9);
        Trim_Flag = Processing_Parameters(ii,10);
        Trim_Field = Processing_Parameters(ii,11);
        Fit_Param_Flag = Processing_Parameters(ii,12);
        FixedBeta_Flag = Processing_Parameters(ii,13);
        FixedBeta_Val = Processing_Parameters(ii,14);
        
        [Processed_Data(ii), Uncorrected_Data(ii), Noise_Data(ii), Fitted_Data(ii), Data_Parameters(ii,:), Processing_Parameters(ii,:), Basis_Coeffs(ii,:)]...
            = Process_Hyst_Data(handles.All_Data(ii), handles.All_Data_Order(ii), 'Offset', Offset_Flag, 'Drift', Drift_Flag, 'Saturation', Saturation_Flag,...
            'SaturationField', Saturation_Field, 'FixedBeta_Flag', FixedBeta_Flag, 'FixedBeta_Val', FixedBeta_Val,...
            'PoleSaturation', Pole_Flag, 'PoleData', [], 'Trim', Trim_Flag, 'TrimField', Trim_Field, 'FitEst', Fit_Param_Flag);
        % keyboard
        
    end
    
else
    error('HystLab:Initial_Processing', 'Unrecognized type flag.');
end


delete(h) % delete the waitbar

handles.Cancel_Processing = 0;

% Check for errors
if Error_Flag == 1
    warndlg(Error_MSG, 'Processing Failed', 'modal');
    handles.Cancel_Processing = 1;
    handles.Data_Loaded = 0;
    return;
end

% Check for user canceling
if Cancel_Flag == 1
    handles.Cancel_Processing = 1;
    handles.Data_Loaded = 0;
    return;
end


handles.All_Processed_Data = Processed_Data;
handles.All_Uncorrected_Data = Uncorrected_Data;
handles.All_Noise_Data = Noise_Data;
handles.All_Fitted_Data = Fitted_Data;
handles.All_Data_Parameters = Data_Parameters;
handles.All_Processing_Parameters = Processing_Parameters;
handles.All_Basis_Coeffs = Basis_Coeffs;


if handles.Data_Loaded == 0
    % Set the starting values
    handles.spec_ind = 1; % Set the specimen index to 1
    handles.Nspec=length(handles.All_Names); % Get the number of specimens
    set(handles.Spec_Num, 'String', handles.All_Names{handles.spec_ind}); % set the index
    handles.Norm_Type = 'None'; % PLot unnormalized data
    % guidata(hObject, handles);
    
    % Set the current data
    handles = Set_Current_Data(handles);
else
    % Data are loaded and we just stay on the current plot
    % Should be ok to do nothing
end


% --- Set the data for the currently loaded sample
function handles = Set_Current_Data(handles)

% Get the data into the current handles
handles.Current_Raw_Loop = handles.All_Data{handles.spec_ind};
handles.Current_Uncorrected_Data = handles.All_Uncorrected_Data{handles.spec_ind};
handles.Current_Processed_Data = handles.All_Processed_Data{handles.spec_ind};
handles.Current_Fitted_Data = handles.All_Fitted_Data{handles.spec_ind};
handles.Current_Noise_Data = handles.All_Noise_Data{handles.spec_ind};
handles.Current_Data_Parameters = handles.All_Data_Parameters(handles.spec_ind,:);
handles.Current_Spec_Mass = handles.All_Masses(handles.spec_ind);

% handles.Current_Processing_Parameters = handles.All_Processing_Parameters(handles.spec_ind,:);

% Update the information handles
spec_ind = handles.spec_ind;

set(handles.Spec_Count, 'String', strcat('Specimen: ', sprintf('%d', spec_ind)));

Slope_Corr_Msg = {'No correction applied', 'Linear high-field applied', 'Approach to saturation applied'};
Drift_Corr_Msg = {'No correction applied', 'Automatic', 'Positive field correction', 'Upper branch correction',...
    'Symmetric averaging', 'Paramagnetic', '', 'Paramagnetic/Postive', 'Paramagnetic/Upper', 'Paramagnetic/Symmetric'};

%                                   1            2                          4             5              6           7
% Processing_Parameters(ii,:) = [Drift_Flag, Drift_Type, Drift_Ratio Temp_Ratio Saturation_Flag, Saturation_Field, Method(ii),
%      8           9        10            11           12            13              14
% Pole_Flag, Offset_Flag, Trim_flag, Trim_Field, Fit_Param_Flag, FixedBeta_Flag, FixedBeta_Val];

set(handles.Drift_Corr_Menu, 'Value', handles.All_Processing_Parameters(spec_ind, 1)+1)
set(handles.Drift_Corr_MSG, 'String', Drift_Corr_Msg{handles.All_Processing_Parameters(spec_ind, 2)});
set(handles.Drift_Ratio, 'String', sprintf('%2.2f',handles.All_Processing_Parameters(spec_ind, 3)), 'Value', handles.All_Processing_Parameters(spec_ind, 3))
set(handles.Temp_Ratio, 'String', sprintf('%1.3f',handles.All_Processing_Parameters(spec_ind, 4)), 'Value', handles.All_Processing_Parameters(spec_ind, 4))

set(handles.Slope_Corr_Field, 'String', sprintf('%3.1f', handles.All_Processing_Parameters(spec_ind, 6)), 'Value', handles.All_Processing_Parameters(spec_ind, 6));
set(handles.Slope_Corr_Menu, 'Value', handles.All_Processing_Parameters(spec_ind, 5)+1);
set(handles.Slope_Corr_MSG, 'String', Slope_Corr_Msg{handles.All_Processing_Parameters(spec_ind, 7)});

set(handles.CB_Pole_Sat, 'Value', handles.All_Processing_Parameters(spec_ind, 8));

set(handles.CB_Offset_Corr, 'Value', handles.All_Processing_Parameters(spec_ind, 9))
set(handles.CB_Trim_Fields, 'Value', handles.All_Processing_Parameters(spec_ind, 10));
set(handles.Trim_Field_Lim, 'String', sprintf('%3.1f', handles.All_Processing_Parameters(spec_ind,11)), 'Value', handles.All_Processing_Parameters(spec_ind, 11));

set(handles.CB_Use_Fit_Params, 'Value', handles.All_Processing_Parameters(spec_ind, 12));

set(handles.CB_Fixed_Beta, 'Value', handles.All_Processing_Parameters(spec_ind, 13));
set(handles.Fixed_Beta_Val, 'String', sprintf('%1.3f', handles.All_Processing_Parameters(spec_ind,14)), 'Value', handles.All_Processing_Parameters(spec_ind, 14));



%                    1                       7             10   11  12
% Data_Parameters = [H0, M0, Q, Qf, Ms, Mrs, Bc, Brh, Bih, Xhf, X0, Shape,...
%  13     14    15     16    17  18    19    20       21   22   23     24
% C_err, Nfit, alpha, beta, p60, p70, p80, pCurrent, F60, F70, F80, FCurrent...
%  25        26   27       28       29       30   31  32    33        34   35  36     37
% S_star, Fit_F, Fit_p, Linear_F, Linear_p, p60, p70, p80, pCurrent, F60, F70, F80, FCurrent];
%  38    39   40     41         42        43      44-47    48-51
%   Qrh  Qih  Bsat   Sat_pct  Noise_RMS Fit_RMS RH_SNR_HN RH_SNR_HL

Data_Stats = handles.Current_Data_Parameters;

% Set the fitting stats
set(handles.Stat_Bo, 'String', sprintf('%3.2f', Data_Stats(1)));
set(handles.Stat_Mo, 'String', sprintf('%2.2e', Data_Stats(2)));
set(handles.Stat_Q, 'String', sprintf('%2.2f', Data_Stats(3)));
set(handles.Stat_Qf, 'String', sprintf('%2.2f', Data_Stats(4)));
set(handles.Stat_Cerr, 'String', sprintf('%2.2e', Data_Stats(13)));
set(handles.Stat_Nfit, 'String', sprintf('%d', Data_Stats(14)));
set(handles.Stat_Fit_F, 'String', sprintf('%3.1f', Data_Stats(26)));
set(handles.Stat_Fit_p, 'String', sprintf('%3.1f', 100.*Data_Stats(27)));

set(handles.Stat_Linear_F, 'String', sprintf('%3.1f', Data_Stats(28)));
set(handles.Stat_Linear_p, 'String', sprintf('%3.1f', 100.*Data_Stats(29)));

set(handles.Stat_Qrh, 'String', sprintf('%2.2f', Data_Stats(38)));
set(handles.Stat_Qih, 'String', sprintf('%2.2f', Data_Stats(39)));



% Set the Magnetic properties
set(handles.Stat_Ms, 'String', sprintf('%2.2e', Data_Stats(5)));
set(handles.Stat_Mrs, 'String', sprintf('%2.2e', Data_Stats(6)));
set(handles.Stat_Bc, 'String', sprintf('%3.1f', Data_Stats(7)));
set(handles.Stat_Brh, 'String', sprintf('%3.1f', Data_Stats(8)));
set(handles.Stat_Bih, 'String', sprintf('%3.1f', Data_Stats(9)));
set(handles.Stat_Xhf, 'String', sprintf('%2.2e', Data_Stats(10)));
set(handles.Stat_Shape, 'String', sprintf('%2.3f', Data_Stats(12)));
set(handles.Stat_Ss, 'String', sprintf('%1.4f', Data_Stats(25)));
set(handles.Stat_Bsat, 'String', sprintf('%3.1f', Data_Stats(40)));
set(handles.Stat_Sat_pct, 'String', sprintf('%3.1f', Data_Stats(41)));

% Set the approach to saturation stats
set(handles.Stat_NL_alpha, 'String', sprintf('%2.2e', Data_Stats(15)));
set(handles.Stat_NL_beta, 'String', sprintf('%2.2f', Data_Stats(16)));

% Check which stats are selected
switch get(handles.Saturation_Stat_Toggle, 'String')
    
    case 'Lack of Fit'
        
        % Set the data
        set(handles.Stat_p_NL_70, 'String', sprintf('%3.1f', 100.*Data_Stats(17)));
        set(handles.Stat_p_NL_80, 'String', sprintf('%3.1f', 100.*Data_Stats(18)));
        set(handles.Stat_p_NL_90, 'String', sprintf('%3.1f', 100.*Data_Stats(19)));
        set(handles.Stat_p_NL_Current, 'String', sprintf('%3.1f', 100.*Data_Stats(20)));
        set(handles.Stat_F_NL_70, 'String', sprintf('%3.1f', Data_Stats(21)));
        set(handles.Stat_F_NL_80, 'String', sprintf('%3.1f', Data_Stats(22)));
        set(handles.Stat_F_NL_90, 'String', sprintf('%3.1f', Data_Stats(23)));
        set(handles.Stat_F_NL_Current, 'String', sprintf('%3.1f', Data_Stats(24)));
        
        % Set the names
        set(handles.Stat_Ca, 'String', 'F')
        set(handles.Stat_70a, 'String', 'F70')
        set(handles.Stat_80a, 'String', 'F80')
        set(handles.Stat_90a, 'String', 'F90')
        set(handles.Stat_Cb, 'String', 'p')
        set(handles.Stat_70b, 'String', 'p70')
        set(handles.Stat_80b, 'String', 'p80')
        set(handles.Stat_90b, 'String', 'p90')
        
        % Set the units
        set(handles.Unit_Ca, 'String', '')
        set(handles.Unit_70a, 'String', '')
        set(handles.Unit_80a, 'String', '')
        set(handles.Unit_90a, 'String', '')
        set(handles.Unit_Cb, 'String', '%')
        set(handles.Unit_70b, 'String', '%')
        set(handles.Unit_80b, 'String', '%')
        set(handles.Unit_90b, 'String', '%')
                
    case 'Model Comparison'
        
        % Set the data
        set(handles.Stat_p_NL_70, 'String', sprintf('%3.1f', 100.*Data_Stats(30)));
        set(handles.Stat_p_NL_80, 'String', sprintf('%3.1f', 100.*Data_Stats(31)));
        set(handles.Stat_p_NL_90, 'String', sprintf('%3.1f', 100.*Data_Stats(32)));
        set(handles.Stat_p_NL_Current, 'String', sprintf('%3.1f', 100.*Data_Stats(33)));
        set(handles.Stat_F_NL_70, 'String', sprintf('%3.1f', Data_Stats(34)));
        set(handles.Stat_F_NL_80, 'String', sprintf('%3.1f', Data_Stats(35)));
        set(handles.Stat_F_NL_90, 'String', sprintf('%3.1f', Data_Stats(36)));
        set(handles.Stat_F_NL_Current, 'String', sprintf('%3.1f', Data_Stats(37)));
                
        % Set the names
        set(handles.Stat_Ca, 'String', 'F')
        set(handles.Stat_70a, 'String', 'F70')
        set(handles.Stat_80a, 'String', 'F80')
        set(handles.Stat_90a, 'String', 'F90')
        set(handles.Stat_Cb, 'String', 'p')
        set(handles.Stat_70b, 'String', 'p70')
        set(handles.Stat_80b, 'String', 'p80')
        set(handles.Stat_90b, 'String', 'p90')
        
        % Set the units
        set(handles.Unit_Ca, 'String', '')
        set(handles.Unit_70a, 'String', '')
        set(handles.Unit_80a, 'String', '')
        set(handles.Unit_90a, 'String', '')
        set(handles.Unit_Cb, 'String', '%')
        set(handles.Unit_70b, 'String', '%')
        set(handles.Unit_80b, 'String', '%')
        set(handles.Unit_90b, 'String', '%')
        
    case 'Loop Closure'
        
        % Set the data
        set(handles.Stat_p_NL_70, 'String', sprintf('%3.1f', Data_Stats(44)));
        set(handles.Stat_p_NL_80, 'String', sprintf('%3.1f', Data_Stats(45)));
        set(handles.Stat_p_NL_90, 'String', sprintf('%3.1f', Data_Stats(46)));
        set(handles.Stat_p_NL_Current, 'String', sprintf('%3.1f', Data_Stats(47)));
        set(handles.Stat_F_NL_70, 'String', sprintf('%3.1f', Data_Stats(48)));
        set(handles.Stat_F_NL_80, 'String', sprintf('%3.1f', Data_Stats(49)));
        set(handles.Stat_F_NL_90, 'String', sprintf('%3.1f', Data_Stats(50)));
        set(handles.Stat_F_NL_Current, 'String', sprintf('%3.1f', Data_Stats(51)));
        
        % Set the names
        set(handles.Stat_Ca, 'String', 'SNR')
        set(handles.Stat_70a, 'String', '70')
        set(handles.Stat_80a, 'String', '80')
        set(handles.Stat_90a, 'String', '90')
        set(handles.Stat_Cb, 'String', 'HAR')
        set(handles.Stat_70b, 'String', '70')
        set(handles.Stat_80b, 'String', '80')
        set(handles.Stat_90b, 'String', '90')
        
        % Set the units
        set(handles.Unit_Ca, 'String', 'dB')
        set(handles.Unit_70a, 'String', 'dB')
        set(handles.Unit_80a, 'String', 'dB')
        set(handles.Unit_90a, 'String', 'dB')
        set(handles.Unit_Cb, 'String', 'dB')
        set(handles.Unit_70b, 'String', 'dB')
        set(handles.Unit_80b, 'String', 'dB')
        set(handles.Unit_90b, 'String', 'dB')
    
end

% Check and deactivate any active boxes and vice versa

if any(get(handles.Slope_Corr_Menu, 'Value')==[1,2])
    set(handles.Slope_Corr_Field, 'Enable', 'off');
else
    set(handles.Slope_Corr_Field, 'Enable', 'on');
end

if get(handles.CB_Trim_Fields, 'Value') == 1
    set(handles.Trim_Field_Lim, 'Enable', 'on');
else
    set(handles.Trim_Field_Lim, 'Enable', 'off');
end


if strcmpi(Slope_Corr_Msg{handles.All_Processing_Parameters(spec_ind, 7)}, 'Approach to saturation applied')
    set(handles.CB_Fixed_Beta, 'Enable', 'on');
else
    set(handles.CB_Fixed_Beta, 'Enable', 'off');
end

% Highlight model estimated params
if get(handles.CB_Use_Fit_Params, 'Value') == 1
    
    set(handles.Stat_Ms, 'FontWeight', 'Bold', 'FontAngle', 'Italic');
    set(handles.Stat_Mrs, 'FontWeight', 'Bold', 'FontAngle', 'Italic');
    set(handles.Stat_Bc, 'FontWeight', 'Bold', 'FontAngle', 'Italic');
    set(handles.Stat_Brh, 'FontWeight', 'Bold', 'FontAngle', 'Italic');
    set(handles.Stat_Bih, 'FontWeight', 'Bold', 'FontAngle', 'Italic');
    set(handles.Stat_Xhf, 'FontWeight', 'Bold', 'FontAngle', 'Italic');
    set(handles.Stat_Shape, 'FontWeight', 'Bold', 'FontAngle', 'Italic');
    set(handles.Stat_Ss, 'FontWeight', 'Bold', 'FontAngle', 'Italic');
    
    set(handles.Stat_Bsat, 'FontWeight', 'Bold', 'FontAngle', 'Italic');
    set(handles.Stat_Sat_pct, 'FontWeight', 'Bold', 'FontAngle', 'Italic');
    
    % Set the approach to saturation stats
    set(handles.Stat_NL_alpha, 'FontWeight', 'Bold', 'FontAngle', 'Italic');
    set(handles.Stat_NL_beta, 'FontWeight', 'Bold', 'FontAngle', 'Italic');
    
    set(handles.Stat_p_NL_70, 'FontWeight', 'Bold', 'FontAngle', 'Italic');
    set(handles.Stat_p_NL_80, 'FontWeight', 'Bold', 'FontAngle', 'Italic');
    set(handles.Stat_p_NL_90, 'FontWeight', 'Bold', 'FontAngle', 'Italic');
    set(handles.Stat_p_NL_Current, 'FontWeight', 'Bold', 'FontAngle', 'Italic');
    set(handles.Stat_F_NL_70, 'FontWeight', 'Bold', 'FontAngle', 'Italic');
    set(handles.Stat_F_NL_80, 'FontWeight', 'Bold', 'FontAngle', 'Italic');
    set(handles.Stat_F_NL_90, 'FontWeight', 'Bold', 'FontAngle', 'Italic');
    set(handles.Stat_F_NL_Current, 'FontWeight', 'Bold', 'FontAngle', 'Italic');
    
    
else
    set(handles.Stat_Ms, 'FontWeight', 'Normal', 'FontAngle', 'Normal');
    set(handles.Stat_Mrs, 'FontWeight', 'Normal', 'FontAngle', 'Normal');
    set(handles.Stat_Bc, 'FontWeight', 'Normal', 'FontAngle', 'Normal');
    set(handles.Stat_Brh, 'FontWeight', 'Normal', 'FontAngle', 'Normal');
    set(handles.Stat_Bih, 'FontWeight', 'Normal', 'FontAngle', 'Normal');
    set(handles.Stat_Xhf, 'FontWeight', 'Normal', 'FontAngle', 'Normal');
    set(handles.Stat_Shape, 'FontWeight', 'Normal', 'FontAngle', 'Normal');
    set(handles.Stat_Ss, 'FontWeight', 'Normal', 'FontAngle', 'Normal');
    
    set(handles.Stat_Bsat, 'FontWeight', 'Normal', 'FontAngle', 'Normal');
    set(handles.Stat_Sat_pct, 'FontWeight', 'Normal', 'FontAngle', 'Normal');
    
    % Set the approach to saturation stats
    set(handles.Stat_NL_alpha, 'FontWeight', 'Normal', 'FontAngle', 'Normal');
    set(handles.Stat_NL_beta, 'FontWeight', 'Normal', 'FontAngle', 'Normal');
    
    set(handles.Stat_p_NL_70, 'FontWeight', 'Normal', 'FontAngle', 'Normal');
    set(handles.Stat_p_NL_80, 'FontWeight', 'Normal', 'FontAngle', 'Normal');
    set(handles.Stat_p_NL_90, 'FontWeight', 'Normal', 'FontAngle', 'Normal');
    set(handles.Stat_p_NL_Current, 'FontWeight', 'Normal', 'FontAngle', 'Normal');
    set(handles.Stat_F_NL_70, 'FontWeight', 'Normal', 'FontAngle', 'Normal');
    set(handles.Stat_F_NL_80, 'FontWeight', 'Normal', 'FontAngle', 'Normal');
    set(handles.Stat_F_NL_90, 'FontWeight', 'Normal', 'FontAngle', 'Normal');
    set(handles.Stat_F_NL_Current, 'FontWeight', 'Normal', 'FontAngle', 'Normal');
    
end


% Check for masses and normalization options

set(handles.Spec_Mass, 'Value', handles.Current_Spec_Mass, 'String', sprintf('%3.2f', handles.Current_Spec_Mass) );

if isnan(handles.Current_Spec_Mass)
    set(handles.RB_Mass_Norm, 'Enable', 'off');
    
    if get(handles.RB_Mass_Norm, 'Value') == 1
        % Switch to default no normalization
        set(handles.RB_No_Norm, 'Value', 1);
        handles.Norm_Type = 'None';
    end
    
else
    % Enable the mass normalization button
    set(handles.RB_Mass_Norm, 'Enable', 'on');
end


% Update the mass specific parameters

Data_Stats = handles.Current_Data_Parameters;

switch handles.Norm_Type
    
    case 'None'
        Normalizer = 1;
        
        % Set the unit text
        set(handles.Txt_Mag_U1, 'String', 'Am2')
        set(handles.Txt_Mag_U2, 'String', 'Am2')
        set(handles.Txt_Mag_U3, 'String', 'Am2')
        set(handles.Txt_Mag_U4, 'String', 'Am2')
        set(handles.Txt_Xhf_U1, 'String', 'm3')
        
        
    case 'Ms'
        Normalizer = Data_Stats(5);
        
        % Set the unit text
        set(handles.Txt_Mag_U1, 'String', '')
        set(handles.Txt_Mag_U2, 'String', '')
        set(handles.Txt_Mag_U3, 'String', '')
        set(handles.Txt_Mag_U4, 'String', '')
        set(handles.Txt_Xhf_U1, 'String', '')
        
    case 'Mass'
        Normalizer = handles.Current_Spec_Mass/1e6; %in kg
        
        % Set the unit text
        set(handles.Txt_Mag_U1, 'String', 'Am2/kg')
        set(handles.Txt_Mag_U2, 'String', 'Am2/kg')
        set(handles.Txt_Mag_U3, 'String', 'Am2/kg')
        set(handles.Txt_Mag_U4, 'String', 'Am2/kg')
        set(handles.Txt_Xhf_U1, 'String', 'm3/kg')
end

set(handles.Stat_Mo, 'String', sprintf('%2.2e', Data_Stats(2)/Normalizer));
set(handles.Stat_Cerr, 'String', sprintf('%2.2e', Data_Stats(13)/Normalizer));
set(handles.Stat_Ms, 'String', sprintf('%2.2e', Data_Stats(5)/Normalizer));
set(handles.Stat_Mrs, 'String', sprintf('%2.2e', Data_Stats(6)/Normalizer));
set(handles.Stat_Xhf, 'String', sprintf('%2.2e', Data_Stats(10)/Normalizer));


try
    % Update the plots
    Update_Plots(handles)
catch
    % May fail if data are not properly loaded
end


% --- Update the plot axes
function Update_Plots(handles)

% Get the colors
Hyst_Color = handles.Default_Hyst_Color;
IR_Color = handles.Default_IR_Color;
Noise_Color = handles.Default_Noise_Color;

% Get the data from the handles for easy management
Raw_Fields = handles.Current_Raw_Loop(:,1);
Raw_Mirh_Fields = handles.Current_Uncorrected_Data(:,1);
Processed_Fields = handles.Current_Processed_Data(:,1:2);
Noise_Fields = handles.Current_Noise_Data(:,1);
Fitted_Fields = handles.Current_Fitted_Data(:,1:2);


switch handles.Norm_Type
    case 'None'
        % Do nothing (or normalize by 1)
        Normalizer = 1;
    case 'Ms'
        % Normalize by Ms
        Normalizer = handles.Current_Data_Parameters(5);
    case 'Mass'
        Normalizer = handles.Current_Spec_Mass .* 1e-6; % Convert to kg
end


% Get the data after normalization
Raw_Loop = handles.Current_Raw_Loop(:,2) ./ Normalizer;
Noise_Moment = handles.Current_Noise_Data(:,2) ./ Normalizer;
Noise_Smooth = handles.Current_Noise_Data(:,3) ./ Normalizer;
Noise_RMS = handles.Current_Data_Parameters(42) ./Normalizer;

Raw_Mih_Moment = handles.Current_Uncorrected_Data(:,5) ./ Normalizer;
Raw_Mrh_Moment = handles.Current_Uncorrected_Data(:,6) ./ Normalizer;
Raw_Noise = handles.Current_Uncorrected_Data(:,7) ./ Normalizer;


Processed_Loop = handles.Current_Processed_Data(:,3:4) ./ Normalizer;
Processed_Mih = handles.Current_Processed_Data(:,5) ./ Normalizer;
Processed_Mrh = handles.Current_Processed_Data(:,6) ./ Normalizer;

Fitted_Loop = handles.Current_Fitted_Data(:,3:4) ./ Normalizer;
Fitted_Mih = handles.Current_Fitted_Data(:,5) ./ Normalizer;
Fitted_Mrh = handles.Current_Fitted_Data(:,6) ./ Normalizer;


% The limits for plotting the origin cross
Cross_YLim = 1.1*max(max( abs([ Raw_Loop; Processed_Loop(:) ]) ));
Cross_XLim = 1.1*max(max( abs([ Raw_Fields; Processed_Fields(:) ]) ));


% Do the plots
FUnits = 'Pixels';
FontSize1 = 12; % 12pt font
FontSize2 = 14; % 14pt font

% Do the loop plot

% reset the axes
% reset(handles.Loop_Axes);
cla(handles.Loop_Axes);
hold(handles.Loop_Axes, 'on')

if get(handles.CB_Origin_Cross, 'Value') == 1
    plot(handles.Loop_Axes, [0, 0], [-Cross_YLim, Cross_YLim], 'k', 'LineWidth', 1);
    plot(handles.Loop_Axes, [-Cross_XLim, Cross_XLim], [0, 0], 'k', 'LineWidth', 1);
    
end


if get(handles.CB_Plot_Raw, 'Value') == 1
    plot(handles.Loop_Axes, Raw_Fields, Raw_Loop, 'Color', Hyst_Color(1,:), 'LineWidth', 1);
end

if get(handles.CB_Plot_Processed, 'Value') == 1
    plot(handles.Loop_Axes, Processed_Fields, Processed_Loop, 'Color', Hyst_Color(2,:), 'LineWidth', 1)
end

if get(handles.CB_Plot_Fitted, 'Value') == 1
    plot(handles.Loop_Axes, Fitted_Fields, Fitted_Loop, 'Color', Hyst_Color(3,:), 'LineWidth', 1);
end

hold(handles.Loop_Axes, 'off')

% Set the labels
set(get(handles.Loop_Axes, 'XLabel'), 'String', 'Field [mT]', 'FontUnits', FUnits, 'FontSize', FontSize1);
set(get(handles.Loop_Axes, 'YLabel'), 'String', 'Moment [Am^2]', 'FontUnits', FUnits, 'FontSize', FontSize1);
set(get(handles.Loop_Axes, 'Title'), 'String', 'Hysteresis Loop', 'FontUnits', FUnits, 'FontSize', FontSize2, 'FontWeight', 'bold', 'FontName', 'Default');
set(handles.Loop_Axes, 'Box', 'on')


% Do the Mih and Mrh plot

% Reset the axes
reset(handles.IR_Axes);
cla(handles.IR_Axes)

if get(handles.CB_irh_Gradient, 'Value') == 1
    % Plot the gradients
    
    % Get the field step converted to A/m
    Fields1 = 1e4 .* Processed_Fields(:,1) ./ (4*pi) ;
    Fields2 =  1e4 .* Fitted_Fields(:,1)  ./ (4*pi) ;
    
    
    hold(handles.IR_Axes, 'on')
    
    if get(handles.CB_Origin_Cross, 'Value') == 1
        plot(handles.IR_Axes, [0, 0], [-Cross_YLim, Cross_YLim], 'k', 'LineWidth', 1);
        plot(handles.IR_Axes, [-Cross_XLim, Cross_XLim], [0, 0], 'k', 'LineWidth', 1);
    end
    
    if get(handles.CB_Plot_Raw, 'Value') == 1
        
        % regularize the field spacing and take the gradient
        %         grad_Mih = gradient(interp1(Raw_Mirh_Fields, Raw_Mih_Moment, Processed_Fields(:,1), 'phicp'), dF1);
        %         grad_Mrh = gradient(interp1(Raw_Mirh_Fields, Raw_Mrh_Moment, Processed_Fields(:,1), 'phicp'), dF1);
        
        grad_Mih = gradient(Interpolate_To_Field(Raw_Mirh_Fields, Raw_Mih_Moment, Processed_Fields(:,1), 'linear'), Fields1);
        grad_Mrh = gradient(Interpolate_To_Field(Raw_Mirh_Fields, Raw_Mrh_Moment, Processed_Fields(:,1), 'linear'), Fields1);
        
        plot(handles.IR_Axes, Processed_Fields(:,1), grad_Mih, 'Color', IR_Color(1,:), 'LineWidth', 1);
        plot(handles.IR_Axes, Processed_Fields(:,1), grad_Mrh, 'Color', IR_Color(1,:), 'LineWidth', 1);
        
    end
    
    if get(handles.CB_Plot_Processed, 'Value') == 1
        plot(handles.IR_Axes, Processed_Fields(:,1), gradient(Processed_Mih, Fields1), 'Color', IR_Color(2,:), 'LineWidth', 1);
        plot(handles.IR_Axes, Processed_Fields(:,1), gradient(Processed_Mrh, Fields1), 'Color', IR_Color(2,:), 'LineWidth', 1);
    end
    
    
    if get(handles.CB_Plot_Fitted, 'Value') == 1
        plot(handles.IR_Axes, Fitted_Fields(:,1), gradient(Fitted_Mih, Fields2), 'Color', IR_Color(3,:), 'LineWidth', 1);
        plot(handles.IR_Axes, Fitted_Fields(:,1), gradient(Fitted_Mrh, Fields2), 'Color', IR_Color(3,:), 'LineWidth', 1);
%         keyboard
    end
    
    hold(handles.IR_Axes, 'off')
    
    % Set the labels
    set(get(handles.IR_Axes, 'XLabel'), 'String', 'Field [mT]', 'FontUnits', FUnits, 'FontSize', FontSize1);
    set(get(handles.IR_Axes, 'YLabel'), 'String', 'Moment [Am^2]', 'FontUnits', FUnits, 'FontSize', FontSize1);
    set(get(handles.IR_Axes, 'Title'), 'String', 'M_{ih} & M_{rh} Gradients', 'FontUnits', FUnits, 'FontSize', FontSize2, 'FontWeight', 'bold');
    set(handles.IR_Axes, 'Box', 'on')
    
else
    % Plot the normal curves
    
    hold(handles.IR_Axes, 'on')
    
    if get(handles.CB_Origin_Cross, 'Value') == 1
        plot(handles.IR_Axes, [0, 0], [-Cross_YLim, Cross_YLim], 'k', 'LineWidth', 1);
        plot(handles.IR_Axes, [-Cross_XLim, Cross_XLim], [0, 0], 'k', 'LineWidth', 1);
    end
    
    if get(handles.CB_Plot_Raw, 'Value') == 1
        plot(handles.IR_Axes, Raw_Mirh_Fields, Raw_Mih_Moment, 'Color', IR_Color(1,:), 'LineWidth', 1);
        plot(handles.IR_Axes, Raw_Mirh_Fields, Raw_Mrh_Moment, 'Color', IR_Color(1,:), 'LineWidth', 1);
    end
    
    if get(handles.CB_Plot_Processed, 'Value') == 1
        plot(handles.IR_Axes, Processed_Fields(:,1), Processed_Mih, 'Color', IR_Color(2,:), 'LineWidth', 1);
        plot(handles.IR_Axes, Processed_Fields(:,1), Processed_Mrh, 'Color', IR_Color(2,:), 'LineWidth', 1);
    end
    
    
    
    if get(handles.CB_Plot_Fitted, 'Value') == 1
        plot(handles.IR_Axes, Fitted_Fields(:,1), Fitted_Mih, 'Color', IR_Color(3,:), 'LineWidth', 1);
        plot(handles.IR_Axes, Fitted_Fields(:,1), Fitted_Mrh, 'Color', IR_Color(3,:), 'LineWidth', 1);
    end
    hold(handles.IR_Axes, 'off')
    
    % Set the labels
    set(get(handles.IR_Axes, 'XLabel'), 'String', 'Field [mT]', 'FontUnits', FUnits, 'FontSize', FontSize1);
    set(get(handles.IR_Axes, 'YLabel'), 'String', 'Moment [Am^2]', 'FontUnits', FUnits, 'FontSize', FontSize1);
    set(get(handles.IR_Axes, 'Title'), 'String', 'M_{ih} & M_{rh} Curves', 'FontUnits', FUnits, 'FontSize', FontSize2, 'FontWeight', 'bold');
    set(handles.IR_Axes, 'Box', 'on')
end


% Do the noise plot
reset(handles.Noise_Axes);
cla(handles.Noise_Axes);

hold(handles.Noise_Axes, 'on')

% Check for plotting the raw noise
if get(handles.CB_Plot_Raw, 'Value') == 1
    plot(handles.Noise_Axes, Raw_Mirh_Fields, Raw_Noise, 'Color', Noise_Color(1,:), 'LineWidth', 1);
end

% Check for plotting procesed noise
if get(handles.CB_Plot_Processed, 'Value') == 1
    plot(handles.Noise_Axes, Noise_Fields, Noise_Moment, 'Color', Noise_Color(2,:), 'LineWidth', 1);
end

% Check for plotting smoothed noise
if get(handles.CB_Plot_Smooth_Noise, 'Value') == 1
    plot(handles.Noise_Axes, Noise_Fields, Noise_Smooth, 'Color', Noise_Color(3,:), 'LineWidth', 2);
end

hold(handles.Noise_Axes, 'off')


set(get(handles.Noise_Axes, 'XLabel'), 'String', 'Field [mT]', 'FontUnits', FUnits, 'FontSize', FontSize1);
set(get(handles.Noise_Axes, 'YLabel'), 'String', 'Moment [Am^2]', 'FontUnits', FUnits, 'FontSize', FontSize1);
set(get(handles.Noise_Axes, 'Title'), 'String', 'Noise Curve', 'FontUnits', FUnits, 'FontSize', FontSize2, 'FontWeight', 'bold');
set(handles.Noise_Axes, 'Box', 'on')


% set field limits for the x-axes
% Raw data, Raw data, Processed, fitted
Max_Field = [max(abs(Raw_Fields)), max(abs(Raw_Mirh_Fields)), max(abs(Processed_Fields)), max(abs(Processed_Fields))];

if get(handles.CB_Plot_Raw, 'Value') == 0
    Max_Field(1:2) = [0,0];
end

if get(handles.CB_Plot_Processed, 'Value') == 0
    Max_Field(3) = 0;
end

if get(handles.CB_Plot_Fitted, 'Value') == 0
    Max_Field(4) = 0;
end

Max_Field = max(Max_Field);

if Max_Field ~= 0 % some data are plotted
    
    %     Field_Lim = 10*ceil((50+Max_Field)/10); % add 50 and round up to nearest 10 mT
    Field_Lim = 1.15*RoundField(Max_Field);
    
    set(handles.Loop_Axes, 'Xlim', [-Field_Lim, Field_Lim], 'FontUnits', FUnits, 'FontName', 'Default', 'FontSize', FontSize1);
    set(handles.IR_Axes, 'Xlim', [-Field_Lim, Field_Lim], 'FontUnits', FUnits, 'FontName', 'Default', 'FontSize', FontSize1);
    set(handles.Noise_Axes, 'Xlim', [-Field_Lim, Field_Lim], 'FontUnits', FUnits, 'FontName', 'Default', 'FontSize', FontSize1);
    
end

% Set moment limits
yscale = 1.05;

Loop_Max = [max(abs(Raw_Loop)), max(abs(Processed_Loop)), max(abs(Fitted_Loop))];

% keyboard
if get(handles.CB_Plot_Raw, 'Value') == 0
    Loop_Max(1) = 0;
end

if get(handles.CB_Plot_Processed, 'Value') == 0
    Loop_Max(3) = 0;
end

if get(handles.CB_Plot_Fitted, 'Value') == 0
    Loop_Max(4) = 0;
end

Loop_Lim = yscale * max(Loop_Max);

if Loop_Lim ~= 0
    set(handles.Loop_Axes, 'Ylim', [-Loop_Lim, Loop_Lim], 'FontUnits', FUnits, 'FontName', 'Default', 'FontSize', FontSize1);
end


% The Mih and Mrh plots
IR_Max = zeros(1,6);
IR_Min = zeros(1,6);

if get(handles.CB_irh_Gradient, 'Value') == 0
    IR_Lim(1) = -Loop_Lim;
    IR_Lim(2) = Loop_Lim;
else
    if get(handles.CB_Plot_Raw, 'Value') == 1
        IR_Min(1:2) = [min(grad_Mih), min(grad_Mrh)];
        IR_Max(1:2) = [max(abs(grad_Mih)), max(abs(grad_Mrh))];
    end
    
    if get(handles.CB_Plot_Processed, 'Value') == 1
        IR_Min(3:4) = [min(gradient(Processed_Mih, Fields1)), min(gradient(Processed_Mrh, Fields1))];
        IR_Max(3:4) = [max(abs(gradient(Processed_Mih, Fields1))), max(abs(gradient(Processed_Mrh, Fields1)))];
    end
    
    
    if get(handles.CB_Plot_Fitted, 'Value') == 1
        IR_Min(5:6) = [min(gradient(Fitted_Mih, Fields2)), min(gradient(Fitted_Mrh, Fields2))];
        IR_Max(5:6) = [max(abs(gradient(Fitted_Mih, Fields2))), max(abs(gradient(Fitted_Mrh, Fields2)))];
    end
    
    IR_Lim(1) = yscale * min(IR_Min);
    IR_Lim(2) = yscale * max(IR_Max);
    
end

if sum(IR_Lim==0) ~=2
    set(handles.IR_Axes, 'Ylim', [IR_Lim(1), IR_Lim(2)], 'FontUnits', FUnits, 'FontName', 'Default', 'FontSize', FontSize1);
end

% The noise plots
Noise_Max = [max(abs(Noise_Moment)), max(abs(Noise_Smooth)), max(abs(Raw_Noise))];

if get(handles.CB_Plot_Smooth_Noise, 'Value') == 0
    Noise_Max(2) = 0;
end

if get(handles.CB_Plot_Raw, 'Value') == 0
    Noise_Max(3) = 0;
end

Noise_Lim = (0.1+yscale) * max(Noise_Max);


if Noise_Lim ~= 0
    set(handles.Noise_Axes, 'Ylim', [-Noise_Lim, Noise_Lim], 'FontUnits', FUnits, 'FontName', 'Default', 'FontSize', FontSize1);
end

% Get the limits for rescaling the plots
% keyboard
% Loop plot
% TL = get(handles.Loop_Axes, 'YTickLabel');
% TV = get(handles.Loop_Axes, 'YTick');
% e = log10(TV(1) ./str2double(TL(1,:)));
% ymax = max(get(handles.Loop_Axes, 'YTick'));
% e = log10(abs(ymax(1)));
% Loop_exp = sign(e).*floor(abs(e));
% 
% Loop_Tick = get(handles.Loop_Axes, 'YTick');
% Loop_Label = [];
% for ii = 1:length(Loop_Tick)
%     Loop_Label = [Loop_Label; sprintf('% 2.1f', Loop_Tick(ii) ./ 10^Loop_exp)]; %#ok<AGROW>
% end
% set(handles.Loop_Axes, 'YTickLabel', Loop_Label);


% Mrh, Mih plot
% ymax = max(get(handles.IR_Axes, 'YTick'));
% e = log10(abs(ymax(1)));
% IR_exp = sign(e).*floor(abs(e));
% 
% IR_Tick = get(handles.IR_Axes, 'YTick');
% IR_Label = [];
% for ii = 1:length(IR_Tick)
%     IR_Label = [IR_Label; sprintf('% 2.1f', IR_Tick(ii) ./ 10^IR_exp)]; %#ok<AGROW>
% end
% set(handles.IR_Axes, 'YTickLabel', IR_Label);


% Noise curve
% ymax = max(get(handles.Noise_Axes, 'YTick'));
% e = log10(abs(ymax(1)));
% Noise_exp = sign(e).*floor(abs(e));
% 
% Noise_Tick = get(handles.Noise_Axes, 'YTick');
% Noise_Label = [];
% for ii = 1:length(Noise_Tick)
%     Noise_Label = [Noise_Label; sprintf('% 2.1f', Noise_Tick(ii) ./ 10^Noise_exp)]; %#ok<AGROW>
% end
% set(handles.Noise_Axes, 'YTickLabel', Noise_Label);




% Re-label for different normalizations
switch handles.Norm_Type
    case 'None'
        % already done
        
       
        
%         set(get(handles.Loop_Axes, 'YLabel'), 'String', strcat('Moment [', sprintf('\\times10^{%d}', Loop_exp), ' Am^2]'),...
%             'FontUnits', FUnits, 'FontSize', FontSize1);
%         set(get(handles.IR_Axes, 'YLabel'), 'String', strcat('Moment [', sprintf('\\times10^{%d}', IR_exp), ' Am^2]'),...
%             'FontUnits', FUnits, 'FontSize', FontSize1);
%         set(get(handles.Noise_Axes, 'YLabel'), 'String', strcat('Moment [', sprintf('\\times10^{%d}', Noise_exp), ' Am^2]'),...
%             'FontUnits', FUnits, 'FontSize', FontSize1);
        
        
        if get(handles.CB_irh_Gradient, 'Value') == 1
            % Am^2 per A/m
            set(get(handles.IR_Axes, 'YLabel'), 'String', 'Moment Gradient [m^3]', 'FontUnits', FUnits, 'FontSize', FontSize1);
        end
        
        % set the noise RMS text
        axes(handles.Noise_Axes);
        x = get(gca, 'xlim');
        x = 0.02*(x(2) - x(1)) + x(1);
        y = get(gca, 'ylim');
        y = y(2) - 0.04*(y(2) - y(1));
        
        text(x, y, ['RMS = ', sprintf('%1.2e', Noise_RMS), ' Am^2'], 'FontSize', 11)
        
        
        
    case 'Ms'
        
        set(get(handles.Loop_Axes, 'YLabel'), 'String', 'Moment / Ms', 'FontUnits', FUnits, 'FontSize', FontSize1);
        set(get(handles.IR_Axes, 'YLabel'), 'String', 'Moment / Ms', 'FontUnits', FUnits, 'FontSize', FontSize1);
        set(get(handles.Noise_Axes, 'YLabel'), 'String', 'Moment / Ms', 'FontUnits', FUnits, 'FontSize', FontSize1);
        
        if get(handles.CB_irh_Gradient, 'Value') == 1
            % Am^2 per A/m
            set(get(handles.IR_Axes, 'YLabel'), 'String', 'Moment Gradient [mA^{-1}]', 'FontUnits', FUnits, 'FontSize', FontSize1);
        end
        
        % set the noise RMS text
        axes(handles.Noise_Axes);
        x = get(gca, 'xlim');
        x = 0.02*(x(2) - x(1)) + x(1);
        y = get(gca, 'ylim');
        y = y(2) - 0.04*(y(2) - y(1));
        
        text(x, y, ['RMS = ', sprintf('%1.2e', Noise_RMS)], 'FontSize', 11)
        
        
    case 'Mass'
        
%         set(get(handles.Loop_Axes, 'YLabel'), 'String', strcat('Moment [', sprintf('\\times10^{%d}', Loop_exp), ' Am^2/kg]'),...
%             'FontUnits', FUnits, 'FontSize', FontSize1);
%                 set(get(handles.IR_Axes, 'YLabel'), 'String', strcat('Moment [', sprintf('\\times10^{%d}', IR_exp), ' Am^2/kg]'),...
%             'FontUnits', FUnits, 'FontSize', FontSize1);
%         set(get(handles.Noise_Axes, 'YLabel'), 'String', strcat('Moment [', sprintf('\\times10^{%d}', Noise_exp), ' Am^2/kg]'),...
%             'FontUnits', FUnits, 'FontSize', FontSize1);

        set(get(handles.Loop_Axes, 'YLabel'), 'String', 'Magnetization [Am^2/kg]', 'FontUnits', FUnits, 'FontSize', FontSize1);
        set(get(handles.IR_Axes, 'YLabel'), 'String', 'Magnetization [Am^2/kg]', 'FontUnits', FUnits, 'FontSize', FontSize1);
        set(get(handles.Noise_Axes, 'YLabel'), 'String', 'Magnetization [Am^2/kg]', 'FontUnits', FUnits, 'FontSize', FontSize1);
        
        if get(handles.CB_irh_Gradient, 'Value') == 1
            % Am^2/kg per A/m
            set(get(handles.IR_Axes, 'YLabel'), 'String', 'Magnetization Gradient [m^3/kg]', 'FontUnits', FUnits, 'FontSize', FontSize1);
        end
        
        % set the noise RMS text
        axes(handles.Noise_Axes);
        x = get(gca, 'xlim');
        x = 0.02*(x(2) - x(1)) + x(1);
        y = get(gca, 'ylim');
        y = y(2) - 0.04*(y(2) - y(1));
        
        text(x, y, ['RMS = ', sprintf('%1.2e', Noise_RMS), ' Am^2/kg'], 'FontSize', 11)
        
        
end

% Reset the button down functions
set(handles.Loop_Axes, 'ButtonDownFcn', {@Loop_Axes_ButtonDownFcn, handles});
set(handles.IR_Axes, 'ButtonDownFcn', {@IR_Axes_ButtonDownFcn, handles});
set(handles.Noise_Axes, 'ButtonDownFcn', {@Noise_Axes_ButtonDownFcn, handles});


if get(handles.CB_Plot_Grid, 'Value') == 1
    grid(handles.Loop_Axes, 'on')
    grid(handles.IR_Axes, 'on')
    grid(handles.Noise_Axes, 'on')
else
    grid(handles.Loop_Axes, 'off')
    grid(handles.IR_Axes, 'off')
    grid(handles.Noise_Axes, 'off')
end


% --- Recalculate the fits after parameter update
function [func_handles, reset_flag] = ReFit_Data(handles)

% disp('Refitting...')

% Get the input parameters
Drift_Flag = get(handles.Drift_Corr_Menu, 'Value') -1;
Sat_Flag = get(handles.Slope_Corr_Menu, 'Value') - 1;
Sat_Field = str2double(get(handles.Slope_Corr_Field, 'String'));
Offset_Flag = get(handles.CB_Offset_Corr, 'Value');
Trim_Flag = get(handles.CB_Trim_Fields, 'Value');
Trim_Field = str2double(get(handles.Trim_Field_Lim, 'String'));
Pole_Sat_Flag = get(handles.CB_Pole_Sat, 'Value');
Fit_Est = get(handles.CB_Use_Fit_Params, 'Value');


FixedBeta_Val = -1.5;
FixedBeta_Flag = 0;

% Check the fixed beta box is ticked and enabled
if get(handles.CB_Fixed_Beta, 'Value') == 1 && strcmpi(get(handles.CB_Fixed_Beta, 'Enable'), 'on')
    %     get(handles.CB_Fixed_Beta, 'Enable')
    FixedBeta_Val = str2double(get(handles.Fixed_Beta_Val, 'String'));
    FixedBeta_Flag = 1;
end

spec_ind = handles.spec_ind;

try
    
    [Processed_Data, Uncorrected_Data, Noise_Data, Fitted_Data, Data_Parameters, Processing_Parameters, Basis_Coeffs, Error_Flag] = Process_Hyst_Data({handles.Current_Raw_Loop}, handles.All_Data_Order(spec_ind),...
        'Offset', Offset_Flag, 'Drift', Drift_Flag, 'Saturation', Sat_Flag, 'SaturationField', Sat_Field, 'FixedBeta_Flag', FixedBeta_Flag, 'FixedBeta_Val', FixedBeta_Val,...
        'PoleSaturation', Pole_Sat_Flag, 'PoleData', [], 'Trim', Trim_Flag, 'TrimField', Trim_Field, 'FitEst', Fit_Est);
    
    
catch
    % TODO - add more sophisticated error catch to allow other data
    % to be used even though only some cause problems
    
    Error_MSG = [ {sprintf('%s could not be correctly processed.', handles.All_Names{spec_ind} )};...
        {'Please check the raw data plot, remove any bad data points, and try reprocessing.'}];
    
    warndlg(Error_MSG, 'Processing Failed', 'modal');
    
    Error_Flag = 1;
%     a = size(handles.All_Data{ii},1);
%     
%     Processed_Data(ii) = {NaN(a, 6)};
%     Uncorrected_Data(ii) = {NaN(floor(a/2), 7)};
%     Noise_Data(ii) = {NaN(floor(a/2), 3)};
%     Fitted_Data(ii) = {NaN(floor(a/2), 6)};
%     Data_Parameters(ii,:) = NaN(1,51);
%     Processing_Parameters(ii,:) = zeros(1,14);
%     Basis_Coeffs(ii,:) = cell(1,2);
%     
%     % Adjust some processing parameters to show no processing
%     Processing_Parameters(ii,2) = 1; % Drift correction
%     Processing_Parameters(ii,7) = 1; % Drift correction
    
end


%
if Error_Flag ~= 0
    % here this means that there was an error
    % We just return
    reset_flag = 1;
    func_handles = handles;
    
    return;
end


handles.All_Processed_Data(spec_ind) = Processed_Data;
handles.All_Uncorrected_Data(spec_ind) = Uncorrected_Data;
handles.All_Noise_Data(spec_ind) = Noise_Data;
handles.All_Fitted_Data(spec_ind) = Fitted_Data;
handles.All_Data_Parameters(spec_ind,:) = Data_Parameters;
handles.All_Processing_Parameters(spec_ind,:) = Processing_Parameters;
handles.All_Basis_Coeffs(spec_ind,:) = Basis_Coeffs;


% keyboard

% Update the current data
handles = Set_Current_Data(handles);

reset_flag = 0;
func_handles = handles;


% --- Executes on selection change in Slope_Corr_Menu.
function Slope_Corr_Menu_Callback(hObject, eventdata, handles)
% hObject    handle to Slope_Corr_Menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Slope_Corr_Menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Slope_Corr_Menu

if handles.Data_Loaded == 0
    % No data loaded so do nothing
    return;
end

% If no or automatic corrections are applied disable the field input
if get(hObject,'Value') == 1 || get(hObject,'Value') == 2
    set(handles.Slope_Corr_Field, 'Enable', 'off');
else
    set(handles.Slope_Corr_Field, 'Enable', 'on');
end

% Disable fixed beta unless approach to saturation is manually selected
if get(hObject,'Value') == 4
    set(handles.CB_Fixed_Beta, 'Enable', 'on');
    
    if get(handles.CB_Fixed_Beta, 'Value') == 1
        set(handles.Fixed_Beta_Val, 'Enable', 'on');
    else
        set(handles.Fixed_Beta_Val, 'Enable', 'off');
    end
    
else
    set(handles.CB_Fixed_Beta, 'Enable', 'off');
    set(handles.Fixed_Beta_Val, 'Enable', 'off');
end

% Call the ReFit function and update handles structure
handles = ReFit_Data(handles);
% Update_Plots(handles);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function Slope_Corr_Menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Slope_Corr_Menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function Slope_Corr_Field_Callback(hObject, eventdata, handles)
% hObject    handle to Slope_Corr_Field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Slope_Corr_Field as text
%        str2double(get(hObject,'String')) returns contents of Slope_Corr_Field as a double

% input should be a number
old_val = get(handles.Slope_Corr_Field, 'Value');

if isnan(str2double(get(hObject,'String') )) || str2double(get(hObject,'String') ) < 0
    set(handles.Slope_Corr_Field, 'Value', old_val, 'String', num2str(old_val) );
else
    set(handles.Slope_Corr_Field, 'Value', str2double(get(hObject,'String') ) );
end

if handles.Data_Loaded == 0
    % No data loaded so do nothing
    return;
end

% Call the ReFit function and update handles structure
[handles, reset_flag] = ReFit_Data(handles);

if reset_flag == 1
    % we have some error - reset to old vale
    set(handles.Slope_Corr_Field, 'Value', old_val, 'String', num2str(old_val) );
end

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function Slope_Corr_Field_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Slope_Corr_Field (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in CB_Fixed_Beta.
function CB_Fixed_Beta_Callback(hObject, eventdata, handles)
% hObject    handle to CB_Fixed_Beta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CB_Fixed_Beta

if handles.Data_Loaded == 0
    % No data loaded so do nothing
    return;
end

if get(hObject,'Value') == 1
    set(handles.Fixed_Beta_Val, 'Enable', 'on');
else
    set(handles.Fixed_Beta_Val, 'Enable', 'off');
end



% Call the ReFit function and update handles structure
handles = ReFit_Data(handles);
guidata(hObject, handles);


function Fixed_Beta_Val_Callback(hObject, eventdata, handles)
% hObject    handle to Fixed_Beta_Val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Fixed_Beta_Val as text
%        str2double(get(hObject,'String')) returns contents of Fixed_Beta_Val as a double

% input should be a number
old_val = get(handles.Fixed_Beta_Val, 'Value');

if isnan(str2double(get(hObject,'String') )) || str2double(get(hObject,'String') ) < -3 || str2double(get(hObject,'String') ) > 0
    set(handles.Fixed_Beta_Val, 'Value', old_val, 'String', num2str(old_val) );
else
    set(handles.Fixed_Beta_Val, 'Value', str2double(get(hObject,'String') ) );
end

if handles.Data_Loaded == 0
    % No data loaded so do nothing
    return;
end

% Call the ReFit function and update handles structure
[handles, reset_flag] = ReFit_Data(handles);

if reset_flag == 1
    % we have some error - reset to old vale
    set(handles.Slope_Corr_Field, 'Value', old_val, 'String', num2str(old_val) );
end

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function Fixed_Beta_Val_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Fixed_Beta_Val (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in Saturation_Stat_Toggle.
function Saturation_Stat_Toggle_Callback(hObject, eventdata, handles)
% hObject    handle to Saturation_Stat_Toggle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


if handles.Data_Loaded == 0
    % No data loaded so do nothing
    return;
end


% Use to switch between different stats
Types = [{'Lack of Fit'}, {'Model Comparison'}, {'Loop Closure'}];
Current = get(hObject, 'String');


switch Current
    
    case Types{1}
        % Set to Types{2}
        set(hObject, 'String', Types{2});
        handles = Set_Current_Data(handles);
        
    case Types{2}
        % Set to Types{3}
        set(hObject, 'String', Types{3});
        handles = Set_Current_Data(handles);
        
    case Types{3}
        % Set to Types{1}
        set(hObject, 'String', Types{1});
        handles = Set_Current_Data(handles);
        
end

% if strcmpi(Current, Types{1})
%     
%     % Lack of fit is currently displayed - switch to model comparison
%     set(hObject, 'String', Types{2});
%     handles = Set_Current_Data(handles);
%     
% else
%     
%     % swith to lack of fit
%     set(hObject, 'String', Types{1});
%     handles = Set_Current_Data(handles);
%     
% end

guidata(hObject, handles);


% --- Executes on button press in CB_Trim_Fields.
function CB_Trim_Fields_Callback(hObject, eventdata, handles)
% hObject    handle to CB_Trim_Fields (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CB_Trim_Fields

if handles.Data_Loaded == 0
    % No data loaded so do nothing
    return;
end

if get(hObject,'Value') == 1
    set(handles.Trim_Field_Lim, 'Enable', 'on');
else
    set(handles.Trim_Field_Lim, 'Enable', 'off');
end


% Call the ReFit function and update handles structure
handles = ReFit_Data(handles);
guidata(hObject, handles);


function Trim_Field_Lim_Callback(hObject, eventdata, handles)
% hObject    handle to Trim_Field_Lim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Trim_Field_Lim as text
%        str2double(get(hObject,'String')) returns contents of Trim_Field_Lim as a double

% input should be a number
Old_val = get(handles.Trim_Field_Lim, 'Value');

if isnan(str2double(get(hObject,'String') )) || str2double(get(hObject,'String') ) < 0
    set(handles.Trim_Field_Lim, 'Value', Old_val, 'String', num2str(Old_val) );
else
    set(handles.Trim_Field_Lim, 'Value', str2double(get(hObject,'String') ) );
end


if handles.Data_Loaded == 0
    % No data loaded so reset and return
    set(handles.Trim_Field_Lim, 'Value', Old_val, 'String', num2str(Old_val) );
    return;
end

% Check there are enough data to process
Trimmed_Fields = handles.Current_Raw_Loop(abs(handles.Current_Raw_Loop(:,1)) < str2double(get(hObject,'String') ),1);
Max_field = 0.9* max(abs(Trimmed_Fields));

if sum(abs(Trimmed_Fields) > Max_field)/2 < 4
    MSG = ['Trimming these fields will leave too few data at ',...
    'high-fields for adequate high-field slope analysis. Please lower the field.'];
    warndlg(MSG, 'Insufficient data')

    set(handles.Trim_Field_Lim, 'Value', Old_val, 'String', num2str(Old_val) );
else
    % Call the ReFit function and update handles structure
    handles = ReFit_Data(handles);
    guidata(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function Trim_Field_Lim_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Trim_Field_Lim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Drift_Corr_Menu.
function Drift_Corr_Menu_Callback(hObject, eventdata, handles)
% hObject    handle to Drift_Corr_Menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Drift_Corr_Menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Drift_Corr_Menu

if handles.Data_Loaded == 0
    % No data loaded so do nothing
    return;
end


% Call the ReFit function and update handles structure
handles = ReFit_Data(handles);
% Update_Plots(handles);
guidata(hObject, handles);


% --- Executes on button press in CB_Offset_Corr.
function CB_Offset_Corr_Callback(hObject, eventdata, handles)
% hObject    handle to CB_Offset_Corr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CB_Offset_Corr

if handles.Data_Loaded == 0
    % No data loaded so do nothing
    return;
end


% Call the ReFit function and update handles structure
handles = ReFit_Data(handles);
guidata(hObject, handles);


% --- Executes on button press in CB_Pole_Sat.
function CB_Pole_Sat_Callback(hObject, eventdata, handles)
% hObject    handle to CB_Pole_Sat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CB_Pole_Sat

if handles.Data_Loaded == 0
    % No data loaded so do nothing
    return;
end


% Call the ReFit function and update handles structure
handles = ReFit_Data(handles);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function Drift_Corr_Menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Drift_Corr_Menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in CB_Use_Fit_Params.
function CB_Use_Fit_Params_Callback(hObject, eventdata, handles)
% hObject    handle to CB_Use_Fit_Params (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CB_Use_Fit_Params

if handles.Data_Loaded == 0
    % No data loaded so do nothing
    return;
end


% warndlg('Sorry, this feature is not yet implemented', 'Still to do', 'modal')
%
% set(handles.CB_Use_Fit_Params, 'Value', 0);
%
% return;


% Call the ReFit function and update handles structure
handles = ReFit_Data(handles);
guidata(hObject, handles);


function Spec_Mass_Callback(hObject, eventdata, handles)
% hObject    handle to Spec_Mass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Spec_Mass as text
%        str2double(get(hObject,'String')) returns contents of Spec_Mass as a double

% input should be a positive number, but allow letters to specify NaN
Old_val = get(handles.Spec_Mass, 'Value');

if str2double(get(hObject,'String') ) <= 0
    set(handles.Spec_Mass, 'Value', Old_val, 'String', num2str(Old_val) );
else
    set(handles.Spec_Mass, 'Value', str2double(get(hObject,'String') ) );
    handles.All_Masses(handles.spec_ind) = str2double(get(hObject,'String'));
    handles.Current_Spec_Mass = handles.All_Masses(handles.spec_ind);
end


% Enable/disable the mass normalization button
if isnan(get(handles.Spec_Mass, 'Value'))
    set(handles.Spec_Mass, 'Value', NaN, 'String', 'NaN' );
    set(handles.RB_Mass_Norm, 'Enable', 'off');
else
    set(handles.RB_Mass_Norm, 'Enable', 'on');
end

% Update the current data
% Used if the mass normalizer is already selected and the mass updated
handles = Set_Current_Data(handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in Prev_Spec.
function Prev_Next_Spec_Callback(hObject, eventdata, handles, str)
% hObject    handle to Prev_Spec (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

try
    % Get the index pointer and the addresses
    index = handles.spec_ind;
    Nspec = handles.Nspec;
catch
    return;
end

% Depending on whether Prev or Next was clicked change the display
switch str
    case 'Prev'
        % Decrease the index by one
        ind = index - 1;
        % If the index is less then one then set it the number of specimens
        % (Nspec)
        if ind < 1
            ind = Nspec;
        end
    case 'Next'
        % Increase the index by one
        ind = index + 1;
        
        % If the index is greater than the snumber of specimens set index
        % to 1
        if ind > Nspec
            ind = 1;
        end
end

handles.spec_ind=ind;
set(handles.Spec_Num, 'String', handles.All_Names{handles.spec_ind});

guidata(hObject,handles);

func_handles=Set_Current_Data(handles);
handles=func_handles;

guidata(hObject,handles);


% --- Executes when selected object is changed in BG_Normalization.
function BG_Normalization_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in BG_Normalization
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

if handles.Data_Loaded == 0
    % No data loaded so do nothing
    return;
end


switch get(eventdata.NewValue,'Tag') % Get Tag of selected object.
    
    case 'RB_No_Norm'
        handles.Norm_Type = 'None';
        
        %         set(handles.RB_No_Norm, 'Value', 1);
        %         set(handles.RB_Ms_Norm, 'Value', 0);
        %         set(handles.RB_Mass_Norm, 'Value', 0);
    case 'RB_Ms_Norm'
        handles.Norm_Type = 'Ms';
        
        %         set(handles.RB_No_Norm, 'Value', 0);
        %         set(handles.RB_Ms_Norm, 'Value', 1);
        %         set(handles.RB_Mass_Norm, 'Value', 0);
    case 'RB_Mass_Norm'
        handles.Norm_Type = 'Mass';
        %         set(handles.RB_No_Norm, 'Value', 0);
        %         set(handles.RB_Ms_Norm, 'Value', 0);
        %         set(handles.RB_Mass_Norm, 'Value', 1);
        
end

handles = Set_Current_Data(handles);

guidata(hObject, handles);


% --- Executes on button press in CB_Plot_Processed.
function CB_Plot_Processed_Callback(hObject, eventdata, handles)
% hObject    handle to CB_Plot_Processed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CB_Plot_Processed

if handles.Data_Loaded == 0
    % No data loaded so do nothing
    return;
end


Update_Plots(handles)


% --- Executes on button press in CB_Plot_Raw.
function CB_Plot_Raw_Callback(hObject, eventdata, handles)
% hObject    handle to CB_Plot_Raw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CB_Plot_Raw

if handles.Data_Loaded == 0
    % No data loaded so do nothing
    return;
end

Update_Plots(handles)


% --- Executes on button press in CB_Plot_Fitted.
function CB_Plot_Fitted_Callback(hObject, eventdata, handles)
% hObject    handle to CB_Plot_Fitted (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CB_Plot_Fitted
if handles.Data_Loaded == 0
    % No data loaded so do nothing
    return;
end

Update_Plots(handles)


% --- Executes on button press in CB_irh_Gradient.
function CB_irh_Gradient_Callback(hObject, eventdata, handles)
% hObject    handle to CB_irh_Gradient (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CB_irh_Gradient
if handles.Data_Loaded == 0
    % No data loaded so do nothing
    return;
end

Update_Plots(handles)


% --- Executes on button press in CB_Plot_Smooth_Noise.
function CB_Plot_Smooth_Noise_Callback(hObject, eventdata, handles)
% hObject    handle to CB_Plot_Smooth_Noise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CB_Plot_Smooth_Noise
if handles.Data_Loaded == 0
    % No data loaded so do nothing
    return;
end

Update_Plots(handles)


% --- Executes on button press in CB_Plot_Grid.
function CB_Plot_Grid_Callback(hObject, eventdata, handles)
% hObject    handle to CB_Plot_Grid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CB_Plot_Grid

if handles.Data_Loaded == 0
    % No data loaded so do nothing
    return;
end

Update_Plots(handles)


% --- Executes on button press in CB_Origin_Cross.
function CB_Origin_Cross_Callback(hObject, eventdata, handles)
% hObject    handle to CB_Origin_Cross (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CB_Origin_Cross

if handles.Data_Loaded == 0
    % No data loaded so do nothing
    return;
end

Update_Plots(handles)


% --- Executes on button press in DB.
function DB_Callback(hObject, eventdata, handles)
% hObject    handle to DB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
keyboard




% --- Executes during object creation, after setting all properties.
function Stat_Mo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Stat_Mo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function Stat_Qf_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Stat_Qf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function Stat_Q_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Stat_Q (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function Stat_Bo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Stat_Bo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function Stat_Bih_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Stat_Bih (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function Stat_Mrs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Stat_Mrs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function Stat_Ms_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Stat_Ms (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function Stat_Bc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Stat_Bc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on mouse press over axes background.
function Loop_Axes_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to Loop_Axes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isdeployed
    PopOutFigure(handles.Loop_Axes, 'Hysteresis')
end


% --- Executes on mouse press over axes background.
function IR_Axes_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to IR_Axes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isdeployed
    PopOutFigure(handles.IR_Axes, 'Remanence')
end

% --- Executes on mouse press over axes background.
function Noise_Axes_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to Noise_Axes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isdeployed
    PopOutFigure(handles.Noise_Axes, 'Noise')
end


%% The Menu bar functions

% --------------------------------------------------------------------
function MB_Save_Plots_Callback(hObject, eventdata, handles)
% hObject    handle to MB_Save_Plots (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.Data_Loaded == 0
    % No data loaded so do nothing
    return;
end


[file,path] = uiputfile(strcat(handles.All_Names{handles.spec_ind}, '.eps'),'Save the specimen plots...');

if ~ischar(file) && file==0
    % User has cancelled
    % Do nothing and...
    return;
end

tmpFig = Publish_Figure(handles.Loop_Axes, handles.IR_Axes, handles.Noise_Axes);

print(tmpFig, '-depsc', strcat(path, file));
close(tmpFig);


% --------------------------------------------------------------------
function MB_Save_Stats_Callback(hObject, eventdata, handles)
% hObject    handle to MB_Save_Stats (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.Data_Loaded == 0
    % No data loaded so do nothing
    return;
end

[file,path] = uiputfile(strcat('Hysteresis_Stats_&_Params.dat'),'Save the specimen statistics...');

if ~ischar(file) && file==0
    % User has cancelled
    % Do nothing and...
    return;
end

Save_HystLab_Stats(handles.All_Names, handles.All_Masses, handles.All_Data_Parameters, handles.All_Processing_Parameters, file, path);


% --------------------------------------------------------------------
function MB_Export_Loop_Callback(hObject, eventdata, handles)
% hObject    handle to MB_Export_Loop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% TODO - Make more flexible
if handles.Data_Loaded == 0
    % No data loaded so do nothing
    return;
end

Spec_Name = handles.All_Names{handles.spec_ind};


% Start with the loop data 
[file,path] = uiputfile(strcat(Spec_Name, '_Loops.dat'),'Save the specimen loops...');

if ~ischar(file) && file==0
    % User has cancelled
    % Do nothing and...
%     return;
else
    % Not cancelled so continue...
    
    switch handles.Norm_Type
        case 'None'
            % Do nothing (or normalize by 1)
            Normalizer = 1;
            Header = [{'Field [mT]'}, {'Moment [Am^2]'}, {'Fitted Moment [Am^2]'}];
        case 'Ms'
            % Normalize by Ms
            Normalizer = handles.Current_Data_Parameters(5);
            Header = [{'Field [mT]'}, {'Moment/Saturation'}, {'Fitted Moment/Saturation'}];
        case 'Mass'
            Normalizer = handles.Current_Spec_Mass .* 1e-6; % Convert to kg
            Header = [{'Field [mT]'}, {'Magnetization [Am^2/kg]'}, {'Fitted Magnetization [Am^2/kg]'}];
    end
    
    
    Fields = [handles.Current_Processed_Data(:,1); handles.Current_Processed_Data(:,2)];
    Mproc = [handles.Current_Processed_Data(:,3); handles.Current_Processed_Data(:,4)] ./ Normalizer;
    Mfit = [handles.Current_Fitted_Data(:,3); handles.Current_Fitted_Data(:,4)] ./ Normalizer;
    
    Data = [num2cell(Fields), num2cell(Mproc), num2cell(Mfit)]';
    
    FID = fopen(strcat(path, file), 'wt');
    
    fprintf(FID, '%s\t%s\t%s\n', Header{:});
    fprintf(FID, '%4.2f\t%1.4e\t%1.4e\n', Data{:});
    
    fclose(FID);
    
end


% Start with the loop data 
[file,path] = uiputfile(strcat(Spec_Name, '_Curves.dat'),'Save the specimen Mrh, Mih curves...');

if ~ischar(file) && file==0
    % User has cancelled
    % Do nothing and...
    return;
else
    % Not cancelled so continue...
    
    switch handles.Norm_Type
        case 'None'
            % Do nothing (or normalize by 1)
            Normalizer = 1;
            Header = [{'Field [mT]'}, {'Mrh [Am^2]'}, {'Mih [Am^2]'}, {'Noise [Am^2]'},...
                {'Fitted Mrh [Am^2]'}, {'Fitted Mih [Am^2]'}, {'Smoothed Noise [Am^2]'}];
        case 'Ms'
            % Normalize by Ms
            Normalizer = handles.Current_Data_Parameters(5);
            Header = [{'Field [mT]'}, {'Mrh/Saturation'}, {'Mih/Saturation'}, {'Noise/Saturation'},...
                {'Fitted Mrh/Saturation'}, {'Fitted Mih/Saturation'}, {'Smoothed Noise/Saturation'}];
        case 'Mass'
            Normalizer = handles.Current_Spec_Mass .* 1e-6; % Convert to kg
            Header = [{'Field [mT]'}, {'Mrh [Am^2/kg]'}, {'Mih [Am^2/kg]'}, {'Noise [Am^2/kg]'},...
                 {'Fitted Mrh [Am^2/kg]'}, {'Fitted Mih [Am^2/kg]'}, {'Smoothed Noise [Am^2/kg]'}];
    end
    
    
    Fields = handles.Current_Processed_Data(:,1);
    Mproc = [handles.Current_Processed_Data(:,6), handles.Current_Processed_Data(:,5), handles.Current_Noise_Data(:,2)] ./ Normalizer;
    Mfit = [handles.Current_Fitted_Data(:,6), handles.Current_Fitted_Data(:,5), handles.Current_Noise_Data(:,3)] ./ Normalizer;
    
    Data = [num2cell(Fields), num2cell(Mproc), num2cell(Mfit)]';
    
    FID = fopen(strcat(path, file), 'wt');
    
    fprintf(FID, '%s\t%s\t%s\t%s\t%s\t%s\t%s\n', Header{:});
    fprintf(FID, '%4.2f\t%1.4e\t%1.4e\t%1.4e\t%1.4e\t%1.4e\t%1.4e\n', Data{:});
    
    fclose(FID);
    
end

% --------------------------------------------------------------------
function MB_MagIC_Callback(hObject, eventdata, handles)
% hObject    handle to MB_MagIC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.Data_Loaded == 0
    % No data loaded so do nothing
    return;
end


% --------------------------------------------------------------------
function MB_Add_Data_Callback(hObject, eventdata, handles)
% hObject    handle to MB_Add_Data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% If data are loaded, check if the user wants to overwrite, append or cancel
if handles.Data_Loaded == 1
    
    Choice = questdlg('Do you wish to append to the current data set or overwrite it?',...
        'Load Data', 'Append', 'Overwrite', 'Cancel', 'Cancel');
    
    
    switch Choice
        case 'Append'
            OW_Flag = 0;
        case 'Overwrite'
            OW_Flag = 1;
        otherwise
            % Assume cancel and return
            return;
    end
else
    % No data are loaded so overwrite
    OW_Flag = 1;
end


[file, path]=uigetfile({'*', 'All Files (*.*)';...
    '*.hys;*.hyst;*.HYS;*.HYST', 'Hys Files (*.hys, *.hyst)'; '*.loop;*.LOOP', 'Loop Files (*.loop)'},...
    'Please select your data file(s).', 'MultiSelect', 'on');


if iscell(file) == 0 & file == 0  %#ok<AND2>
    % User cancelled the open dialog box - Do nothing
    % But old data still exist
    %     handles.Data_Loaded = 0;
    guidata(hObject, handles);
    
elseif size(file,2) ~=0 % load and process
    
    % Create a cell if only a single file is selected
    if ~iscell(file)
        file = {file};
    end
    
    % Load the options window to get details on the machine type
    Options = HystLab_Load_Options('Main_Window_Call', handles.HystLab_Fig);
    
    % Check if the user canceled
    if Options.LoadStatus == 0
        % Canceled, so just return
        return;
    end
    
    % Read the data
    try
        [Specimen_Names, Data, Masses, Data_Order] = Read_Hyst_Files(path, file, Options.File_Type);
    catch err
        %         MSG = [{'Data files could not be loaded.'}; {'Please check the file formats.'};...
        %             {'If the format is currently unsupported, please contact Greig Paterson.'}];
        % TODO - Add some additional error information for the various possibilities
        MSG = [{'Data files could not be loaded.'}; {err.message}];
        warndlg(MSG, 'Load Data Failure');
        return;
    end
    
    % Warning of specimens with missing data
    if any(cellfun(@isempty, Data))
        Error_MSG = [{'The following specimens have no readable data:'};...
            Specimen_Names(cellfun(@isempty, Data));...
            {'These specimens cannot be displayed. Please check the data files.'}];
        waitfor(warndlg(Error_MSG, 'Reading Failed', 'modal'));
    end
    
    % Remove from the read data
    Specimen_Names(cellfun(@isempty, Data)) = [];
    Masses(cellfun(@isempty, Data)) = [];
    Data_Order(cellfun(@isempty, Data)) = [];
    Data(cellfun(@isempty, Data)) = [];
    
    if isempty(Data)
        % All data removed so send a cancel processing flag to return
        handles.Cancel_Processing = 1;
    end
    
    
    % Check for cancel flag and return exiting data if needed
    if handles.Cancel_Processing == 1
        % User canceled the processing
        % Reset the cancel processing dlag
        handles.Cancel_Processing = 0;
        guidata(hObject, handles);
        return;
    end
 
    
    % Create/append the data
    if OW_Flag == 1
        % Overwrite existing data
        
        tmp_handles = handles;
        tmp_handles.Data_Loaded = 0; % Reset for the call to Set_Initial
        tmp_handles.All_Names = Specimen_Names;
        tmp_handles.All_Data = Data;
        tmp_handles.All_Masses = Masses;
        tmp_handles.All_Data_Order = Data_Order;
        tmp_handles.All_Processing_Parameters = [];%DataTransfer.All_Processing_Parameters;
        tmp_handles.Nspec = length(tmp_handles.All_Names);
        tmp_handles = Set_Initial_Data(tmp_handles, 0);
        

        
        %         handles.Data_Loaded = 0; % Reset for the call to Set_Initial
        %         handles.All_Names = Specimen_Names;
        %         handles.All_Data = Data;
        %         handles.All_Masses = Masses;
        %         handles.All_Processing_Parameters = [];%DataTransfer.All_Processing_Parameters;
        %         handles.Nspec = length(handles.All_Names);
        
        handles = tmp_handles;
        
        
    else
        % Append new data to existing data
       
        % Put into temporary handle and get the default processing
        tmp.Data_Loaded = 1; % Reset for the call to Set_Initial
        tmp.All_Names = Specimen_Names;
        tmp.All_Data = Data;
        tmp.All_Masses = Masses;
        tmp.All_Data_Order = Data_Order;
        tmp.All_Processing_Parameters = [];
        tmp.Nspec = length(tmp.All_Names);
        tmp.spec_ind = handles.spec_ind;
        tmp = Set_Initial_Data(tmp, 0);
        
        handles.Cancel_Processing = 0;
        
        % Append the temporary handles to the main handles
        handles.All_Names = [handles.All_Names; tmp.All_Names];
        handles.All_Data = [handles.All_Data; tmp.All_Data];
        handles.All_Masses = [handles.All_Masses; tmp.All_Masses];
        handles.All_Data_Order = [handles.All_Data_Order; tmp.All_Data_Order];

        handles.All_Processed_Data = [handles.All_Processed_Data; tmp.All_Processed_Data];
        handles.All_Uncorrected_Data = [handles.All_Uncorrected_Data; tmp.All_Uncorrected_Data];
        handles.All_Noise_Data = [handles.All_Noise_Data; tmp.All_Noise_Data];
        handles.All_Fitted_Data = [handles.All_Fitted_Data; tmp.All_Fitted_Data];
        handles.All_Data_Parameters = [handles.All_Data_Parameters; tmp.All_Data_Parameters];
        handles.All_Processing_Parameters = [handles.All_Processing_Parameters; tmp.All_Processing_Parameters];
        handles.All_Basis_Coeffs = [handles.All_Basis_Coeffs; tmp.All_Basis_Coeffs];
        
        % update the number of specimens
        handles.Nspec = length(handles.All_Names);
        
    end
    
    
    if handles.Cancel_Processing ~= 1
        
        % Set data loaded status
        handles.Data_Loaded = 1;
        
        % save the updated handles
        guidata(hObject, handles);
    end
    
    
end


% --------------------------------------------------------------------
function MB_Remove_Data_Callback(hObject, eventdata, handles)
% hObject    handle to MB_Remove_Data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


if handles.Data_Loaded == 0
    % No data so return
    return;
end

% Only one specimn so don't remove, return
if handles.Nspec == 1
    MSG = {'Only 1 specimn is currently loaded.';...
        'Please load new data or quite the browser'};
    
    warndlg(MSG, '1 specimen loaded', 'modal')
    return;
end


Title = ['Remove ', char(handles.All_Names{handles.spec_ind}), '?'];

Choice = questdlg('Do you wish to remove the current specimen?',...
    Title, 'Yes', 'No', 'No');


% User didn't choose "Yes" so return
if ~strcmpi(Choice, 'Yes')
    return;
end

% Passed all the checks, so remove the specimen...

% Get the specimen index
spec_ind = handles.spec_ind;

% Remove the specimen from all relevant handles
handles.All_Names(spec_ind) = [];
handles.All_Data(spec_ind) = [];
handles.All_Masses(spec_ind) = [];

handles.All_Processed_Data(spec_ind) = [];
handles.All_Uncorrected_Data(spec_ind) = [];
handles.All_Noise_Data(spec_ind) = [];
handles.All_Fitted_Data(spec_ind) = [];
handles.All_Data_Parameters(spec_ind,:) = [];
handles.All_Processing_Parameters(spec_ind,:) = [];
handles.All_Basis_Coeffs(spec_ind,:) = [];

% Update the number of specimens
handles.Nspec = length(handles.All_Names);


% reset the current specimen
if spec_ind-1 < 1
    handles.spec_ind = 1;
else
    handles.spec_ind = spec_ind - 1;
end

if handles.Nspec == 0
    % No data left
    % Reset data loaded status
    handles.Data_Loaded = 0;
    
    % Reset plots etc
    
else
    % Still have data
    % Update plots etc
    
    handles = Set_Current_Data(handles);
    set(handles.Spec_Num, 'String', handles.All_Names{handles.spec_ind});
    
end

% save the updated handles
guidata(hObject, handles);


% --------------------------------------------------------------------
function MB_Save_Session_Callback(hObject, eventdata, handles)
% hObject    handle to MB_Save_Session (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.Data_Loaded == 0
    % No data loaded so do nothing
    return;
end


Save_HystLab_Session(handles, 'Hys');


% --------------------------------------------------------------------
function MB_Load_Session_Callback(hObject, eventdata, handles)
% hObject    handle to MB_Load_Session (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Load the saved session

[file,path] = uigetfile('*.mat','Load a saved session...');

if ~ischar(file) && file==0
    % User has cancelled
    % Do nothing and...
    return;
end

% try loading the file
load(strcat(path, file), 'Session_handles', 'Type');

% Check to see if load threw any warning about missing variables
if strcmpi(lastwarn, 'Variable ''Type'' not found.') || strcmpi(lastwarn, 'Variable ''Session_handles'' not found.')
    % Warn about invalid mat files and return
    warndlg('This is not a valid session file. Please try another.', 'Invalid files');
    lastwarn(''); % clear that last warning
    return;
end

% Catch non-hysteresis session files
if ~strcmpi(Type, 'Hys')
    % Warn about non-hysteresis session file and return
    warndlg('This is not a hysteresis session file. Please try another.', 'Non-Hysteresis Session');
    return;
end


tmp_handles = Load_HystLab_Session(handles, Session_handles);

if isempty(tmp_handles) || tmp_handles.Session_Fail == 1
        % An error has ocurred in the loading so simply return
        warndlg('The session failed to load. Please see the MATLAB command window for details.', 'Cannot load session');
    return;
end

% Save the loaded version for debugging etc
tmp_handles.HystLab_Version_Loaded = tmp_handles.HystLab_Version;

% Update the HystLab version to the current version
tmp_handles.HystLab_Version = handles.HystLab_Version;

handles = tmp_handles;

% Reset the normalization
handles.Norm_Type = 'None';
set(handles.RB_No_Norm, 'Value', 1)

% Reset the Saturation stats display
set(handles.Saturation_Stat_Toggle, 'String', 'Lack of Fit');

% save the updated handles
guidata(hObject, handles);

% Update the data
handles = Set_Current_Data(handles);

% Set the the specimen browser label
set(handles.Spec_Num, 'String', handles.All_Names{handles.spec_ind});

if handles.Data_Loaded == 0
    handles.Data_Loaded = 1;
end

% save the updated handles
guidata(hObject, handles);


% --------------------------------------------------------------------
function MB_Pole_Sat_Callback(hObject, eventdata, handles)
% hObject    handle to MB_Pole_Sat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function MB_Redo_All_Callback(hObject, eventdata, handles)
% hObject    handle to MB_Redo_All (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Function to re-process all loaded data using the set data parameters
% This is useful if the processing functions have been updated
% Required data to be loaded.

if handles.Data_Loaded == 0
    % No data so return
    return;
end

handles = Set_Initial_Data(handles, 1);
guidata(hObject, handles);


% --------------------------------------------------------------------
function MB_Set_Colors_Callback(hObject, eventdata, handles)
% hObject    handle to MB_Set_Colors (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Return_Data = Set_HystLab_Colors(handles.HystLab_Fig);

if Return_Data.CancelFlag == 1
    return;
end

% Update the defaults
Defaults = Return_Data.Defaults;
setappdata(handles.HystLab_Fig, 'Defaults', Defaults);

handles.Default_Hyst_Color = Defaults.Hyst_Plot_Color;
handles.Default_IR_Color = Defaults.IR_Plot_Color;
handles.Default_Noise_Color = Defaults.Noise_Plot_Color;


guidata(hObject,handles);

% Update the plots
Update_Plots(handles)


% --------------------------------------------------------------------
function MB_Set_Default_Proc_Callback(hObject, eventdata, handles)
% hObject    handle to MB_Set_Default_Proc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function MB_Open_Manual_Callback(hObject, eventdata, handles)
% hObject    handle to MB_Open_Manual (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get the current path of the main m-file
S = mfilename('fullpath');
name_len = length(mfilename());
MyPath = S(1:end-name_len);

SPLT = strsplit(handles.HystLab_Version, '.');

% file_name = strcat('HystLab_Manual_v', strjoin([SPLT(1:2), {'x'}], '.'), '.pdf');
file_name = strcat('HystLab_Manual_v', strjoin(SPLT(1:3), '.'), '.pdf');
file_path = strcat(MyPath, './Documents/', file_name);

try
    open(file_path);
catch
    warndlg([file_name, ' not found.'], 'Manual Not Found');
end


function MB_Open_Paper_Callback(hObject, eventdata, handles)
% hObject    handle to MB_Open_Manual (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get the current path of the main m-file
S = mfilename('fullpath');
name_len = length(mfilename());
MyPath = S(1:end-name_len);

file_name = 'Paterson et al., 2018, Measuring, Processing, and Analyzing Hysteresis Data.pdf';
file_path = strcat(MyPath, './Documents/', file_name);

try
    open(file_path);
catch
    warndlg([file_name, ' not found.'], 'Paper Not Found');
end


% --------------------------------------------------------------------
function MB_About_Callback(hObject, eventdata, handles)
% hObject    handle to MB_About (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

About_HystLab('Version', handles.HystLab_Version, 'Date', handles.HystLab_Date);


% --------------------------------------------------------------------
function MB_BiPLot_Callback(hObject, eventdata, handles)
% hObject    handle to MB_BiPLot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


try
    Transfer.Data = handles.All_Data_Parameters;
    Transfer.Names = handles.All_Names;
catch
    warndlg('No data currently loaded.', 'No Data', 'modal')
    return;
end

Hysteresis_BiPlot('DataTransfer', Transfer, handles.HystLab_Fig);
