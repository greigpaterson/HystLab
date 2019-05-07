function varargout = HystLab_Load_Options(varargin)
% HYSTLAB_LOAD_OPTIONS MATLAB code for HystLab_Load_Options.fig
%      HYSTLAB_LOAD_OPTIONS, by itself, creates a new HYSTLAB_LOAD_OPTIONS or raises the existing
%      singleton*.
%
%      H = HYSTLAB_LOAD_OPTIONS returns the handle to a new HYSTLAB_LOAD_OPTIONS or the handle to
%      the existing singleton*.
%
%      HYSTLAB_LOAD_OPTIONS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in HYSTLAB_LOAD_OPTIONS.M with the given input arguments.
%
%      HYSTLAB_LOAD_OPTIONS('Property','Value',...) creates a new HYSTLAB_LOAD_OPTIONS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before HystLab_Load_Options_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to HystLab_Load_Options_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help HystLab_Load_Options

% Last Modified 2019/05/07
%

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @HystLab_Load_Options_OpeningFcn, ...
    'gui_OutputFcn',  @HystLab_Load_Options_OutputFcn, ...
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


% --- Executes just before HystLab_Load_Options is made visible.
function HystLab_Load_Options_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to HystLab_Load_Options (see VARARGIN)


% Choose default command line output for HYSTLAB_LOAD_OPTIONS
handles.output = hObject;

if isempty(varargin)
    error('HystLab_Load_Options:Opening', 'Standalone calls are not supported.')
end

dontOpen = false;
mainGuiInput = find(strcmp(varargin, 'Main_Window_Call'));

if (isempty(mainGuiInput)) ...
        || (length(varargin) <= mainGuiInput) ...
        || (~ishandle(varargin{mainGuiInput+1}))
    
    dontOpen = true;
    
else
    % Set the title
    set(hObject, 'Name', 'Data file load options');
    % Remember the handle, and adjust our position
    handles.MainWindow = varargin{mainGuiInput+1};
    
    % Position to be relative to parent:
    parentPosition = get(handles.MainWindow, 'Position'); %getpixelposition(handles.MainWindow)
    currentPosition = get(hObject, 'Position');
    
    newX = parentPosition(1) + (parentPosition(3)/2 - currentPosition(3)/2);
    newY = parentPosition(2) + (parentPosition(4)/2 - currentPosition(4)/2);
    newW = currentPosition(3);
    newH = currentPosition(4);
    
    set(hObject, 'Position', [newX, newY, newW, newH]);
    
end


% Set the default values to hysteresis data and MicroMag
handles.File_Type_Flag = 1;


% Update handles structure
guidata(hObject, handles);

if dontOpen
    disp('------------------------------------------');
    disp('If you see this, something has gone wrong.')
    disp('      Please contact Greig Paterson       ')
    disp('------------------------------------------');
else
    % Set the uiwait
    uiwait(handles.Load_Options_Figure);
end



% --- Outputs from this function are returned to the command line.
function varargout = HystLab_Load_Options_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

delete(hObject);


% --- Executes when selected object is changed in File_Format_Panel.
function File_Format_Panel_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in File_Format_Panel
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

switch get(eventdata.NewValue,'Tag') % Get Tag of selected object.
    case 'MicroMag_File_Input'
        handles.File_Type_Flag = 1;
        set(handles.Other_File_Menu, 'Enable', 'off');
    case 'MPMS_File_Input'
        handles.File_Type_Flag = 2;
        set(handles.Other_File_Menu, 'Enable', 'off');
    case 'VFTB_File_Input'
        handles.File_Type_Flag = 3;
        set(handles.Other_File_Menu, 'Enable', 'off');
    case 'LakeShore_File_Input'
        handles.File_Type_Flag = 4;
        set(handles.Other_File_Menu, 'Enable', 'off');
%     case 'MagIC_File_Input'
%         handles.File_Type_Flag = 5;
%         set(handles.Other_File_Menu, 'Enable', 'off');
    case 'Other_File_Input'
        set(handles.Other_File_Menu, 'Enable', 'on');
        handles.File_Type_Flag = 'Other';
end

guidata(hObject, handles);


% --- Executes on button press in OK_Button.
function OK_Button_Callback(hObject, eventdata, handles)
% hObject    handle to OK_Button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Set the fit status to zero and return
Return_struct.LoadStatus = 1;

if strcmpi(handles.File_Type_Flag, 'Other')
    % Get the file format value from the pull down menu
    str = get(handles.Other_File_Menu, 'String');
    val = get(handles.Other_File_Menu, 'Value');
    
    if iscell(str)
        str = str{val};
    end
    
    switch str % Get Tag of selected object.
        case 'MicroSense VSM'
            handles.File_Type_Flag = 6;
        case 'Generic 2 Column'
            handles.File_Type_Flag = 7;
        case 'MolSpin VSM'
            handles.File_Type_Flag = 8;
        otherwise
            disp('Still to add');
            return;
    end
end

% Get the options
Return_struct.File_Type = handles.File_Type_Flag;

handles.output = Return_struct;
guidata(hObject,handles);

% Return to the main GUI window
uiresume(handles.Load_Options_Figure);


% --- Executes on button press in Cancel_Button.
function Cancel_Button_Callback(hObject, eventdata, handles)
% hObject    handle to Cancel_Button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Set the fit status to zero and return
Return_struct.LoadStatus = 0;

handles.output = Return_struct;
guidata(hObject,handles);

% Return to the main GUI window
uiresume(handles.Load_Options_Figure);


% --- Executes when user attempts to close Load_Options_Figure.
function Load_Options_Figure_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to Load_Options_Figure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Set the fit status to zero and return
Return_struct.LoadStatus = 0;

handles.output = Return_struct;
guidata(hObject,handles);

% Return to the main GUI window
uiresume(handles.Load_Options_Figure);


% --- Executes on button press in DB.
function DB_Callback(hObject, eventdata, handles)
% hObject    handle to DB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
keyboard


% --- Executes during object creation, after setting all properties.
function Other_File_Menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Other_File_Menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Other_File_Menu.
function Other_File_Menu_Callback(hObject, eventdata, handles)
% hObject    handle to Other_File_Menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Other_File_Menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Other_File_Menu
