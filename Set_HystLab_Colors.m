function varargout = Set_HystLab_Colors(varargin)
% Set_HystLab_Colors MATLAB code for Set_HystLab_Colors.fig
%      Set_HystLab_Colors, by itself, creates a new Set_HystLab_Colors or raises the existing
%      singleton*.
%
%      H = Set_HystLab_Colors returns the handle to a new Set_HystLab_Colors or the handle to
%      the existing singleton*.
%
%      Set_HystLab_Colors('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in Set_HystLab_Colors.M with the given input arguments.
%
%      Set_HystLab_Colors('Property','Value',...) creates a new Set_HystLab_Colors or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Set_HystLab_Colors_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Set_HystLab_Colors_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Set_HystLab_Colors

% Last Modified 2019/05/07
%

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @Set_HystLab_Colors_OpeningFcn, ...
    'gui_OutputFcn',  @Set_HystLab_Colors_OutputFcn, ...
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


% --- Executes just before Set_HystLab_Colors is made visible.
function Set_HystLab_Colors_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Set_HystLab_Colors (see VARARGIN)

% Choose default command line output for Set_HystLab_Colors
handles.output = hObject;


% DataTransfer = find(strcmp(varargin, 'DataTransfer'));
if nargin < 2
    
    disp('------------------------------------------');
    disp('If you see this, something has gone wrong.')
    disp('      Please contact Greig Paterson       ')
    disp('------------------------------------------');
    
else
    handles.MainWindow = varargin{1};
    handles.Defaults = getappdata(varargin{1}, 'Defaults');
end

% keyboard
% Position to be relative to parent:
parentPosition = get(handles.MainWindow, 'Position');
currentPosition = get(hObject, 'Position');

newX = parentPosition(1) + (parentPosition(3)/1.2 - currentPosition(3)/2);
newY = parentPosition(2) + (parentPosition(4)/1.25 - currentPosition(4)/2);
newW = currentPosition(3);
newH = currentPosition(4);
set(hObject, 'Position', [newX, newY, newW, newH]);

try
    S = mfilename('fullpath');
    name_len = length(mfilename());
    MyPath = S(1:end-name_len);
    load(strcat(MyPath,'HystLab_Color_Plot_Data.mat'), 'Example_Data');
    handles.Data = Example_Data;
catch
    error('Set_HystLab_Colors:Data', 'Example data not found.');
end


% Update handles structure
guidata(hObject, handles);

% Make the plot
Update_Plot(handles)

% UIWAIT makes Set_HystLab_Colors wait for user response (see UIRESUME)
uiwait(handles.Set_HystLab_Colors_Fig);


% --- Outputs from this function are returned to the command line.
function varargout = Set_HystLab_Colors_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

delete(hObject);


% --- Executes on button press in OK.
function OK_Callback(hObject, eventdata, handles)
% hObject    handle to OK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Create a return structure
Return.CancelFlag = 0;
Return.Defaults = handles.Defaults;

handles.output = Return;
guidata(hObject,handles);

uiresume(handles.Set_HystLab_Colors_Fig);


% --- Executes on button press in Cancel.
function Cancel_Callback(hObject, eventdata, handles)
% hObject    handle to Cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Return.CancelFlag = 1;
handles.output = Return;
guidata(hObject,handles);

uiresume(handles.Set_HystLab_Colors_Fig);


% --- Executes when user attempts to close Set_HystLab_Colors_Fig.
function Set_HystLab_Colors_Fig_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to Set_HystLab_Colors_Fig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
% delete(hObject);

Return.CancelFlag = 1;
handles.output = Return;
guidata(hObject,handles);

uiresume(handles.Set_HystLab_Colors_Fig);


function Update_Plot(handles)

% PLot params
% width = 0.1;
LineWidth = 1;

Defaults = handles.Defaults;

% Get the colors
Hyst_Color = Defaults.Hyst_Plot_Color;
IR_Color = Defaults.IR_Plot_Color;
Noise_Color = Defaults.Noise_Plot_Color;


% Get the Data
Data = handles.Data;


hold(handles.Plot_Axes_1, 'on')
plot(handles.Plot_Axes_1, Data{1}(:,1)-15, Data{1}(:,2), 'Color', Hyst_Color(1,:), 'LineWidth', LineWidth);
plot(handles.Plot_Axes_1, Data{2}(:,1), 0.9.*Data{2}(:,2), 'Color', Hyst_Color(2,:), 'LineWidth', LineWidth);
plot(handles.Plot_Axes_1, Data{3}(:,1)+15, Data{3}(:,2), 'Color', Hyst_Color(3,:), 'LineWidth', LineWidth);
hold(handles.Plot_Axes_1, 'off')

set(handles.Plot_Axes_1, 'Xlim', [-10, 1e3], 'Ylim', [0 2])

hold(handles.Plot_Axes_2, 'on')
plot(handles.Plot_Axes_2, Data{4}(:,1)-15, Data{4}(:,2), 'Color', IR_Color(1,:), 'LineWidth', LineWidth);
plot(handles.Plot_Axes_2, Data{5}(:,1), 0.95.*Data{5}(:,2), 'Color', IR_Color(2,:), 'LineWidth', LineWidth);
plot(handles.Plot_Axes_2, Data{6}(:,1)+15, 1.05.*Data{6}(:,2), 'Color', IR_Color(3,:), 'LineWidth', LineWidth);

plot(handles.Plot_Axes_2, Data{7}(:,1)-15, 0.8.*Data{7}(:,2), 'Color', IR_Color(1,:), 'LineWidth', LineWidth);
plot(handles.Plot_Axes_2, Data{8}(:,1), Data{8}(:,2), 'Color', IR_Color(2,:), 'LineWidth', LineWidth);
plot(handles.Plot_Axes_2, Data{9}(:,1)+25, 1.1.*Data{9}(:,2), 'Color', IR_Color(3,:), 'LineWidth', LineWidth);
hold(handles.Plot_Axes_2, 'off')

set(handles.Plot_Axes_2, 'Xlim', [-10, 1e3], 'Ylim', [0 2])

hold(handles.Plot_Axes_3, 'on')
plot(handles.Plot_Axes_3, Data{10}(:,1), Data{10}(:,2) + 1.2e-6.*Data{10}(:,1), 'Color', Noise_Color(1,:), 'LineWidth', LineWidth);
plot(handles.Plot_Axes_3, Data{10}(:,1), 0.6.*Data{10}(:,2), 'Color', Noise_Color(2,:), 'LineWidth', LineWidth);
plot(handles.Plot_Axes_3, Data{11}(:,1), Data{11}(:,2), 'Color', Noise_Color(3,:), 'LineWidth', LineWidth);
hold(handles.Plot_Axes_3, 'off')

set(handles.Plot_Axes_3, 'Xlim', [-1e3, 1e3], 'Ylim', [-2e-3 2e-3])



% --- Executes on button press in DB_ME.
function DB_ME_Callback(hObject, eventdata, handles)
% hObject    handle to DB_ME (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
keyboard



% --- Executes on button press in Set_Default.
function Set_Default_Callback(hObject, eventdata, handles)
% hObject    handle to Set_Default (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get the defaults, which are already saved in the function calls
Defaults = handles.Defaults;

% Save them to the user defaul config file
Save_HystLab_User_Defaults(Defaults);

% Save them to the main window appdata
setappdata(handles.MainWindow, 'Defaults', Defaults);





% --- Executes on button press in Set_Hys_Raw_Color.
function Set_Hys_Raw_Color_Callback(hObject, eventdata, handles)
% hObject    handle to Set_Hys_Raw_Color (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get the new color
handles.Defaults.Hyst_Plot_Color(1,:) = uisetcolor(handles.Defaults.Hyst_Plot_Color(1,:));
guidata(hObject, handles);

% Update the plot
Update_Plot(handles)



% --- Executes on button press in Set_Hys_Processed_Color.
function Set_Hys_Processed_Color_Callback(hObject, eventdata, handles)
% hObject    handle to Set_Hys_Processed_Color (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get the new color
handles.Defaults.Hyst_Plot_Color(2,:) = uisetcolor(handles.Defaults.Hyst_Plot_Color(2,:));
guidata(hObject, handles);

% Update the plot
Update_Plot(handles)



% --- Executes on button press in Set_Hys_Fit_Color.
function Set_Hys_Fit_Color_Callback(hObject, eventdata, handles)
% hObject    handle to Set_Hys_Fit_Color (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get the new color
handles.Defaults.Hyst_Plot_Color(3,:) = uisetcolor(handles.Defaults.Hyst_Plot_Color(3,:));
guidata(hObject, handles);

% Update the plot
Update_Plot(handles)



% --- Executes on button press in Set_IR_Raw_Color.
function Set_IR_Raw_Color_Callback(hObject, eventdata, handles)
% hObject    handle to Set_IR_Raw_Color (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get the new color
handles.Defaults.IR_Plot_Color(1,:) = uisetcolor(handles.Defaults.IR_Plot_Color(1,:));
guidata(hObject, handles);

% Update the plot
Update_Plot(handles)


% --- Executes on button press in Set_IR_Processed_Color.
function Set_IR_Processed_Color_Callback(hObject, eventdata, handles)
% hObject    handle to Set_IR_Processed_Color (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get the new color
handles.Defaults.IR_Plot_Color(2,:) = uisetcolor(handles.Defaults.IR_Plot_Color(2,:));
guidata(hObject, handles);

% Update the plot
Update_Plot(handles)


% --- Executes on button press in Set_IR_Fit_Color.
function Set_IR_Fit_Color_Callback(hObject, eventdata, handles)
% hObject    handle to Set_IR_Fit_Color (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get the new color
handles.Defaults.IR_Plot_Color(3,:) = uisetcolor(handles.Defaults.IR_Plot_Color(3,:));
guidata(hObject, handles);

% Update the plot
Update_Plot(handles)


% --- Executes on button press in Set_Noise_Raw_Color.
function Set_Noise_Raw_Color_Callback(hObject, eventdata, handles)
% hObject    handle to Set_Noise_Raw_Color (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get the new color
handles.Defaults.Noise_Plot_Color(1,:) = uisetcolor(handles.Defaults.Noise_Plot_Color(1,:));
guidata(hObject, handles);

% Update the plot
Update_Plot(handles)


% --- Executes on button press in Set_Noise_Color.
function Set_Noise_Color_Callback(hObject, eventdata, handles)
% hObject    handle to Set_Noise_Color (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get the new color
handles.Defaults.Noise_Plot_Color(2,:) = uisetcolor(handles.Defaults.Noise_Plot_Color(2,:));
guidata(hObject, handles);

% Update the plot
Update_Plot(handles)


% --- Executes on button press in Set_Smoothed_Color.
function Set_Smoothed_Color_Callback(hObject, eventdata, handles)
% hObject    handle to Set_Smoothed_Color (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get the new color
handles.Defaults.Noise_Plot_Color(3,:) = uisetcolor(handles.Defaults.Noise_Plot_Color(3,:));
guidata(hObject, handles);

% Update the plot
Update_Plot(handles)




