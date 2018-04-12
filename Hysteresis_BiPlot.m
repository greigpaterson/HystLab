function varargout = Hysteresis_BiPlot(varargin)
% HYSTERESIS_BIPLOT MATLAB code for Hysteresis_BiPlot.fig
%      HYSTERESIS_BIPLOT, by itself, creates a new HYSTERESIS_BIPLOT or raises the existing
%      singleton*.
%
%      H = HYSTERESIS_BIPLOT returns the handle to a new HYSTERESIS_BIPLOT or the handle to
%      the existing singleton*.
%
%      HYSTERESIS_BIPLOT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in HYSTERESIS_BIPLOT.M with the given input arguments.
%
%      HYSTERESIS_BIPLOT('Property','Value',...) creates a new HYSTERESIS_BIPLOT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Hysteresis_BiPlot_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Hysteresis_BiPlot_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Hysteresis_BiPlot

% Last Modified by GUIDE v2.5 06-Sep-2017 17:02:33

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @Hysteresis_BiPlot_OpeningFcn, ...
    'gui_OutputFcn',  @Hysteresis_BiPlot_OutputFcn, ...
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


% --- Executes just before Hysteresis_BiPlot is made visible.
function Hysteresis_BiPlot_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Hysteresis_BiPlot (see VARARGIN)

% Choose default command line output for Hysteresis_BiPlot
handles.output = hObject;

DataTransfer = find(strcmp(varargin, 'DataTransfer'));
if (isempty(DataTransfer)) ...
        || (length(varargin) <= DataTransfer)% ...
    
    disp('------------------------------------------');
    disp('If you see this, something has gone wrong.')
    disp('      Please contact Greig Paterson       ')
    disp('------------------------------------------');
    
else
    DataTransfer = varargin{DataTransfer+1};
    handles.MainWindow = varargin{3};
    handles.Defaults = getappdata(varargin{3}, 'Defaults');
end

% Position to be relative to parent:
parentPosition = get(handles.MainWindow, 'Position'); %getpixelposition(handles.MainWindow)
currentPosition = get(hObject, 'Position');
% Set x and y to be directly in the middle
newX = parentPosition(1) + (parentPosition(3)/2 - currentPosition(3)/2);
newY = parentPosition(2) + (parentPosition(4)/2 - currentPosition(4)/2);
newW = currentPosition(3);
newH = currentPosition(4);
set(hObject, 'Position', [newX, newY, newW, newH]);

% Get the version
Ver = ver('MATLAB');
handles.Version = str2double(Ver.Version);

% Process the defaults
handles.Plot_Symbol = handles.Defaults.BiPlotSymbol;
handles.Symbol_Color = handles.Defaults.BiPlot_Color;
handles.Symbol_Size = handles.Defaults.BiPlotSymbolSize;

if strcmpi(handles.Defaults.BiPlotFaceColor, 'filled')
    handles.Face_Color = handles.Defaults.BiPlot_Color;
else
    handles.Face_Color = 'none';
end


% Get the data from previous window
handles.Data = DataTransfer.Data;
handles.All_Names = DataTransfer.Names;
handles.nData = size(handles.Data, 1);

% 5 - Ms
% 6 - Mrs
% 7 - Bc
% 8 - Brh
handles.Ms = handles.Data(:,5);
handles.Mrs = handles.Data(:,6);
handles.Bc = handles.Data(:,7);
handles.Brh = handles.Data(:,8);
handles.Shape = handles.Data(:,12);

handles.X1 = handles.Brh ./ handles.Bc;
handles.Y = handles.Mrs ./ handles.Ms;

% Set the logic for the table check boxes
handles.Table_Logic = ones(handles.nData, 1);

table_vals = [cellfun(@(x) sprintf('%3.1f', x), num2cell(handles.Brh), 'UniformOutput', false),...
    cellfun(@(x) sprintf('%3.1f', x), num2cell(handles.Bc), 'UniformOutput', false),...
    cellfun(@(x) sprintf('%3.1f', x), num2cell(handles.X1), 'UniformOutput', false),...
    cellfun(@(x) sprintf('%1.3e', x), num2cell(handles.Mrs), 'UniformOutput', false),...
    cellfun(@(x) sprintf('%1.3e', x), num2cell(handles.Ms), 'UniformOutput', false),...
    cellfun(@(x) sprintf('%1.3f', x), num2cell(handles.Y), 'UniformOutput', false)];


% Set the table
set(handles.Spec_Table, 'Data', [handles.All_Names, table_vals]);


% Set the plot control panel
handles.All_Symbols = {'^', 'o', 'd', '*', 's', '.', '+'};

% Get the input symbol propoerties to set the deafults
% set the radiobutton
Symb_ind=find(strcmpi(handles.Plot_Symbol, handles.All_Symbols));
Symb_Name = strcat('Sym', sprintf('%d', Symb_ind));
set(handles.(Symb_Name), 'Value', 1);

% Set the filled checkbox
if strcmpi(handles.Face_Color, 'none')
    set(handles.CB_Filled, 'Value', 0);
else
    set(handles.CB_Filled, 'Value', 1);
end

% Set the symbol size
set(handles.Size_Select, 'Value', handles.Symbol_Size - 2 );

guidata(hObject, handles);

% Set the current data to all the data and then plot
Update_Plots(handles);



% Update handles structure
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = Hysteresis_BiPlot_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in Plot_Type.
function Plot_Type_Callback(hObject, eventdata, handles)
% hObject    handle to Plot_Type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Plot_Type contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Plot_Type

Update_Plots(handles);

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function Plot_Type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Plot_Type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Save_Plot.
function Save_Plot_Callback(hObject, eventdata, handles)
% hObject    handle to Save_Plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[file,path] = uiputfile('Hysteresis_BiPlot.eps','Save the bi-plot...');

if ~ischar(file) && file==0
    % User has cancelled
    % Do nothing and...
    return;
end

tmpFig=figure('Visible', 'off', 'Units', 'Centimeters');
oldPos=get(tmpFig, 'Position');
set(tmpFig, 'Position', [oldPos(1), oldPos(2), 7.5, 7.5]); % make the figure bigger than needed (300x300)

% Copy and adjust the axes
newAxes=copyobj(handles.Plot_Axes, tmpFig);
set(newAxes, 'Units', 'Centimeters');
axis(newAxes, 'square');
set(newAxes, 'FontUnits', 'Points', 'FontSize', 9)
set(get(newAxes, 'XLabel'), 'FontUnits', 'Points', 'FontSize', 10)
set(get(newAxes, 'YLabel'), 'FontUnits', 'Points', 'FontSize', 10);
set(get(newAxes, 'Title'), 'FontUnits', 'Points', 'FontSize', 11);

% Readjust the x-axis scale and tickmarks
% set(newAxes, 'Xlim', get(handles.Plot_Axes, 'Xlim'))
% set(newAxes, 'XTick', get(handles.Plot_Axes, 'XTick'))
% set(newAxes, 'XTickLabel', get(handles.Plot_Axes, 'XTickLabel'))


NewPos = [1.5, 1.5, 4.5, 4.5];
set(newAxes, 'Position', NewPos, 'XColor', [1,1,1], 'YColor', [1,1,1], 'Box', 'off', 'TickDir', 'Out');

% Place a new set of axes on top to create the box
if handles.Version <= 8.3 % 2014a and before
    h0 = axes('Units', 'Centimeters', 'Position', NewPos);
    set(h0, 'box', 'on', 'XTick', [], 'YTick', [], 'color', 'none');
else% 2014b and later
    h0=copyobj(newAxes, tmpFig);
    cla(h0);
    set(h0, 'box', 'on', 'XTick', [], 'YTick', [], 'color', 'none');
    set(h0, 'Title', [], 'XLabel', [], 'YLabel', []);
end

print(tmpFig, '-depsc', strcat(path, file));
close(tmpFig);


function Update_Plots(handles)


Color = handles.Symbol_Color;
Size = handles.Symbol_Size;
FaceColor = handles.Face_Color;

if strcmpi(FaceColor, 'none')
        FaceColor = 'none';
else
    FaceColor = Color;
end

cla(handles.Symbol_Axes);

hold(handles.Symbol_Axes, 'on')
plot(handles.Symbol_Axes, 0.5, 7, '^', 'Color', Color, 'MarkerSize', Size, 'MarkerFaceColor', FaceColor);
plot(handles.Symbol_Axes, 0.5, 6, 'o', 'Color', Color, 'MarkerSize', Size, 'MarkerFaceColor', FaceColor);
plot(handles.Symbol_Axes, 0.5, 5, 'd', 'Color', Color, 'MarkerSize', Size, 'MarkerFaceColor', FaceColor);
plot(handles.Symbol_Axes, 0.5, 4, '*', 'Color', Color, 'MarkerSize', Size, 'MarkerFaceColor', FaceColor);
plot(handles.Symbol_Axes, 0.5, 3, 's', 'Color', Color, 'MarkerSize', Size, 'MarkerFaceColor', FaceColor);
plot(handles.Symbol_Axes, 0.5, 2, '.', 'Color', Color, 'MarkerSize', Size, 'MarkerFaceColor', FaceColor);
plot(handles.Symbol_Axes, 0.5, 1, '+', 'Color', Color, 'MarkerSize', Size, 'MarkerFaceColor', FaceColor);
hold(handles.Symbol_Axes, 'off')
set(handles.Symbol_Axes, 'Xlim', [0, 1], 'Ylim', [0.5, 7.5]);


% Do the data plot
Plot_Type = get(handles.Plot_Type, 'Value');
Plot_Style = get(handles.Plot_Style, 'Value');


FUnits = 'Pixels';
FontSize1 = 12;
% FontSize2 = 14;

switch Plot_Style
    case 1 % Bc
        Xplot = handles.Bc;
        Xstring = 'B_{c} [mT]';
    case 2 % Brh/Bc
        Xplot = handles.X1;
        Xstring = 'B_{rh} / B_{c}';
    case 3 % Bcr/Bc
        % Not yet implemented
        %                 Xplot = handles.X3;
        %         Xstring = 'B_{cr} / B_{c}';
end


plot(handles.Plot_Axes, Xplot, handles.Y, handles.Plot_Symbol,...
    'color', Color, 'MarkerSize', Size, 'MarkerFaceColor', FaceColor);

set(get(handles.Plot_Axes, 'XLabel'), 'String', Xstring, 'FontUnits', FUnits, 'FontSize', FontSize1);
set(get(handles.Plot_Axes, 'YLabel'), 'String', 'M_{rs} / M_{s}', 'FontUnits', FUnits, 'FontSize', FontSize1);
%         set(get(handles.Plot_Axes, 'Title'), 'String', 'Hysteresis Bi-Plot', 'FontUnits', FUnits, 'FontSize', FontSize2);


switch Plot_Type
    case 1 % Linear Scale
        set(handles.Plot_Axes, 'XScale', 'Linear', 'YScale', 'Linear')
        
        
    case 2 % Log-Linear Scale
        set(handles.Plot_Axes, 'XScale', 'Log', 'YScale', 'Linear')
        
    case 3 %Log-Log Scale
        set(handles.Plot_Axes, 'XScale', 'Log', 'YScale', 'Log')
        
end

% Reset the button down function
set(handles.Plot_Axes, 'ButtonDownFcn', {@Plot_Axes_ButtonDownFcn, handles});

% guidata(hObject, handles);



% --- Executes on button press in DB_ME.
function DB_ME_Callback(hObject, eventdata, handles)
% hObject    handle to DB_ME (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
keyboard


% --- Executes on mouse press over axes background.
function Plot_Axes_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to Plot_Axes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

PopOutFigure(handles.Plot_Axes, 'Bi-Plot')


% --- Executes on button press in Set_Color.
function Set_Color_Callback(hObject, eventdata, handles)
% hObject    handle to Set_Color (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get the new color
handles.Symbol_Color = uisetcolor(handles.Symbol_Color);
guidata(hObject, handles);

% Update the plot
Update_Plots(handles)


% --- Executes on selection change in Size_Select.
function Size_Select_Callback(hObject, eventdata, handles)
% hObject    handle to Size_Select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Size_Select contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Size_Select


contents = cellstr(get(hObject,'String'));
handles.Symbol_Size = str2double(contents{get(hObject,'Value')});

guidata(hObject, handles);

% Update the plot
Update_Plots(handles)


% --- Executes during object creation, after setting all properties.
function Size_Select_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Size_Select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in CB_Filled.
function CB_Filled_Callback(hObject, eventdata, handles)
% hObject    handle to CB_Filled (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CB_Filled

if get(hObject,'Value') == 1
    handles.Face_Color = 'filled';
else
    handles.Face_Color = 'none';
end

guidata(hObject, handles);

% Update the plot
Update_Plots(handles)


% --- Executes on button press in Set_Default.
function Set_Default_Callback(hObject, eventdata, handles)
% hObject    handle to Set_Default (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Call the main defaults and update them
Defaults = handles.Defaults;

Defaults.BiPlotFaceColor = handles.Face_Color;
Defaults.BiPlotSymbol = handles.Plot_Symbol;
Defaults.BiPlotSymbolSize = handles.Symbol_Size;
Defaults.BiPlot_Color = handles.Symbol_Color;


% Save them to a user config file
Save_HystLab_User_Defaults(Defaults);

% Save them to the main window appdata
setappdata(handles.MainWindow, 'Defaults', Defaults);


% --- Executes when selected object is changed in Symbol_Select.
function Symbol_Select_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in Symbol_Select
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

switch get(eventdata.NewValue,'Tag') % Get Tag of selected object.
    case 'Sym1'
        handles.Plot_Symbol = '^';
    case 'Sym2'
        handles.Plot_Symbol = 'o';
    case 'Sym3'
        handles.Plot_Symbol = 'd';
    case 'Sym4'
        handles.Plot_Symbol = '*';
    case 'Sym5'
        handles.Plot_Symbol = 's';
    case 'Sym6'
        handles.Plot_Symbol = '.';
    case 'Sym7'
        handles.Plot_Symbol = '+';
end

guidata(hObject, handles);

Update_Plots(handles);


% --- Executes on button press in Save_Data.
function Save_Data_Callback(hObject, eventdata, handles)
% hObject    handle to Save_Data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


Header = [{'Specimen'}, {'Brh [mT]'}, {'Bc [mT]'},  {'Brh/Bc'},...
    {'Mrs [Am^2]'}, {'Ms [Am^2]'}, {'Mrs/Ms'}];

Data = [handles.All_Names,...
    num2cell([handles.Brh, handles.Bc, handles.X1, handles.Mrs, handles.Ms, handles.Y])]';

% keyboard


[file,path] = uiputfile(strcat('Hysteresis_Bi-Plot_Data.dat'), 'Save bi-plot data...');

if ~ischar(file) && file==0
    % User has cancelled
    % Do nothing and...
    return;
end


FID = fopen(strcat(path, file), 'wt');

fprintf(FID, '%s\t%s\t%s\t%s\t%s\t%s\t%s\n', Header{:});
fprintf(FID, '%s\t%3.1f\t%3.1f\t%3.1f\t%1.3e\t%1.3e\t%1.3f\n', Data{:});

fclose(FID);




% --- Executes on selection change in Plot_Style.
function Plot_Style_Callback(hObject, eventdata, handles)
% hObject    handle to Plot_Style (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Plot_Style contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Plot_Style


Update_Plots(handles);

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function Plot_Style_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Plot_Style (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
