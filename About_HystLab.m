function varargout = About_HystLab(varargin)
% ABOUT_HYSTLAB MATLAB code for About_HystLab.fig
%      ABOUT_HYSTLAB, by itself, creates a new ABOUT_HYSTLAB or raises the existing
%      singleton*.
%
%      H = ABOUT_HYSTLAB returns the handle to a new ABOUT_HYSTLAB or the handle to
%      the existing singleton*.
%
%      ABOUT_HYSTLAB('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ABOUT_HYSTLAB.M with the given input arguments.
%
%      ABOUT_HYSTLAB('Property','Value',...) creates a new ABOUT_HYSTLAB or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before About_HystLab_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to About_HystLab_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help About_HystLab

% Last Modified by GUIDE v2.5 21-Aug-2017 17:50:49

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @About_HystLab_OpeningFcn, ...
                   'gui_OutputFcn',  @About_HystLab_OutputFcn, ...
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


% --- Executes just before About_HystLab is made visible.
function About_HystLab_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to About_HystLab (see VARARGIN)

% Choose default command line output for About_HystLab
handles.output = hObject;

Version = find(strcmp(varargin, 'Version'));
Date = find(strcmp(varargin, 'Date'));

if isempty(Version) || isempty(Date)
    error('About_HystLab:Input', 'Not input version number or date provided');
end

Version = varargin{Version+1};
Date = varargin{Date+1};

Logo = imread('HL_Logo.png');
image(Logo,'Parent',handles.Logo_Axes)
axis(handles.Logo_Axes, 'off');

% Set the text
set(handles.txt_Title, 'string', ['HystLab v' Version]);
set(handles.txt_Date, 'string', Date);
set(handles.txt_MSG, 'string',...
    'Thank you for using HystLab. If you found it useful and you use it in your work, we would be very grateful if you cited the following reference: ');
set(handles.txt_Ref, 'string', 'Paterson, G. A., Zhao, X., Jackson, M., & Heslop, D. (2018). Measuring, processing, and analyzing hysteresis data. Geochemistry, Geophysics, Geosystems, 19. doi: 10.1029/2018GC007620');
set(handles.txt_URL, 'string', [{'The latest version of HystLab is avaiable at:'}, {'https://github.com/greigpaterson/HystLab'}]);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes About_HystLab wait for user response (see UIRESUME)
% uiwait(handles.About_HystLab_Fig);


% --- Outputs from this function are returned to the command line.
function varargout = About_HystLab_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
