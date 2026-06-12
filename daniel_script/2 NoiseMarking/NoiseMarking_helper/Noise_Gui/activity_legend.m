function varargout = activity_legend(varargin)
% ACTIVITY_LEGEND MATLAB code for activity_legend.fig
%      ACTIVITY_LEGEND, by itself, creates a new ACTIVITY_LEGEND or raises the existing
%      singleton*.
%
%      H = ACTIVITY_LEGEND returns the handle to a new ACTIVITY_LEGEND or the handle to
%      the existing singleton*.
%
%      ACTIVITY_LEGEND('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ACTIVITY_LEGEND.M with the given input arguments.
%
%      ACTIVITY_LEGEND('Property','Value',...) creates a new ACTIVITY_LEGEND or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before activity_legend_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to activity_legend_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help activity_legend

% Last Modified by GUIDE v2.5 26-Sep-2020 14:58:28

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @activity_legend_OpeningFcn, ...
                   'gui_OutputFcn',  @activity_legend_OutputFcn, ...
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


% --- Executes just before activity_legend is made visible.
function activity_legend_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to activity_legend (see VARARGIN)

% Choose default command line output for activity_legend
figure(handles.figure1)
movegui(gcf,'east');
hold on;
handles.output = hObject;
colors = [256, 0, 0; % red artifact
          256, 0, 256; % Hypopnea / unsure ( <50% ubstruction), apnea
          128, 128, 128; % movement
          0, 256, 256; % desaturation
          0, 0, 256; % other
          ];
colors = colors / 256;
h = zeros(5, 1);

h(1) = scatter(NaN,NaN,5,colors(1,:), 'square', 'filled');
h(2) = scatter(NaN,NaN,5,colors(2,:), 'square', 'filled');
h(3) = scatter(NaN,NaN,5,colors(3,:), 'square', 'filled');
h(4) = scatter(NaN,NaN,5,colors(4,:), 'square', 'filled');
h(5) = scatter(NaN,NaN,5,colors(5,:), 'square', 'filled');
lgd = legend({'Artifact', 'Hypopnea/apena', 'Movement', 'SpO2 Desaturation', 'Other'});
title(lgd,'Activity Markings')
disableDefaultInteractivity(gca)
set(gca, 'XTick',[])
set(gca, 'YTick',[])

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes activity_legend wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = activity_legend_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
