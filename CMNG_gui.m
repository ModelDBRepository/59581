function varargout = CMNG_gui(varargin)
%CMNG_GUI M-file for CMNG_gui.fig
%      CMNG_GUI, by itself, creates a new CMNG_GUI or raises the existing
%      singleton*.
%
%      H = CMNG_GUI returns the handle to a new CMNG_GUI or the handle to
%      the existing singleton*.
%
%      CMNG_GUI('Property','Value',...) creates a new CMNG_GUI using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to CMNG_gui_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      CMNG_GUI('CALLBACK') and CMNG_GUI('CALLBACK',hObject,...) call the
%      local function named CALLBACK in CMNG_GUI.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help CMNG_gui

% Last Modified by GUIDE v2.5 26-Oct-2005 17:33:08

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @CMNG_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @CMNG_gui_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
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


% --- Executes just before CMNG_gui is made visible.
function CMNG_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for CMNG_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes CMNG_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% Set data to predefined values
initialize_gui(handles);


% --- Outputs from this function are returned to the command line.
function varargout = CMNG_gui_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function sim_time_Callback(hObject, eventdata, handles)
% hObject    handle to sim_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sim_time as text
%        str2double(get(hObject,'String')) returns contents of sim_time as a double


% --- Executes during object creation, after setting all properties.
function sim_time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sim_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dt_Callback(hObject, eventdata, handles)
% hObject    handle to dt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dt as text
%        str2double(get(hObject,'String')) returns contents of dt as a double


% --- Executes during object creation, after setting all properties.
function dt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function alpha_CreateFcn(hObject, eventdata, handles)
% hObject    handle to alpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in run_button.
function run_button_Callback(hObject, eventdata, handles)
% hObject    handle to run_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get parameter values
% simulation parameters
simp.tmax = str2double(get(handles.sim_time, 'String'));
simp.dt = str2double(get(handles.dt, 'String'));;
simp.N = str2double(get(handles.N, 'String'));
simp.kmax = str2double(get(handles.kmax, 'String'));
simp.mc = str2double(get(handles.mc, 'String'));
simp.ml = str2double(get(handles.ml, 'String'));
simp.datat = str2double(get(handles.datat, 'String'));
% model parameters
modp.D = str2double(get(handles.D, 'String'));
modp.a = str2double(get(handles.a, 'String'));
modp.g = str2double(get(handles.g, 'String'));
modp.rg = str2double(get(handles.rg, 'String'));
modp.sg = str2double(get(handles.sg, 'String'));
modp.e0 = str2double(get(handles.e0, 'String'));
GS = str2double(get(handles.GS, 'String'));
modp.l0 = str2double(get(handles.l0, 'String'));
modp.c0 = str2double(get(handles.c0, 'String'));
theta = str2double(get(handles.theta, 'String'));
modp.el = str2double(get(handles.el, 'String'));
modp.el = str2double(get(handles.el, 'String'));
alphath = modp.g*modp.sg/(modp.e0*modp.c0*modp.rg*modp.a);  % alpha_twid_h
modp.er = theta*modp.e0;        % soma tubulin autoregulation
modp.rdt = 0;                   % autoregulation time delay
modp.el = GS*modp.rg;           % growth cone flux-sink rate
modp.zl = GS*modp.sg;           % growth cone flux-source rate
set(handles.alpha, 'String', alphath);
set(handles.el, 'String', modp.el);
set(handles.zl, 'String', modp.el);

% Simulation with linear ICs
% calculated parameters
[calcp] = CMNG_calcparams(simp, modp);
% run model for jmax time steps
[C1, C01, CN1, l1] = CMNG_run(simp, modp, calcp, -1, modp);
[t, C1, C01, CN1, l1] = CMNG_dimen(simp, modp, C1, C01, CN1, l1);  % dimensionalise
Ca1 = [C01 C1 CN1];

% Plot results
plot(handles.axes1,t,l1,'k-');
hold on;
title(handles.axes1,'Length');
xlabel(handles.axes1,'Time');
ylabel(handles.axes1,'Length');

space=0:1/simp.N:1;
surf(handles.axes2,space,t,Ca1);
title(handles.axes2,'Concentration over Unit Space and Time');
xlabel(handles.axes2,'Space');
ylabel(handles.axes2,'Time');
zlabel(handles.axes2,'Concentration');
view(handles.axes2,90-37.5,30);


% --- Executes on button press in reset_button.
function reset_button_Callback(hObject, eventdata, handles)
% hObject    handle to reset_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Reset data to predefined values
initialize_gui(handles);


function N_Callback(hObject, eventdata, handles)
% hObject    handle to N (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of N as text
%        str2double(get(hObject,'String')) returns contents of N as a double


% --- Executes during object creation, after setting all properties.
function N_CreateFcn(hObject, eventdata, handles)
% hObject    handle to N (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function kmax_Callback(hObject, eventdata, handles)
% hObject    handle to kmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of kmax as text
%        str2double(get(hObject,'String')) returns contents of kmax as a double


% --- Executes during object creation, after setting all properties.
function kmax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to kmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mc_Callback(hObject, eventdata, handles)
% hObject    handle to mc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mc as text
%        str2double(get(hObject,'String')) returns contents of mc as a double


% --- Executes during object creation, after setting all properties.
function mc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ml_Callback(hObject, eventdata, handles)
% hObject    handle to ml (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ml as text
%        str2double(get(hObject,'String')) returns contents of ml as a double


% --- Executes during object creation, after setting all properties.
function ml_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ml (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function datat_Callback(hObject, eventdata, handles)
% hObject    handle to datat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of datat as text
%        str2double(get(hObject,'String')) returns contents of datat as a double


% --- Executes during object creation, after setting all properties.
function datat_CreateFcn(hObject, eventdata, handles)
% hObject    handle to datat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function D_Callback(hObject, eventdata, handles)
% hObject    handle to D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of D as text
%        str2double(get(hObject,'String')) returns contents of D as a double

% Update calculated parameter values
update_params(handles);


% --- Executes during object creation, after setting all properties.
function D_CreateFcn(hObject, eventdata, handles)
% hObject    handle to D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function a_Callback(hObject, eventdata, handles)
% hObject    handle to a (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of a as text
%        str2double(get(hObject,'String')) returns contents of a as a double

% Update calculated parameter values
update_params(handles);


% --- Executes during object creation, after setting all properties.
function a_CreateFcn(hObject, eventdata, handles)
% hObject    handle to a (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function g_Callback(hObject, eventdata, handles)
% hObject    handle to g (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of g as text
%        str2double(get(hObject,'String')) returns contents of g as a double

% Update calculated parameter values
update_params(handles);


% --- Executes during object creation, after setting all properties.
function g_CreateFcn(hObject, eventdata, handles)
% hObject    handle to g (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function rg_Callback(hObject, eventdata, handles)
% hObject    handle to rg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rg as text
%        str2double(get(hObject,'String')) returns contents of rg as a double

% Update calculated parameter values
update_params(handles);


% --- Executes during object creation, after setting all properties.
function rg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sg_Callback(hObject, eventdata, handles)
% hObject    handle to sg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sg as text
%        str2double(get(hObject,'String')) returns contents of sg as a double

% Update calculated parameter values
update_params(handles);


% --- Executes during object creation, after setting all properties.
function sg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function e0_Callback(hObject, eventdata, handles)
% hObject    handle to e0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of e0 as text
%        str2double(get(hObject,'String')) returns contents of e0 as a double

% Update calculated parameter values
update_params(handles);


% --- Executes during object creation, after setting all properties.
function e0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to e0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function GS_Callback(hObject, eventdata, handles)
% hObject    handle to GS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GS as text
%        str2double(get(hObject,'String')) returns contents of GS as a double

% Update calculated parameter values
update_params(handles);


% --- Executes during object creation, after setting all properties.
function GS_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function l0_Callback(hObject, eventdata, handles)
% hObject    handle to l0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of l0 as text
%        str2double(get(hObject,'String')) returns contents of l0 as a double

% Update calculated parameter values
update_params(handles);


% --- Executes during object creation, after setting all properties.
function l0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to l0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function c0_Callback(hObject, eventdata, handles)
% hObject    handle to c0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of c0 as text
%        str2double(get(hObject,'String')) returns contents of c0 as a double

% Update calculated parameter values
update_params(handles);


% --- Executes during object creation, after setting all properties.
function c0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to c0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function el_CreateFcn(hObject, eventdata, handles)
% hObject    handle to el (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function zl_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



function theta_Callback(hObject, eventdata, handles)
% hObject    handle to theta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of theta as text
%        str2double(get(hObject,'String')) returns contents of theta as a double

% Update calculated parameter values
update_params(handles);


% --- Executes during object creation, after setting all properties.
function theta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to theta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- AUTHOR-DEFINED UTILITY FUNCTIONS

% -- Sets up default parameter values
function initialize_gui(handles)

% simulation parameters
simp.dt = 0.01;                 % time step
simp.tmax = 100;                % simulation time
simp.datat = 100;               % data collection time step
simp.N = 100;                   % number of spatial points
simp.kmax = 10000;              % maximum corrector steps
simp.mc = 0.0001;               % tolerance on C;
simp.ml = 0.0001;               % tolerance on l;
set(handles.sim_time, 'String', simp.tmax);
set(handles.dt, 'String', simp.dt);
set(handles.datat, 'String', simp.datat);
set(handles.N, 'String', simp.N);
set(handles.kmax, 'String', simp.kmax);
set(handles.mc, 'String', simp.mc);
set(handles.ml, 'String', simp.ml);

% model parameters
modp.c0 = 10;                   % concentration scale
modp.l0 = 0.01;                 % initial (min) length;
modp.D = 30000;                 % diffusion constant
modp.a = 100;                   % active transport rate
modp.g = 0.002;                 % decay rate
modp.rg = 10;                   % growth rate constant
modp.sg = 100;                  % growth rate set point (threshold)
GS = 0.00001;
modp.e0 = 0.000002;             % soma flux-source rate (small growth)
alphath = modp.g*modp.sg/(modp.e0*modp.c0*modp.rg*modp.a);  % alpha_twid_h
theta = 0;                      % fractional autoregulation
modp.er = theta*modp.e0;        % soma tubulin autoregulation
modp.rdt = 0;                   % autoregulation time delay
modp.el = GS*modp.rg;           % growth cone flux-sink rate
modp.zl = GS*modp.sg;           % growth cone flux-source rate
set(handles.D, 'String', modp.D);
set(handles.a, 'String', modp.a);
set(handles.g, 'String', modp.g);
set(handles.rg, 'String', modp.rg);
set(handles.sg, 'String', modp.sg);
set(handles.GS, 'String', GS);
set(handles.alpha, 'String', alphath);
set(handles.e0, 'String', modp.e0);
set(handles.el, 'String', modp.el);
set(handles.zl, 'String', modp.el);
set(handles.l0, 'String', modp.l0);
set(handles.c0, 'String', modp.c0);



% -- Updates calculated parameter values
function update_params(handles)

% Get parameter values
% simulation parameters
simp.tmax = str2double(get(handles.sim_time, 'String'));
simp.dt = str2double(get(handles.dt, 'String'));;
simp.N = str2double(get(handles.N, 'String'));
simp.kmax = str2double(get(handles.kmax, 'String'));
simp.mc = str2double(get(handles.mc, 'String'));
simp.ml = str2double(get(handles.ml, 'String'));
simp.datat = str2double(get(handles.datat, 'String'));
% model parameters
modp.D = str2double(get(handles.D, 'String'));
modp.a = str2double(get(handles.a, 'String'));
modp.g = str2double(get(handles.g, 'String'));
modp.rg = str2double(get(handles.rg, 'String'));
modp.sg = str2double(get(handles.sg, 'String'));
modp.e0 = str2double(get(handles.e0, 'String'));
GS = str2double(get(handles.GS, 'String'));
modp.l0 = str2double(get(handles.l0, 'String'));
modp.c0 = str2double(get(handles.c0, 'String'));
theta = str2double(get(handles.theta, 'String'));
modp.el = str2double(get(handles.el, 'String'));
modp.el = str2double(get(handles.el, 'String'));
alphath = modp.g*modp.sg/(modp.e0*modp.c0*modp.rg*modp.a);  % alpha_twid_h
modp.er = theta*modp.e0;        % soma tubulin autoregulation
modp.rdt = 0;                   % autoregulation time delay
modp.el = GS*modp.rg;           % growth cone flux-sink rate
modp.zl = GS*modp.sg;           % growth cone flux-source rate
set(handles.alpha, 'String', alphath);
set(handles.el, 'String', modp.el);
set(handles.zl, 'String', modp.el);
