function varargout = optical_remote_sensing(varargin)
% OPTICAL_REMOTE_SENSING MATLAB code for optical_remote_sensing.fig
%      OPTICAL_REMOTE_SENSING, by itself, creates a new OPTICAL_REMOTE_SENSING or raises the existing
%      singleton*.
%
%      H = OPTICAL_REMOTE_SENSING returns the handle to a new OPTICAL_REMOTE_SENSING or the handle to
%      the existing singleton*.
%
%      OPTICAL_REMOTE_SENSING('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in OPTICAL_REMOTE_SENSING.M with the given input arguments.
%
%      OPTICAL_REMOTE_SENSING('Property','Value',...) creates a new OPTICAL_REMOTE_SENSING or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before optical_remote_sensing_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to optical_remote_sensing_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

%profile on
%profiler

% Edit the above text to modify the response to help optical_remote_sensing

% Last Modified by GUIDE v2.5 05-Nov-2016 18:22:27

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @optical_remote_sensing_OpeningFcn, ...
                   'gui_OutputFcn',  @optical_remote_sensing_OutputFcn, ...
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


% --- Executes just before optical_remote_sensing is made visible.
function optical_remote_sensing_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to optical_remote_sensing (see VARARGIN)

% Choose default command line output for optical_remote_sensing
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes optical_remote_sensing wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% legend axis
leg = handles.legend_axis;
text(leg,.77,1.1,'Legend','FontSize',13)

hold(leg,'on');

% plot world horizon on legend axis
x_world = -1:.1:3;
y_world = -.03*(x_world-1).^2+.3;
plot(leg,x_world,y_world,'LineWidth',2)

% plot FOR lines on legend axis
x_for = 0:.1:2.5;
y_for = -.9*abs(x_for-1)+1;
plot(leg,x_for,y_for,'--r')
x_forlabel = .81:.01:1.19;
y_forlabel = (x_forlabel-1).^2+.8;
plot(leg,x_forlabel,y_forlabel,'r')
text(leg,.5,.83,'FOR','Color','r','FontSize',13)
plot(leg,[0 2.5], [.1 .1],'r')

% plot IFOV and FOV lines on legend axis
x_fov = .943:.001:1.057;
y_fov = -16*abs(x_fov-1)+1;
plot(leg,x_fov,y_fov,'k');
x_fovlabel = .97:.01:1.03;
y_fovlabel = (x_fovlabel-1).^2+.48;
plot(leg,x_fovlabel,y_fovlabel,'k');
text(leg,1.13,.48,'IFOV','Color','k','FontSize',13)
x1 = .6;
twist = .16;
width =.8;
y1 = 0.01;
height = .18;
x2 = x1+width;
y2 = y1+height;
plot(leg,[x1 x1+twist x2 x2-twist x1],[y1 y2 y2 y1 y1],'m');

% plot pixels IFOV
po =  .07;
plot(leg,[x1+po x1+twist+.01-po x2-po x2+po-twist-.01 x1+po],[y1+po y2-po y2-po y1+po y1+po],'k');
pwidth = (width-twist)/5;
px = [x1+po+pwidth x1+pwidth-po+twist+.01];
py = [y1+po y2-po];
plot(leg,px, py,'k')
px2 = [x1+po+pwidth*2 x1+pwidth*2-po+twist+.01];
plot(leg,px2, py,'k');
px3 = [x1+po+pwidth*3 x1+pwidth*3-po+twist+.01];
plot(leg,px3, py,'k');
px4 = [x1+po+pwidth*4 x1+pwidth*4-po+twist+.01];
plot(leg,px4, py,'k');



%FOV label
x_fov = .68:.01:1.32;
y_fov = -2.8*abs(x_fov-1)+1;
plot(leg,x_fov,y_fov,'m')
x_forlabel = .79:.01:1.21;
y_forlabel = .5*(x_forlabel-1).^2+.35;
plot(leg,x_forlabel,y_forlabel,'m')
text(leg,.41,.38,'FOV','Color','m','FontSize',13)


% last of legend plotting
axis(leg,[0,2,0,1],'off')
set(leg,'Color','White')

plot(leg,1,1,'sk','MarkerFaceColor','k','MarkerSize',20)


%setting up map info
nasa = wmsfind('nasa','SearchField','serverurl');
layer = nasa.refine('bluemarbleng','SearchField','layername','MatchType','exact');
[A,R] = wmsread(layer);

A= (A+30)*1.5;
A_orig = A;
handles.A_orig = A_orig;
handles.R = R;
%    C = (A+30)*1.5;
C=A;
C(:,:,2) = C(:,:,1)+100;

handles.C = C;
%    C(:,:,3) = C(:,:,3)+80;

% globe on new axis
earth = referenceEllipsoid('earth','m');
%im not sure why but the current axis is messing up here. I got stuck on
%figureing out why, I'll fix it later.
axes(handles.globe_axis)
%handles.globe_axis = axes('Units','Pixels'); 

handles.globe_map = axesm ('globe','Grid', 'on','Geoid',earth);
% parenty = get(handles.globe_map,'Parent')
axis(handles.globe_axis,'vis3d');
view(handles.globe_axis,60,60);
axis(handles.globe_axis,'off');
handles.hgeo = geoshow(handles.globe_axis,A,R);

% create swath (terrain) axis
axes(handles.swath_axis);
%handles.swath_axis = axes('Units','Pixels');
%handles.swath_map = axesm('ortho','Origin',[0 0],'FLatLimit',[-Inf 0], ...
%         'frame','off','grid','off');

box(handles.swath_axis,'off')
set(handles.swath_axis,'color','none');
axis(handles.swath_axis,'on')
axis(handles.swath_axis,'vis3d')
axis(handles.swath_axis,'on')
handles.swathgeo = geoshow(handles.swath_axis,A,R);


% settting up spacecraft axis
axes(handles.spacecraft_axis);
handles.spacecraft_map = axesm('ortho','Origin',[0 0],'FLatLimit',[-Inf 0], ...
         'frame','off','grid','off');
view(handles.spacecraft_axis,60,27);
handles.spacecraftgeo = geoshow(handles.spacecraft_axis,A,R);
% axis tight
zoom(handles.spacecraft_axis,2);
box(handles.spacecraft_axis,'off');
axis(handles.spacecraft_axis,'off');


% create a variable that contains 1 if animation loop has been called
handles.donethat = 0;

% Update handles structure
guidata(hObject, handles);


% Make the UI visible.
%f.Visible = 'on';


% --- Outputs from this function are returned to the command line.
function varargout = optical_remote_sensing_OutputFcn(~, ~, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



% Diameter
function diameter_edit_Callback(~, ~, ~) %#ok<DEFNU>
% hObject    handle to diameter_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of diameter_edit as text
%        str2double(get(hObject,'String')) returns contents of diameter_edit as a double






% --- Executes during object creation, after setting all properties.
% Diameter
function diameter_edit_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to diameter_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% Altitude
function altitude_edit_Callback(~, ~, ~) %#ok<DEFNU>
% hObject    handle to altitude_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of altitude_edit as text
%        str2double(get(hObject,'String')) returns contents of altitude_edit as a double


% --- Executes during object creation, after setting all properties.
function altitude_edit_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to altitude_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% FOV
function FOR_edit_Callback(~,~, ~) %#ok<DEFNU>
% hObject    handle to FOR_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FOR_edit as text
%        str2double(get(hObject,'String')) returns contents of FOR_edit as a double

% --- Executes during object creation, after setting all properties.





% FOV
function FOR_edit_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to FOR_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% Wavelength
function lambda_edit_Callback(~, ~, ~) %#ok<DEFNU>
% hObject    handle to lambda_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lambda_edit as text
%        str2double(get(hObject,'String')) returns contents of lambda_edit as a double


% --- Executes during object creation, after setting all properties.
%Wavelength
function lambda_edit_CreateFcn(hObject,~, ~) %#ok<DEFNU>
% hObject    handle to lambda_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% Ny
function Ny_edit_Callback(~, ~, ~) %#ok<DEFNU>
% hObject    handle to Ny_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Ny_edit as text
%        str2double(get(hObject,'String')) returns contents of Ny_edit as a double


% --- Executes during object creation, after setting all properties.
% Ny
function Ny_edit_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to Ny_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% Nx
function Nx_edit_Callback(~, ~, ~) %#ok<DEFNU>
% hObject    handle to Nx_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Nx_edit as text
%        str2double(get(hObject,'String')) returns contents of Nx_edit as a double


% --- Executes during object creation, after setting all properties.
% Nx
function Nx_edit_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to Nx_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in whiskbroom_button.
function whiskbroom_button_Callback(~, ~, handles) %#ok<DEFNU>
% hObject    handle to whiskbroom_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of whiskbroom_button
set(handles.broom_panel,'Title','Whiskbroom');
set(handles.Ny_edit,'String','12');
set(handles.FOR_text,'Visible','on');
set(handles.FOR_edit,'Visible','on');
set(handles.scanrate_edit,'Visible','on');
set(handles.scan_text,'Visible','on');
set(handles.scanrate_slider,'Visible','on');

% --- Executes on button press in pushbroom_button.
function pushbroom_button_Callback(~, ~, handles) %#ok<DEFNU>
% hObject    handle to pushbroom_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of pushbroom_button
set(handles.broom_panel,'Title','Pushbroom');
set(handles.Ny_edit,'String','6');
set(handles.FOR_text,'Visible','off');
set(handles.FOR_edit,'Visible','off');






function i_edit_Callback(~, ~, ~) %#ok<DEFNU>
% hObject    handle to i_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of i_edit as text
%        str2double(get(hObject,'String')) returns contents of i_edit as a double


% --- Executes during object creation, after setting all properties.
function i_edit_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to i_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function e_edit_Callback(~, ~, ~) %#ok<DEFNU>
% hObject    handle to e_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of e_edit as text
%        str2double(get(hObject,'String')) returns contents of e_edit as a double


% --- Executes during object creation, after setting all properties.
function e_edit_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to e_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function RAAN_edit_Callback(~, ~, ~) %#ok<DEFNU>
% hObject    handle to RAAN_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of RAAN_edit as text
%        str2double(get(hObject,'String')) returns contents of RAAN_edit as a double


% --- Executes during object creation, after setting all properties.
function RAAN_edit_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to RAAN_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function w_edit_Callback(~, ~, ~) %#ok<DEFNU>
% hObject    handle to w_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of w_edit as text
%        str2double(get(hObject,'String')) returns contents of w_edit as a double


% --- Executes during object creation, after setting all properties.
function w_edit_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to w_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function anom_edit_Callback(~, ~, ~) %#ok<DEFNU>
% hObject    handle to anom_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of anom_edit as text
%        str2double(get(hObject,'String')) returns contents of anom_edit as a double




% --- Executes during object creation, after setting all properties.
function anom_edit_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to anom_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function dt_slider_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to dt_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
dt = 20+get(hObject,'Value')*600;
set(handles.dt_edit,'String',num2str(dt));

% --- Executes during object creation, after setting all properties.
function dt_slider_CreateFcn(hObject, ~,~) %#ok<DEFNU>
% hObject    handle to dt_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function dt_edit_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to dt_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dt_edit as text
%        str2double(get(hObject,'String')) returns contents of dt_edit as a double
dt = get(hObject,'String');
dt = str2double(dt);

% for robustness, can put in whatever values you want
if dt<20
    dt = 20;
elseif dt>620
    dt = 620;
end

dt_scaled = (dt-20)/600;
set(handles.dt_slider,'Value',dt_scaled);


% --- Executes during object creation, after setting all properties.
function dt_edit_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to dt_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in anim_button.
function anim_button_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to anim_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if hObject.Value == 1
    set(handles.anim_button,'String','Stop')
    
    %unpacking + defining variables
    
    %radius of earth
    Re = 6.371E6;
    w = get(handles.w_edit,'String');
    w = str2double(w);
    i = get(handles.i_edit,'String');
    i = str2double(i);
    e = get(handles.e_edit,'String');
    e = str2double(e);
    RAAN = get(handles.RAAN_edit,'String');
    RAAN = str2double(RAAN);
    % alt in m
    h = get(handles.altitude_edit,'String');
    h = str2double(h) * 1e3;
    lambda = get(handles.lambda_edit,'String');
    lambda = str2double(lambda) * 1E-9;
    D = get(handles.diameter_edit,'String');
    D = str2double(D);
    dt = get(handles.dt_edit,'String');
    dt = str2double(dt);
    Nb = get(handles.channel_edit,'String');
    Nb = str2double(Nb);
    % Because scan rate is in scan/min
    scanrate = get(handles.scanrate_edit,'String');
    scan_rate = str2double(scanrate);
    % selva says anamoly doesn't need to be an input for this gui
    anom =0;
    bits = get(handles.bits_edit,'String');
    bits = str2double(bits);
    Nx = get(handles.Nx_edit,'String');
    Nx = str2double(Nx);
    Ny = get(handles.Ny_edit,'String');
    Ny = str2double(Ny);
    FOR = get(handles.FOR_edit,'String');
    FOR = deg2rad(str2double(FOR));
    psize = get(handles.p_edit,'String');
    p = 1e-6*str2double(psize);
    focal = get(handles.focal_edit,'String');
    focal = .01*str2double(focal);

    a = h + Re;
    
    %initialize maps, because I can't figure out how to pass maps
    axes(handles.swath_axis);
    swath_map = axesm('ortho','Origin',[0 0],'FLatLimit',[-Inf 0], ...
         'frame','off','grid','off');
    axes(handles.spacecraft_axis);
    spacecraft_map = axesm('ortho','Origin',[0 0],'FLatLimit',[-Inf 0], ...
        'frame','off','grid','off');
    
    % deleteing handles of objects if animate has already been pushed.
    % Would like a better way of doing this but idk if its possible.
    if handles.donethat==1
           delete(handles.hsat);
           delete(handles.hanim_orb);
           delete(handles.circ_swath);
           delete(handles.cone_anim);
           delete(handles.sat_z);
           delete(handles.cone_z);
           delete(handles.diam);
           delete(handles.diamtext);
           delete(handles.gtrack_z);
           set(handles.hgeo,'CData',handles.A_orig);
           set(handles.swathgeo,'CData',handles.A_orig);
           set(handles.spacecraftgeo,'CData',handles.A_orig);
    end
    
    

    % could be changed later
    num_P = 1;
    % calculates lat and lon of 10 orbital periods
    [lat, lon, alt, P,vel] = Sat_orbit( a,e,i,RAAN,anom,w,10*num_P,0,dt);
    
    % calculate longitude, latitude of groundtrack. alt, P, v not used
    [latg, long] = Sat_orbit( a,e,i,RAAN,anom,w,10*num_P,1,dt);
   
    % spatial resolution
    del_xprime = 1.22*h*lambda/D;
    set(handles.del_xp_output,'String',strcat(num2str(del_xprime),' m'));
    
    % ground velocity
    alt=alt*Re;
    v_g = Re*vel(1) / (Re+alt(1));
    v_g_string = strcat(num2str(v_g,'%10.2g'),' m/s');
    set(handles.vg_output,'String',v_g_string);
    
    
    
    
    % rotation of earth in deg/sec
    w_earth = 360/(24*3600);

    % setting up rotation transform
    t1 = hgtransform(handles.globe_axis);
    set(handles.hgeo,'Parent',t1);
    
    % instantaneous field of view
    IFOV = p/focal;
    
    % IFOV output
    IFOV_string = strcat(num2str(IFOV*1e3,'%10.2g'),' mRad');
    set(handles.IFOV_output,'String',IFOV_string);
    
    % ground pizel size
    g_pixelsize = p*h/focal;
    % integration time and data rate
    push = get(handles.pushbroom_button,'Value');
    if push == 1
        % PUSHBROOM
        % calculate int time, data rate
        int_time = g_pixelsize/v_g;
        data_rate = (Nb*Nx*bits)/(g_pixelsize/v_g);
        % update int time output
        int_time = strcat(num2str(int_time,'%10.2g'), ' s');
        set(handles.int_time_output,'String',int_time);
        % update data rate output
        data_string = strcat(num2str(data_rate/1e6,'%10.2f'),' Mbps');
        set(handles.data_rate_output,'String',data_string);
        % calculate swath area
        swathy = 2*h*tan(Nx*IFOV/2);
        set(handles.swath_output,'String',strcat(num2str(swathy,'%10.2g'),' m'));
    else
        % WHISKBROOM
        %calculate int time and data rate
        int_time = Ny*g_pixelsize/(2*Nx*v_g);
        data_rate = (Nb*Nx*Ny*bits)/(g_pixelsize/v_g);
        % update int time output
        int_time = strcat(num2str(int_time,'%10.2g'), ' s');
        set(handles.int_time_output,'String',int_time);
        % update data rate output
        data_string = strcat(num2str(data_rate/1e6,'%10.2f'),' Mbps');
        set(handles.data_rate_output,'String',data_string);
        % calculate swath area
        swathy = 2*h*tan(FOR/2);
        set(handles.swath_output,'String',strcat(num2str(swathy,'%10.2g'),' m'));
    end
    
    r = swathy;
    
     % pixel ground size
    pgs_string = strcat(num2str(g_pixelsize,'%10.2g'),' m');
    set(handles.pixel_g_size_output,'String',pgs_string);
    
    
    % orbit x y z pre calculated
    earth = referenceEllipsoid('earth','m');
    
    axes(handles.globe_axis)
    %handles.globe_axis = axes('Units','Pixels'); 
    handles.globe_map = axesm ('globe','Grid', 'on','Geoid',earth);
    mstruct = gcm(handles.globe_map);
    [ox, oy, oz] = mfwdtran(mstruct,lat,lon,alt);
    % set up satellite
    handles.hsat = plot3(handles.globe_axis,ox(1), oy(1), oz(1),'Marker','o',...
      'MarkerFaceColor','k','MarkerSize',8);
   
    % plots animated orbit
    handles.hanim_orb = plot3(handles.globe_axis,ox(1), oy(1),oz(1),'Color','m','LineWidth',3);
%     hgeo.CData = handles.A_orig;


    % calculates and runs conical swath animation on globe
    FOV = Nx*IFOV;
    
    [xECEF, yECEF, zECEF] = conical_animation(FOV*180/pi/2,alt(1),latg(1),long(1),lat(1),lon(1));

    handles.cone_anim = surf(handles.globe_axis,xECEF,yECEF,zECEF,...
        'Facecolor','m','Facealpha',.4,'LineStyle','none','Visible','off');

    t=0:dt:round(10*num_P*P);
    del_N =  num_P*P/3600*15;
    x = del_N/360;
    % nodal precession doesnt perfectly close up groundtack on itself at
    % larger altitudes (higher P)
    % Introducing fudgefactor which should be improved. this ff is much smaller than
    % the calculate one. not sure why that is.
    ff = round(a/6e6);
    fudgefactor = -4+ff;
    tail_length = round((num_P*P+x*num_P*P)/dt)+fudgefactor;
    if num_P ==1
        tail_length = tail_length+5;
    end
    


    % plot satellite orbit on globe invisibly so axis doesnt shift
    % during animation
    endy = round(P/dt);
    plot3(handles.globe_axis,ox(1:endy),oy(1:endy),oz(1:endy),...
       'LineStyle','none','Marker','none');
    axis(handles.globe_axis,'tight');

   % determine heading
   heading = zeros(length(latg)-1);
   for x=1:length(latg)-1
        slat = [latg(x) latg(x+1)];
        slon = [long(x) long(x+1)];
        [heading(x), ~] = legs(slat,slon,'gc');
   end



   % setting up swath
    ang=0:0.1:2*pi+.1;
    handles.circ_swath = plot(handles.swath_axis,0,0,'m','Visible','off','LineWidth',4);
    handles.diam = plot(handles.swath_axis,0,0,'m','Visible','off','LineWidth',3);
    handles.diamtext = text(handles.swath_axis,0,0,'','Color','m','Visible','off','FontSize',15,...
    'FontWeight','bold');

    % setting up zoom axis
    axis(handles.spacecraft_axis,'tight');
    handles.sat_z = plot3(handles.spacecraft_axis,0,0,0,'ks','Markerfacecolor','k','MarkerSize',15,'Visible','off');
    rix = [0 0;0 0];
    hold on
    handles.gtrack_z = plot3(handles.spacecraft_axis,0,0,0,'c--','LineWidth',3);
    handles.cone_z = surf(handles.spacecraft_axis,rix,rix,rix,'Visible','off','Facecolor','m',...
        'Facealpha',.4,'LineStyle','none');
    zlim(handles.spacecraft_axis, [0 .1]);

    
    % converting cross and along track angles 
    cross = FOR*180/pi;
    r_cross = tand(cross)*alt(1);
    arc_cross = atand(r_cross/Re);
    
    
    
    k=0;
    increasing=1;
    swing=0;
    
    %scan rate stuff before loop
    %convert scan_rate to m/s
    if scan_rate == 0
            % not scanning, stationary
            arc_cross=0;
            anim_step=1;
    else
            anim_step=1/(scan_rate/60*dt);
    end

   
        
    % loop of animation
    while k <=length(t) && hObject.Value==1
        k=k+1;
        %animating=1;

        % n is total radians rotated by earth at time
        % earth rotation
        n = w_earth/180*pi*t(k);
        Txy = makehgtform('zrotate',n);
        set(t1,'Matrix',Txy)

        % if beginning
        if k<tail_length 
            start=1;
        % else if middle
        else
            start=start+1;
        end
        stop=k;

        %%%%%%%% globe


        % plot satellite on globe
        set(handles.hsat,'XData',ox(stop),'YData',oy(stop),...
            'ZData',oz(stop));

        %plot animated orbit on globe
        set(handles.hanim_orb,'XData',ox(start:stop),...
            'YData',oy(start:stop),'ZData',oz(start:stop));
        

        % changing the colormap to represent whats been seen
%             if val == 1
%                 %conical
%                 temp_r = r_nadir + abs(r_actview);
%                 [circle_lat2, circle_lon2] = scircle1(latg(stop),long(stop),temp_r);
%             elseif val == 2
%                 %cross
%                 temp_r = arc_cross;
%                 [circle_lat2, circle_lon2] = scircle1(latg(stop),long(stop),temp_r);
%             else
%                 temp_r = r;
%                 [circle_lat2, circle_lon2] = scircle1(latg(stop),long(stop),temp_r);
%             end
%             [rmin,rmax,cmin,cmax] = coverage_lineECEF( circle_lat2,circle_lon2,temp_r,A,R);
%             
%             if val==2 || val==3
%                 [~,~,cmin,cmax] = coverage_lineECEF( circle_lat2,circle_lon2,arc_along,A,R);
%             end

%             
%             %This is an approx way. i'm thinking to make this better, use
%             % scircle, find rows and cols at 4 quartes of the circle of
%             % view?
%             % or convert to xyz, find min and max  of x,y then convert to
%             % rows,cols?
%              if abs(cmin-cmax) > 512/2;
%                 % if at  edge of image
%                 if latg(stop)+r >= 90
%                     %if at north pole
% %                     disp('north')
%                     hgeo.CData(1:rmax,cmin:cmax,:) = C(1:rmax,cmin:cmax,:);
%                     hgeo.CData(1:rmin,cmin:cmax,:) = C(1:rmin,cmin:cmax,:);
% %                     A(1:rmin,:,:) = C(1:rmin,:,:);
%                 elseif latg(stop)-r <= -90
% %                         disp('south')
%                 else
%                     % else crossing prime meridian
%                     hgeo.CData(rmin:rmax,cmax:512,:) = C(rmin:rmax,cmax:512,:);
%                     hgeo.CData(rmin:rmax,1:cmin,:) = C(rmin:rmax,1:cmin,:);
%                     
%                 end
%             
%              else
%                 hgeo.CData(rmin:rmax,cmin:cmax,:) = C(rmin:rmax,cmin:cmax,:);
%              end



       

        % approximate altitude
        zeta = alt(stop)/1e7;

        zlim(handles.spacecraft_axis, [0 zeta]);

        


        % animation on globe
        if swing >= arc_cross/2
            increasing=0;
        elseif swing <= -arc_cross/2
            increasing=1;
        end
        
        
         

        if increasing == 1
            swing = swing+arc_cross/anim_step;
        else
            swing = swing-arc_cross/anim_step;
        end

        if push == 1
            % if pushbroom
            % change direction of animation
            directi = 180;
            swing = 0;
        else
            % if whiskbroom
            directi = 90;
        end

        [latout,lonout] = reckon(lat(stop),lon(stop),swing,directi);

        %without real stuff it spazzes at the poles
        [xsurf, ysurf, zsurf] = conical_animation(FOV*180/pi,...
            alt(stop),real(latout),real(lonout),...
            real(lat(stop)),real(lon(stop)));

        xsurf = real(xsurf);
        ysurf = real(ysurf);
        zsurf = real(zsurf);

        % flat map of swath
        origin = [round(latout) round(lonout) -heading(stop)];
        %swath_map = axesm('ortho','Origin',origin,'FLatLimit',[-Inf FOV*180/pi/2], ...
        % 'frame','off','grid','off');


        % back to globe
        set(handles.cone_anim,'Visible','on','XData',xsurf,...
            'YData',ysurf,'ZData',zsurf);

        % plot cone on zoom
        originz = [round(latg(stop)) round(long(stop)) -heading(stop)];
        setm(swath_map,'Origin',origin,'FLatLimit',[-Inf FOV*180/pi/5*3]);
        setm(spacecraft_map,'Origin',originz,'FLatLimit',[-Inf FOV*180/pi*2]);

        axis(handles.spacecraft_axis,'off')
        
         %%%%%%%% hzoom axis %%%%%%%%%%%%%

        %%%%% zoom axis
        %handles.spacecraft_axis.Visible = 'on' ;

        xlimits2 = xlim(handles.spacecraft_axis);
        x4 = xlimits2(1);
        x6 = xlimits2(2);
        sc_axislength = x6-x4;
        x5 = x4+(sc_axislength)/2;

        ylimits2 = ylim(handles.spacecraft_axis);
        y4 = ylimits2(1);
        y6 = ylimits2(2);
        y5 = y4+(y6-y4)/2;
        
        % updating sat on spacecraft axis
        set(handles.sat_z,'XData',x5,'YData',y5,'ZData',zeta,'Marker','s',...
            'Visible','on');

        % updating groundtrack on spacecraft axis
        gy = y5-sc_axislength/3+.01:.1:y5+sc_axislength/3+.005;
        gz = zeros(size(gy));
        set(handles.gtrack_z,'XData',gz,'YData',gy,'ZData',gz,'Visible','on');

        % xyz on orthogonal projection are not the same as on the
        % globe. Scaled theglobe ECEF xyz down by a factor of 40
        % and it looks roughly correct
        
        if push == 1
            xcone = FOV*180/pi/80*cos(ang);
            ycone = swing/40 + FOV*180/pi/80*sin(ang);
        else
            xcone= swing/40 + FOV/80*180/pi * cos(ang);
            ycone=  FOV/80*180/pi * sin(ang);
        end
        xcirc2 = xcone+x5;
        ycirc2 = ycone+y5;
        xsurfz = zeros(2,length(xcirc2));
        xsurfz(1,:) = xcirc2; 
        xsurfz(2,:) = x5;

        ysurfz = zeros(2,length(ycirc2));
        ysurfz(1,:) = ycirc2;
        ysurfz(2,:) = y5;

        zsurfz = zeros(2,length(xcirc2));
        zsurfz(1,:) = 0;
        zsurfz(2,:) = zeta;

        set(handles.cone_z,'Visible','on',...
            'XData',xsurfz,'YData',ysurfz,'ZData',zsurfz);




        % getting the limits of the swath
        xlimits = xlim(handles.swath_axis);
        x1 = xlimits(1);
        x2 = x1+(xlimits(2)-xlimits(1))/2;
        x3 = xlimits(2);
        rad = x2-x1;
        ylimits = ylim(handles.swath_axis);
        y1 = ylimits(1);
        y2 = y1+(ylimits(2)-ylimits(1))/2;
        y3 = ylimits(2);

        % setting the lat ticks and tick labels for swath
        set(handles.swath_axis,'XTick',[x1 x2 x3],'YTick',[y1 y2 y3]);
        lattext = strcat(num2str(round(latg(stop))),'{\circ}');
        lat1 = round(latg(stop)-FOR*180/2/pi);
        if lat1<-90
            lat1 = lat1+180;
        end
        lat3 = round(latg(stop)+FOR*180/2/pi);
        if lat3>90
            lat3 =180-lat3;
        end
        lat1text = strcat(num2str(lat1),'{\circ}');
        lat3text = strcat(num2str(lat3),'{\circ}');

        % setting the lon ticks and tick labels for swath
        lontext = strcat(num2str(round(long(stop))),'{\circ}');
        lon1 = round(long(stop)-FOR*180/2/pi);
        if lon1<-180
            lon1 = lon1+360;
        end
        lon3 = round(long(stop)+FOR*180/2/pi);
        if lon3>180
            lon3 =360-lon3;
        end
        lon1text = strcat(num2str(lon1),'{\circ}');
        lon3text = strcat(num2str(lon3),'{\circ}');
        set(swath_map,'XTickLabel',{lon1text,lontext,lon3text},...
            'YTickLabel',{lat1text,lattext,lat3text}); 


        % this fixes the fact that the frame of the ortho projection 
        % is not filled with pixels.
        % diy frame, slightly smaller than the matlab frame would be
        xp=(rad-rad/r)*cos(ang);
        yp=(rad-rad/r)*sin(ang);
        circ_swath.XData = x2+xp;
        circ_swath.YData = y2+yp;
        circ_swath.Visible = 'on';


        % plot diameter and diameter text
        xdiam = x1+1/r*rad:.001:x3-1/r*rad;
        ydiam = ones(size(xdiam))*y2;
        diam.XData = xdiam;
        diam.YData = ydiam;
        diam.Visible='on';

        % updating diameter km text
        str = strcat(num2str(round(abs(r*1E-3))),' km');
        diamtext.Visible = 'on';
        diamtext.String = str;
        diamtext.Position = [x2-x1/8 y2-y1/4];

%             swathgeo.CData(rmin:rmax,cmin:cmax,:) = C(rmin:rmax,cmin:cmax,:);
        %set(handles.hgeo,'CData',handles.A_orig);
        %set(handles.swathgeo,'CData',handles.A_orig);



        drawnow
        
        
    end




else
    set(handles.anim_button,'String','Animate')
end


handles.donethat = 1;

% Update handles structure
guidata(hObject, handles);




function channel_edit_Callback(~, ~, ~) %#ok<DEFNU>
% hObject    handle to channel_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of channel_edit as text
%        str2double(get(hObject,'String')) returns contents of channel_edit as a double


% --- Executes during object creation, after setting all properties.
function channel_edit_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to channel_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bits_edit_Callback(~, ~, ~) %#ok<DEFNU>
% hObject    handle to bits_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bits_edit as text
%        str2double(get(hObject,'String')) returns contents of bits_edit as a double


% --- Executes during object creation, after setting all properties.
function bits_edit_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to bits_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function spacecraft_axis_CreateFcn(~, ~, ~) %#ok<DEFNU>
% hObject    handle to spacecraft_axis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate spacecraft_axis



% --- Executes during object creation, after setting all properties.
function legend_axis_CreateFcn(~, ~, ~) %#ok<DEFNU>
% hObject    handle to legend_axis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% Hint: place code in OpeningFcn to populate legend_axis

function figure1_CreateFcn(~, ~, ~) %#ok<DEFNU>


% --- Executes during object creation, after setting all properties.
function inc_text_CreateFcn(~, ~, ~) %#ok<DEFNU>
% hObject    handle to inc_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



function p_edit_Callback(~, ~, ~) %#ok<DEFNU>
% hObject    handle to p_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of p_edit as text
%        str2double(get(hObject,'String')) returns contents of p_edit as a double


% --- Executes during object creation, after setting all properties.
function p_edit_CreateFcn(hObject,~, ~) %#ok<DEFNU>
% hObject    handle to p_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function focal_edit_Callback(~,~, ~) %#ok<DEFNU>
% hObject    handle to focal_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of focal_edit as text
%        str2double(get(hObject,'String')) returns contents of focal_edit as a double


% --- Executes during object creation, after setting all properties.
function focal_edit_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to focal_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function scanrate_edit_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to scanrate_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of scanrate_edit as text
%        str2double(get(hObject,'String')) returns contents of scanrate_edit as a double
scanrate = get(hObject,'String');
scanrate = str2double(scanrate);

sr_scaled = scanrate*10;
set(handles.scanrate_edit,'Value',sr_scaled);


% --- Executes during object creation, after setting all properties.
function scanrate_edit_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to scanrate_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function scanrate_slider_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to scanrate_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


scanrate = get(hObject,'Value')*6;
set(handles.scanrate_edit,'String',num2str(scanrate));


% --- Executes during object creation, after setting all properties.
function scanrate_slider_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
% hObject    handle to scanrate_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
