function varargout = Sun_Synchronous_GUI(varargin)
% SUNSYNCH MATLAB code for sunsynch.fig
%      SUNSYNCH, by itself, creates a new SUNSYNCH or raises the existing
%      singleton*.
%
%      H = SUNSYNCH returns the handle to a new SUNSYNCH or the handle to
%      the existing singleton*.
%
%      SUNSYNCH('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SUNSYNCH.M with the given input arguments.
%
%      SUNSYNCH('Property','Value',...) creates a new SUNSYNCH or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before sunsynch_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to sunsynch_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help sunsynch

% Last Modified by GUIDE v2.5 12-Nov-2016 13:40:35

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @sunsynch_OpeningFcn, ...
                   'gui_OutputFcn',  @sunsynch_OutputFcn, ...
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





% --- Executes just before sunsynch is made visible.
function sunsynch_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to sunsynch (see VARARGIN)

% Choose default command line output for sunsynch
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes sunsynch wait for user response (see UIRESUME)
% uiwait(handles.figure1);

%%%%%%%%%%%%%%%%%%%%%%%%%  Earth (Right) Axis  %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load map data
if exist('earthymap','file')==2
    earthmap = imread('earthymap.png');
else
    load topo
    imwrite(topo,'earthymap.png')
    earthmap = topo;
end

% this if/else statement improves calc time by a facor of 10
% as compared to just "load topo"

% earth = referenceEllipsoid('earth','m');
% % Setting up Earth axis
% axes(handles.earth_axes) 
% handles.earth_map = axesm ('globe','Grid', 'on','Geoid',earth);
% axis(handles.earth_axes,'vis3d');
% view(handles.earth_axes,90,27);
% axis(handles.earth_axes,'off');
% earth_surf = meshm(topo, topolegend, size(topo));
% demcmap(topo);

% setting up Earth axis on right w less maps
axes(handles.earth_axes)
[x_er,y_er,z_er] = sphere(handles.earth_axes,80);          % create a sphere 
Re = 6.371E6;
earth_surf = surf(x_er*Re,y_er*Re,z_er*Re);            % plot spherical surface
zlimits = [min(earthmap(:)) max(earthmap(:))];  %set the min and max values of the color map
demcmap(zlimits);
earth_surf.CData = earthmap;                % set color data to topographic data
earth_surf.FaceColor = 'texturemap';    % use texture mapping
earth_surf.EdgeColor = 'none';          % remove edges
earth_surf.FaceLighting = 'gouraud';    % preferred lighting for curved surfaces

% lighting the earth
hold on
x_light = -1;
y_light = .3;
light(handles.earth_axes,'Position',[x_light y_light .24],...
    'Color','white','Style','infinite');
earth_surf.AmbientStrength = 0.12;
earth_surf.SpecularStrength = .2;
earth_surf.DiffuseStrength = 1;

axis(handles.earth_axes,'vis3d','off'); % makes it rotate without deforming
view(handles.earth_axes,-20,0);

% setting up rotation transform
t1 = hgtransform(handles.earth_axes);
set(earth_surf,'Parent',t1);

%plot3(7e6,7e6,7e6,'MarkerSize',8,'MarkerFaceColor','y')

% rotation of earth in deg/sec
w_earth = 360/(24*3600);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Sun and Earth Axis %%%%%%%%%%%%%%%%%%%%%%%%%%

% Setting up Sun/Earth axis
axes(handles.solar_axes)
axis(handles.solar_axes,'vis3d')
rotate3d on;
axis equal

[x_sphere,y_sphere,z_sphere] = sphere();

% the sun
sun_factor = 8;
x_sun = x_sphere*sun_factor;
y_sun = y_sphere*sun_factor;
z_sun = z_sphere*sun_factor;

handles.sun = surf(handles.solar_axes,x_sun,y_sun,z_sun);
set(handles.sun,'FaceColor','y','EdgeColor','none');

% the earth
r_earth = 20;
x_e = x_sphere+r_earth;
y_e = y_sphere+r_earth;
z_e = z_sphere;
hold on;
handles.earth = surf(handles.solar_axes,x_e,y_e,z_e);
set(handles.earth,'FaceColor','b','EdgeColor','none');

set(handles.solar_axes,'XLim', [-r_earth*1.5 r_earth*1.5],...
    'YLim', [-r_earth*1.5 r_earth*1.5], 'ZLim',[-r_earth*1.5 r_earth*1.5])

% setting up rotation transform for earth around sun
t2 = hgtransform(handles.solar_axes);
set(handles.earth,'Parent',t2);

%noon vectors
earth_quiv = quiver(handles.earth_axes,0,0,x_light*1e7,y_light*1e7); % for right axis
earth_quiv.Color = [0.8500 0.3250 0.0980];


solar_quiv = quiver(handles.solar_axes,r_earth,r_earth,-5,-5);
solar_quiv.Color = [0.8500 0.3250 0.0980];
set(solar_quiv,'Parent',t2);
solar_quiv.MaxHeadSize = 1;


% calculates lat and lon of sat orbit
handles.i = str2double(get(handles.i_edit,'String'));
handles.a = str2double(get(handles.a_edit,'String'))*1000+Re;
handles.e = str2double(get(handles.e_edit,'String'));
handles.RAAN = str2double(get(handles.RAAN_edit,'String'));
handles.w = str2double(get(handles.perigee_edit,'String'));

handles.num_P=20;
handles.dt=20;

[lat, lon, alt, P,vel] = Sat_orbit( handles.a,handles.e,handles.i,...
    handles.RAAN,0,handles.w,10*handles.num_P,0,handles.dt);

alt=alt*Re;

earth = referenceEllipsoid('earth','m');
[handles.ox,handles.oy,handles.oz] = geodetic2ecef(earth,lat,lon,alt+Re);

% set up satellite
handles.hsat = plot3(handles.earth_axes,handles.ox(1), handles.oy(1),...
    handles.oz(1),'Marker','o','MarkerFaceColor','k','MarkerSize',8);

% plot satellite orbit on globe invisibly so axis doesnt shift
    % during animation
endy = round(P/handles.dt);
plot3(handles.earth_axes,handles.ox(1:endy),handles.oy(1:endy),...
    handles.oz(1:endy),'LineStyle','none','Marker','none');
axis(handles.earth_axes,'tight');

  
% time
t = 0:handles.dt:round(10*handles.num_P*P);
del_N = handles.num_P*P/3600*15;
x = del_N/360;

ff = round(handles.a/6e6);
fudgefactor = -4+ff;
tail_length = round((handles.num_P*P+x*handles.num_P*P)/handles.dt)+fudgefactor;
if handles.num_P ==1
    tail_length = tail_length+5;
end
    
% k indexes time
k=0;
% save handles
guidata(hObject,handles)
while k <=length(t)
    k=k+1;
    %%%%%% Earth (Right) Axis %%%%
    % make earth rotate
    n = w_earth/180*pi*t(k);
    Txy = makehgtform('zrotate',n);
    set(t1,'Matrix',Txy);
    
    % if beginning
    if k<tail_length 
        start=1;
    % else if middle
    else
        start=start+1;
    end
    stop=k;
    
    % plot satellite
    set(handles.hsat,'XData',handles.ox(stop),'YData',handles.oy(stop),...
        'ZData',handles.oz(stop));
    
    %%%%%% Sun/ Earth (Left) Axis %%%%%
    n2 = 1.9910e-07*t(k)*1e4;
    Txy2 = makehgtform('zrotate',n2);
    set(t2,'Matrix',Txy2);
    axis(handles.solar_axes,'off')

    
    drawnow
end


% --- Outputs from this function are returned to the command line.
function varargout = sunsynch_OutputFcn(~, ~, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function i_edit_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to scanrate_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of scanrate_edit as text
%        str2double(get(hObject,'String')) returns contents of scanrate_edit as a double
i = str2double(get(hObject,'String'));
handles.i = i;
[lat, lon, alt, P,vel] = Sat_orbit( handles.a,handles.e,handles.i,...
    handles.RAAN,0,handles.w,10*handles.num_P,0,handles.dt);
Re= 6.371E6;
alt=alt*Re;

earth = referenceEllipsoid('earth','m');
[handles.ox,handles.oy,handles.oz] = geodetic2ecef(earth,lat,lon,alt+Re);

guidata(hObject,handles)




