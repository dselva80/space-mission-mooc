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

%setting up map info for Earth
% nasa = wmsfind('nasa','SearchField','serverurl');
% layer = nasa.refine('bluemarbleng','SearchField','layername','MatchType','exact');
% [A,R] = wmsread(layer);

load topo;

earth = referenceEllipsoid('earth','m');

[x_sphere,y_sphere,z_sphere] = sphere();


% Setting up Earth axes
axes(handles.earth_axes) 
handles.earth_map = axesm ('globe','Grid', 'on','Geoid',earth);
axis(handles.earth_axes,'vis3d');
view(handles.earth_axes,90,27);
axis(handles.earth_axes,'off');
% showing earths surface w/ map
%handles.hgeo = geoshow(handles.earth_axes,A,R);
earth_surf = meshm(topo, topolegend, size(topo));
demcmap(topo);

% setting up rotation transform
t1 = hgtransform(handles.earth_axes);
set(earth_surf,'Parent',t1);

plot3(7e6,7e6,7e6,'MarkerSize',8,'MarkerFaceColor','y')
%lighting the earth
lightm(0,-90,24933,'color','white','style','infinite');
earth_surf.AmbientStrength = 0.25;
earth_surf.SpecularStrength = .2;
earth_surf.DiffuseStrength = 1;

% rotation of earth in deg/sec
w_earth = 360/(24*3600);
    


% Setting up Sun/Earth axis
axes(handles.solar_axes)
axis(handles.solar_axes,'vis3d')
rotate3d on;
axis equal





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
earth_quiv = quiver(handles.earth_axes,0,0,0,-1e7);

solar_quiv = quiver(handles.solar_axes,r_earth,r_earth,-5,-5);
solar_quiv.Color = [0.8500 0.3250 0.0980];
set(solar_quiv,'Parent',t2);
solar_quiv.MaxHeadSize = 1;

% time
dt=20;
t=0:dt:100000;

% k indexes time
k=0;
while k <=length(t)-1
    k=k+1;
    n = w_earth/180*pi*t(k);
    Txy = makehgtform('zrotate',n);
    set(t1,'Matrix',Txy);
    
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
