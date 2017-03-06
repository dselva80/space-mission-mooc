function varargout = Sun_Synchronous_GUI(varargin)
% SUNSYNCH MATLAB code for sunsynch.fig
%      SUNSYNCH, by itself, creates a new SUNSYNCH or raises the existing
%      singleton*.
%
%      H = SUNSYNCH returns the handle to a new SUNSYNCH or the handle to
%      the existing singleton*.
%
%      SUNSYNCH('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in .M with the given input arguments.
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

set(gcf,'Name','Sun Synchronous Orbit Design');

%%%%%%%%%%%%%%%%%%%%%%%%%  Earth (Right) Axis  %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input from Button
handles.planet = get(handles.planet_selection,'Value');
handles = Update_Planet_Constants(hObject,handles);
% Draw Planet on Right Axis

% setting up Earth axis on right w less maps
axes(handles.earth_axes)

[x,y,z] = sphere(handles.earth_axes,80); % create a sphere, 80 for smoothness

handles.planet_surf = surf(x*handles.R,y*handles.R,z*handles.R,'FaceColor','texturemap',...
    'EdgeColor','none','FaceLighting','gouraud'); % plot spherical surface




% update map and lighting for planet
handles = Draw_Planet(hObject,handles);

% lighting the planet
hold on
light(handles.earth_axes,'Position',handles.light_pos,...
    'Color','white','Style','infinite');
set(handles.planet_surf,'AmbientStrength',.12,'SpecularStrength',.2,...
    'DiffuseStrength',1);

axis(handles.earth_axes,'vis3d','off'); % makes it rotate without deforming
view(handles.earth_axes,-30,24);



% Inputs from Buttons
handles.i = str2double(get(handles.i_edit,'String')); % kept as degrees
handles.a = str2double(get(handles.a_edit,'String'))*1000+handles.R; % converted from alt to semi major axis
handles.e = str2double(get(handles.e_edit,'String'));
handles.w = str2double(get(handles.perigee_edit,'String')); % degrees
handles.RAAN = str2double(get(handles.RAAN_edit,'String')); % kept as degrees

% % setting up rotation transform
t1 = hgtransform(handles.earth_axes);
set(handles.planet_surf,'Parent',t1);

%plot3(7e6,7e6,7e6,'MarkerSize',8,'MarkerFaceColor','y')

% rotation of earth around sun in deg/sec
w_earth = 360/(24*3600);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Sun and Earth Axis %%%%%%%%%%%%%%%%%%%%%%%%%%

% Setting up Sun/Earth axis
axes(handles.solar_axes)
% whitebg([0 0 0])
axis(handles.solar_axes,'vis3d')
rotate3d on;
axis equal

[x_sphere,y_sphere,z_sphere] = sphere();

% the sun
sun_factor = handles.R/1E6;
x_sun = x_sphere*sun_factor;
y_sun = y_sphere*sun_factor;
z_sun = z_sphere*sun_factor;

handles.sun = surf(handles.solar_axes,x_sun,y_sun,z_sun);
set(handles.sun,'FaceColor','y','EdgeColor','none');

% the earth
r_earth = 30;
x_e = x_sphere+r_earth;
y_e = y_sphere+r_earth;
z_e = z_sphere;
hold on;
handles.earth = surf(handles.solar_axes,x_e,y_e,z_e);
set(handles.earth,'FaceColor',handles.colorz,'EdgeColor','none');

set(handles.solar_axes,'XLim', [-r_earth*1.5 r_earth*1.5],...
    'YLim', [-r_earth*1.5 r_earth*1.5], 'ZLim',[-r_earth*1.5 r_earth*1.5])

% setting up rotation transform for earth around sun
t2 = hgtransform(handles.solar_axes);
set(handles.earth,'Parent',t2);

% noon vector on right axis
handles = draw_noon_vectors(hObject,handles,1);

% noon vector left axis
handles.solar_quiv = quiver(handles.solar_axes,r_earth,r_earth,-5,-5,...
    'Color',[0.8500 0.3250 0.0980],'MaxHeadSize',1); 

handles.solar_orb = quiver(handles.solar_axes,r_earth,r_earth,-5,-5,...
    'Color',[0 0 0],'MaxHeadSize',1); 

handles = yearly_sat_orbit(hObject,handles);


set(handles.solar_quiv,'Parent',t2);
set(handles.solar_orb,'Parent',t2);


% calculates lat and lon of sat orbit
handles.num_P=1;
handles.dt=20;

updated_handles2 = Sat_orbit(hObject, handles);
handles=updated_handles2;

  
% time
handles = calculate_t(hObject,handles);


% set up satellite
handles.hsat = plot3(handles.earth_axes,handles.ox(1), handles.oy(1),...
    handles.oz(1),'Marker','o','MarkerFaceColor','k','MarkerSize',8);

% plot satellite orbit on globe invisibly so axis doesnt shift
    % during animation
endy = handles.tail_length;
plot3(handles.earth_axes,handles.ox(1:endy),handles.oy(1:endy),...
    handles.oz(1:endy),'LineStyle','none','Marker','none');
       
axis(handles.earth_axes,'tight');



%set up orbital plane
handles = draw_orbital_plane(hObject,handles,1);

% year time in rad, index of yr = day
handles.yr = 0:2*pi/365:2*pi*20;
    
% k indexes time
handles.k=0;
% save handles
guidata(hObject,handles)
while handles.k <=length(handles.t)
    
    handles = guidata(hObject);
    handles.k = handles.k+1;
    
    %%%%%% Earth (Right) Axis %%%%
    % make earth rotate
    n = w_earth/180*pi*handles.t(handles.k);
    Txy = makehgtform('zrotate',n);
    set(t1,'Matrix',Txy);
    
    %leaving this bit in here in case i decide to plot a tail to the sat
    %(just add start:stop)
    % if beginning
    if handles.k<handles.tail_length 
        start=1;
    % else if middle
    else
        start=start+1;
    end
    stop=handles.k;
    
    % plot satellite
    set(handles.hsat,'XData',handles.ox(stop),'YData',handles.oy(stop),...
        'ZData',handles.oz(stop));
    handles = draw_orbital_plane(hObject,handles,0);
    
    %%%%%% Sun/ Earth (Left) Axis %%%%%
    n2 = handles.yr(handles.k);
    Txy2 = makehgtform('zrotate',n2);
    set(t2,'Matrix',Txy2);
    axis(handles.solar_axes,'off')
    
    h=7;
    newu = -h*sin(handles.Omega_t(handles.k));
    newv = -h*cos(handles.Omega_t(handles.k));
    set(handles.solar_orb,'UData',newu,'VData',newv);
    
    
    
    
    
    % save handles
    guidata(hObject,handles)

    pause(.01)
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
handles = Sat_orbit(hObject,handles);
handles.k=0; %restart loop
guidata(hObject,handles)

function a_edit_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to scanrate_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of scanrate_edit as text
%        str2double(get(hObject,'String')) returns contents of scanrate_edit as a double
a = str2double(get(hObject,'String'));
handles.a = a*1E3+handles.R;
handles = Sat_orbit(hObject,handles);
handles = calculate_t(hObject,handles);
handles.k=1; %restart loop
guidata(hObject,handles)

function e_edit_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to scanrate_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of scanrate_edit as text
%        str2double(get(hObject,'String')) returns contents of scanrate_edit as a double
e = str2double(get(hObject,'String'));
handles.e = e;
handles = Sat_orbit(hObject,handles);
handles.k=0; %restart loop
guidata(hObject,handles)

function w_edit_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to scanrate_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of scanrate_edit as text
%        str2double(get(hObject,'String')) returns contents of scanrate_edit as a double
w = str2double(get(hObject,'String'));
handles.w = w;
handles = Sat_orbit(hObject,handles);
handles.k=0; %restart loop
guidata(hObject,handles)

function RAAN_edit_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to scanrate_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of scanrate_edit as text
%        str2double(get(hObject,'String')) returns contents of scanrate_edit as a double
RAAN = str2double(get(hObject,'String'));
handles.RAAN = RAAN;
handles = Sat_orbit(hObject,handles);
handles.k=0; %restart loop
guidata(hObject,handles)

function Planet_Selection_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns contents...
%        contents{get(hObject,'Value')} returns selected item...
%items = get(hObject,'String');
index_selected = get(hObject,'Value');
handles.planet = index_selected;
handles = Update_Planet_Constants(hObject,handles);
handles = Draw_Planet(hObject,handles);

% update orbit of satellite to account for change in radius of planet
handles.a = str2double(get(handles.a_edit,'String'))*1000+handles.R;
handles = Sat_orbit(hObject,handles);
set(handles.earth,'FaceColor',handles.colorz);

%update time
handles = calculate_t(hObject,handles);

% update solar (left) axis
[xs,ys] = sphere();
r_p = 10*handles.planet;
if r_p>40
    xs=xs*3;
    ys=ys*3;
    set(handles.earth,'FaceColor',[.2 .2 .2]);
end
set(handles.earth,'XData',xs+r_p,'YData',ys+r_p)
% set(handles.solar_quiv,'XData',r_p,'YData',r_p,'UData',-r_p*.2,'VData',-r_p*.2)
set(handles.solar_quiv,'XData',r_p,'YData',r_p)
handles = draw_noon_vectors(hObject,handles,0);

get_MLTAN(hObject,handles);

% limits get wonky
szfctr = 1.5;
low = -r_p*szfctr;
hi = r_p*szfctr;
xlim(handles.solar_axes,[low hi]);
ylim(handles.solar_axes,[low hi]);
zlim(handles.solar_axes,[low hi]);

%restart loop
%handles.k=0;

guidata(hObject,handles)

%this function changes the radius and rotation values held in handles to
%match the planet selection
function updated_handles = Update_Planet_Constants(hObject,handles)
planet = handles.planet;
if planet == 3
    %earth
    handles.mu = 3.986E14;
    handles.J2 = 0.00108263;
    handles.R = 6.371E6;
    handles.name = 'earth';
    % rotation of the Earth [rad/s]
    handles.rot = 7.2921e-5;
    handles.colorz = [0 0 160/255];
elseif planet == 1
    % Mercury
    handles.mu = 2.2032e13;
    handles.J2 = 0.00006;
    handles.R = 2.439E6;
    handles.rot = 1.24e-06;
    handles.colorz = [.6 .6 .6];
    handles.name = 'mercury';
elseif planet == 2
    % Venus
    handles.mu = 3.24859E14;
    handles.J2 = 0.000027;
    handles.R = 6.051E6;
    handles.rot = 2.99E-7;
    handles.colorz = [250/255, 220/255, 121/255];
    handles.name = 'venus';
elseif planet == 4
    % Mars
    handles.mu = 4.282837e13;
    handles.J2 = 0.001964;
    handles.R = 3.396E6;
    handles.rot = 7.088e-05;
    handles.colorz = [.9 .3 0];
    handles.name = 'mars';
elseif planet == 5
    % Jupiter
    handles.mu = 1.26686534E17;
    handles.J2 = 0.01475;
    handles.R = 71.492E6;
    handles.rot = 1.77E-4;
    handles.colorz = [210/255, 200/255, 150/255];
    handles.name = 'jupiter';
elseif planet == 6
    % Saturn
    handles.mu = 3.7931187E16;
    handles.J2 = 0.01645;
    handles.R = 60.268E6;
    handles.rot =1.63E-4;
    handles.map = [210/255, 200/255, 150/255];
    handles.name = 'saturn';
elseif planet == 7
    % Uranus
    handles.mu = 5.794e15;
    handles.J2 = 0.012;
    handles.R = 25.559E6;
    handles.rot = -1.04E-4;
    handles.colorz = [209/255 231/255 231/255];
    handles.name = 'uranus';
elseif planet ==8
    % Neptune
    handles.mu = 6.809e15;
    handles.J2 = 0.0004;
    handles.R = 24.764E6;
    handles.rot = 1.08E-4;
    handles.colorz = [63/255 84/255 186/255];
    handles.name = 'neptune';
else
    %Pluto
    handles.mu = 8.71E11;
    handles.J2 = 0;
    handles.R = 1.195E6;
    handles.rot = -1.29E-5;
    handles.colorz = [.8 .7 .5];
    handles.name = 'pluto';
end
updated_handles=handles;
guidata(hObject,handles)

function updated_handles = Draw_Planet(hObject,handles)


if handles.planet == 1
    % mercury
    if exist('mercury_map','file')==2
        mercury_map = imread('mercury_map.png');
    else
        try
            mercury_map =imread('http://images.spaceref.com/news/2011/ooglobal_view_MESSENGER.jpg','jpg');
        catch
            disp('Issue loading planet surface image from web.');
            mercury_map = rand(1,1);
            colormap(handles.colorz);
        end
        imwrite(mercury_map,'mercury_map.png') 
    end
   map = mercury_map;
    
elseif handles.planet == 2
    % venus
    if exist('venus_map','file')==2
        venus_map = imread('venus_map.png');
    else
        try
            venus_map =imread('http://maps.jpl.nasa.gov/pix/ven0mss2.jpg','jpg');
        catch
            disp('Issue loading planet surface image from web.');
            venus_map = rand(1,1);
            colormap(handles.colorz);
        end
        imwrite(venus_map,'venus_map.png') 
    end
   map = venus_map;
    
elseif handles.planet == 3
    % earth
    if exist('earthymap','file')==2
        earthmap = imread('earthymap.png');
    else
        load topo
        imwrite(topo,'earthymap.png')
        earthmap = topo;
    end
    zlimits = [min(earthmap(:)) max(earthmap(:))];  %set the min and max values of the color map
    demcmap(zlimits);
    map=earthmap;
    
elseif handles.planet == 4
    % mars
    if exist('mars_map','file')==2
        mars_map = imread('mars_map.png');
    else
        try
            mars_map =imread('http://www.vendian.org/mncharity/dir3/planet_globes/TemporaryURL/mars1_src.jpg','jpg');
        catch
            disp('Issue loading planet surface image from web.');
            mars_map = rand(1,1);
            colormap(handles.colorz);
        end
        imwrite(mars_map,'mars_map.png') 
    end
   map = mars_map;
elseif handles.planet == 5
    % jupiter
    if exist('jupiter_map','file')==2
        jupiter_map = imread('jupiter_map.png');
    else
        try
            jupiter_map =imread('http://webglbasic.com/compiler/textures/solarsystem/jupiter2_1k.jpg','jpg');
        catch
            disp('Issue loading planet surface image from web.');
            jupiter_map = rand(1,1);
            colormap(handles.colorz);
        end
        imwrite(jupiter_map,'jupiter_map.png') 
    end
   map = jupiter_map;
elseif handles.planet == 6
    % saturn
    if exist('saturn_map','file')==2
        saturn_map = imread('saturn_map.png');
    else
        try
            saturn_map =imread('https://s-media-cache-ak0.pinimg.com/736x/75/c1/ec/75c1ec29945eb9d20285c458853fa9bc.jpg','jpg');
        catch
            disp('Issue loading planet surface image from web.');
            saturn_map = rand(1,1);
            colormap(handles.colorz);
        end
        imwrite(saturn_map,'saturn_map.png') 
    end
   map = saturn_map;
elseif handles.planet == 7
    % uranus
    if exist('uranus_map','file')==2
        uranus_map = imread('uranus_map.png');
    else
        try
            uranus_map =imread('http://www.physics.sfasu.edu/astro/courses/phy315/images/uranusmap.jpg','jpg');
        catch
            disp('Issue loading planet surface image from web.');
            uranus_map = rand(1,1);
            colormap(handles.colorz);
        end
        imwrite(uranus_map,'uranus_map.png') 
    end
   map = uranus_map;
elseif handles.planet == 8
    % neptune
    if exist('neptune_map','file')==2
        neptune_map = imread('neptune_map.png');
    else
        try
            neptune_map =imread('http://orbits.wthr.us/img/planets_small/neptune.jpg','jpg');
        catch
            disp('Issue loading planet surface image from web.');
            neptune_map = rand(1,1);
            colormap(handles.colorz);
        end
        imwrite(neptune_map,'neptune_map.png') 
    end
   map = neptune_map;
elseif handles.planet == 9
    %pluto
    if exist('pluto_map','file')==2
        pluto_map = imread('pluto_map.png');
    else
        try
            pluto_map =imread('http://pre15.deviantart.net/5635/th/pre/f/2016/163/7/8/pluto_texture_map___2016_by_fargetanik-da3hjtb.png','jpg');
        catch
            disp('Issue loading planet surface image from web.');
            pluto_map = rand(1,1);
            colormap(handles.colorz);
        end
        imwrite(pluto_map,'pluto_map.png') 
    end
   map = pluto_map;

end

[x,y,z] = sphere(80);
x = x*handles.R;
y = y*handles.R;
z = z*handles.R;
set(handles.planet_surf,'XData',x,'YData',y,'ZData',z,'CData',map);



handles.light_pos = [-1 .3 .26]; % tilt of the earth normalized
% maybe change later for planets

updated_handles=handles;
guidata(hObject,handles)

function handles = Sat_orbit(hObject,handles)

% unpacking handles for readability
a = handles.a;
e = handles.e;
i = deg2rad(handles.i);
Omega = deg2rad(handles.RAAN);
v = 0;
num_P = handles.num_P*20;
dt = handles.dt;
mu = handles.mu;
J2 = handles.J2;
R = handles.R;
w = deg2rad(handles.w);

% ORBITAL VARIABLES CALCULATIONS

% Mean motion
n = sqrt(mu/a^3);

% Eccentric anomoly
E = 2*atan( sqrt((1-e)/(1+e)) * tan(v/2) );

% Initial mean anomoly
Mo = E - e*sin(E);

% period
handles.P = 2*pi/n;

strg = strcat(num2str(round(handles.P/60)),'min');
set(handles.p_output,'String',strg)


% time vector
tf = round(handles.P*num_P);
t = 0:dt:tf ;

% RAAN
dot_omega = 3/2*(1-e^2)^2*n*J2*(R/a)^2*cos(i);

% initialize time dependent vectors before theyre assigned in loop
M_t = zeros(size(t));
v_t = zeros(size(t));
Omega_t = zeros(size(t));



for j=1:length(t)
    
   
   % Mean anomoly [rad]
   M_t(j) = Mo + n*t(j);
   
   % True anomoly [rad]
   v_t(j) = M_t(j) + (2*e - e^3/4+e^5*e^5/96).*sin(M_t(j));
   
   % time dependent RAAN
   Omega_t(j) = Omega+(dot_omega*t(j));
   
end

% orbital radius vector
r = a.*(1-e^2)./(1+e.*cos(v_t));

%vel = sqrt( mu*( 2./r - 1/a ));

% rotation matrix divided into rows
first = cos(Omega_t).*cos(w+v_t) - sin(Omega_t).*sin(w+v_t).*cos(i);
second = sin(Omega_t).*cos(w+v_t) + cos(Omega_t).*sin(w+v_t).*cos(i);
third = sin(w+v_t).*sin(i);


% conversion to latitude longitude from rotation matrix
lat = rad2deg(asin(third));
longi = rad2deg( acos(first./cos(deg2rad(lat))) );
lon = longi .* sign(second);

alt = r-handles.R;
planet = referenceEllipsoid(handles.name,'m');
[ox,oy,oz] = geodetic2ecef(planet,lat,lon,alt);

% we don't care about imaginary part

handles.ox = real(ox);
handles.oy = real(oy);
handles.oz = real(oz);


% update mltan output
handles = get_MLTAN(hObject,handles);
strn = strcat( handles.eqtime,'AM-',handles.eqtime,'PM');
set(handles.eqtime_output,'String',strn)


guidata(hObject,handles)

function handles = yearly_sat_orbit(hObject,handles)
a = handles.a;
e = handles.e;
i = deg2rad(handles.i);
Omega = deg2rad(handles.RAAN);
J2 = handles.J2;
R = handles.R;
mu = handles.mu;

year_t = 0:3600*24:3600*24*365*20;


% Mean motion
n = sqrt(mu/a^3);

% RAAN
dot_omega = 3/2*(1-e^2)^2*n*J2*(R/a)^2*cos(i);

% initialize time dependent vectors before theyre assigned in loop
handles.Omega_t = zeros(size(year_t));



for j=1:length(year_t)
   % time dependent RAAN
   handles.Omega_t(j) = Omega+(dot_omega*year_t(j));
   
end
guidata(hObject,handles)

function handles = draw_orbital_plane(hObject,handles,first)

xfill = handles.ox(1:handles.tail_length);
yfill = handles.oy(1:handles.tail_length);
zfill = handles.oz(1:handles.tail_length);

if first == 1
    plane_color = [1 .6 0];
    hold on
    handles.oplane = fill3(handles.earth_axes,xfill, yfill, zfill,...
        plane_color,'EdgeColor','none','FaceAlpha',.5);
else
    set(handles.oplane,'XData',xfill,'YData',yfill,'ZData',zfill);
end

guidata(hObject,handles)

function handles = calculate_t(hObject,handles)
%time has to be recalculated for each planet selection

% time 
handles.t = 0:handles.dt:round(20*handles.num_P*handles.P);


% tail_length is one orbit length of indeces 
handles.tail_length = round(handles.P/handles.dt);

% updates handles
guidata(hObject,handles)

function handles = get_MLTAN(hObject,handles)
%given an orbit, figures out the local mean solar time of the ascending node
%(which is the crossing time of the equatorial plane)

% noon vector components
u = get(handles.earth_quiv,'UData');
v = get(handles.earth_quiv,'VData');
n = [u v 0];

% ascending node components assuming sat starts at AN
x = handles.ox(1);
y = handles.oy(1);
an = [x y 0];


angle = atan2(norm(cross(an,n)),dot(an,n))*180/pi;

% convert that angle to time
% 6 hours corresponds to 90 degrees i think
first_hr = num2str(round(angle*6/90));

handles.eqtime = first_hr;

%disp(strcat( first_hr,'AM-',first_hr,'PM'))



% updates handles
guidata(hObject,handles)

function handles = draw_noon_vectors(hObject, handles,first)

x_light = handles.light_pos(1)*handles.R*1.5;
y_light = handles.light_pos(2)*handles.R*1.5;
z_light = handles.light_pos(3)*handles.R*1.5;
if first ==1
    handles.earth_quiv = quiver3(handles.earth_axes,0,0,0,x_light,y_light,...
        z_light,'Color',[0.8500 0.3250 0.0980]); % for right axis
else
    set(handles.earth_quiv,'UData',x_light,'VData',y_light','WData',z_light)
end






