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
%       Sun synchronous GUI demonstrates sun synchronous orbits. In manual
%       mode, the user made edit any of the orbit parameters and the orbit
%       will be displayed. In auto mode, when the user adjusts the orbit
%       parameters the GUI will automatically adjust the orbit to be sun
%       synchronous
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

set(gcf,'Name','Sun Synchronous Orbit Design','CloseRequestFcn',@clothes);

%%%%%%%%%%%%%%%%%%%%%%%%%  Earth (Right) Axis  %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input from Button
handles.planet = get(handles.planet_selection,'Value');
handles.step = 200;
handles.earthdays=0;
handles.righttime=0;
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
handles.lighty = light(handles.earth_axes,'Position',handles.light_pos,...
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Sun and Earth Axis %%%%%%%%%%%%%%%%%%%%%%%%%%

% Setting up Sun/Earth axis
axes(handles.solar_axes)
axis(handles.solar_axes,'vis3d','off')

rotate3d on;
axis equal


[x_sphere,y_sphere,z_sphere] = sphere();

% the sun
sun_factor = 10;
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
set(handles.earth,'FaceColor',handles.colorz,'EdgeColor','none');

% earth's orbit
theta = 0:pi/20:2*pi;
r_o = r_earth+sun_factor-1.5;
plot(handles.solar_axes,r_o*cos(theta),r_o*sin(theta),'--')

q = 1.8;
set(handles.solar_axes,'XLim', [-r_earth*q r_earth*q],...
    'YLim', [-r_earth*q r_earth*q], 'ZLim',[-r_earth*q r_earth*q])

% setting up rotation transform for earth around sun
t2 = hgtransform(handles.solar_axes);
set(handles.earth,'Parent',t2);

% noon vector on right axis
handles = draw_noon_vectorR(hObject,handles,1);

% noon vector left axis
handles.solar_quiv = quiver(handles.solar_axes,r_earth,r_earth,-5,-5,...
    'Color',[0.8500 0.330 0.10],'MaxHeadSize',.8); 

handles.solar_orb = quiver(handles.solar_axes,r_earth,r_earth,-30,-30,...
    'Color',[0 0 0],'MaxHeadSize',1); 

handles = yearly_sat_orbit(hObject,handles);


set(handles.solar_quiv,'Parent',t2);
set(handles.solar_orb,'Parent',t2);


% calculates lat and lon of sat orbit
handles.num_P=20;
handles.dt=20;

handles = Sat_orbit(hObject, handles);

  
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

% draw equator on planet
handles = plot_equator_poles(hObject,handles,1);

% view top down left axis
view(handles.solar_axes,2)

%set default mode
handles.mode=0;
%j2 output
set(handles.J2readout,'String',strcat('Using J2 = ',num2str(handles.J2)));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Graph  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  

handles = plot_lat_MLTAN(hObject,handles,1);

% x y labels
set(get(handles.graph, 'xlabel'),'String','Mean Local Time [hr]');
set(get(handles.graph, 'ylabel'),'String','Lattitude [deg]');
% x y limits hard set
xlim(handles.graph,[0 24]);
ylim(handles.graph,[-90 90]);
% x y ticks
set(handles.graph,'Xtick',[0 2 4 6 8 10 12 14 16 18 20 22 24])
% this line below only works for 2017
% xticks(handles.graph,[0 2 4 6 8 10 12 14 16 18 20 22 24])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Animation Loop  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% k indexes time
handles.k=0;
% save handles
guidata(hObject,handles)
while ishandle(hObject) && handles.k <=length(handles.t)
    
    handles = guidata(hObject);
    handles.k = handles.k+1;
    
    %%%%%% Earth (Right) Axis %%%%
    % make earth rotate
    n = handles.rot*(handles.t(handles.k));
    Txy = makehgtform('zrotate',n);
    set(t1,'Matrix',Txy);
    
    %leaving this bit in here in case i decide to plot a tail to the sat
    %(just add start:stop)
    % if beginning
%     if handles.k<handles.tail_length 
%         start=1;
%     % else if middle
%     else
%         start=start+1;
%     end
    stop=handles.k;
    
    % plot satellite
    set(handles.hsat,'XData',handles.ox(stop),'YData',handles.oy(stop),...
        'ZData',handles.oz(stop));
    
    % only do this once every orbit
    
    handles = draw_orbital_plane(hObject,handles,0);
   
    
    %%%%%% Sun/ Earth (Left) Axis %%%%%
    n2 = handles.yr(handles.k);
    Txy2 = makehgtform('zrotate',n2);
    set(t2,'Matrix',Txy2);
    axis(handles.solar_axes,'off')
    
    arrowsize=10;
    % the pi/4 offset is because the noon vector does not line up with the
    % axis. noon vector is [-5 -5]
%     angle = (-handles.Omega_t(handles.k)+pi/4)*180/pi;
%     newu1 = -arrowsize*sin(-handles.Omega_t(handles.k)+pi/4);
%     newv1 = -arrowsize*cos(-handles.Omega_t(handles.k)+pi/4);
    
%     earthpos = [28.5*cos(n2+pi/4) 28.5*sin(n2+pi/4)];
% % %     plot(handles.solar_axes,earthpos(1),earthpos(2),'.','MarkerSize',40)
    

    newu1=handles.year_x(stop);
    newv1=handles.year_y(stop);

    sun_angle = handles.tilt/180*pi;
    rotm = [cos(sun_angle) 0 sin(sun_angle); 0 1 0; -sin(sun_angle) 0 cos(sun_angle)];
    answer = rotm*[newu1;newv1;0];
    magn = sqrt(answer(1)^2+answer(2)^2);
    newu = answer(1)/magn*arrowsize;
    newv = answer(2)/magn*arrowsize;


    
    set(handles.solar_orb,'UData',newu,'VData',newv);
    
    handles.earthdays = handles.earthdays+handles.P_year/handles.step;
    
    if handles.earthdays <= 365
        str = ['Earth Days Passed ',num2str(floor(handles.earthdays))];
    else
        str = ['Earth Years Passed ',num2str(round(handles.earthdays/365,1))];
    end
    set(handles.left_day,'String',str)
    
    
    % add a string for right axis
    handles.righttime = handles.righttime+handles.P/handles.step;
    if handles.righttime < 3600
        str2 = ['Minutes Passed ',num2str(round(handles.righttime/60))];
    else
        str2 = ['Hours Passed ',num2str(round(handles.righttime/3600,1))];
    end
    set(handles.righttime_text,'String',str2)
    
    
    %output
    
    % if crossing equator, update MLTAN
    if handles.oz(stop) < 1E5 && handles.oz(stop) > -1E5
        eqtime = round(get_MLT(hObject,handles,stop));
        if eqtime <=12
            secondhr = eqtime+12;
        else
            secondhr = eqtime-12;
        end
        str = strcat( num2str(eqtime),'-',num2str(secondhr), '[hr]');
        set(handles.eqtime_output,'String',str)
    end
    
    % restart loop at the end
    % at the close of figure, i set handles.k = length(handles.t)*2
    if handles.k>=length(handles.yr)-1 && handles.k~=length(handles.t)*2
        handles.k = 0;
    end
    
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
old_i = handles.i;

if isinf(i) || isnan(i)
    set(handles.i_edit,'String',num2str(old_i)) % put old i back in input box
    set(handles.errormsg,'String','Inclination input not valid')
% if auto mode
elseif handles.mode==1
    % if successful input
    set(handles.errormsg,'String','') % reset error message
    % if change inclination in auto mode, change a
    e = handles.e;
    i = deg2rad(i);
    J2 = handles.J2;
    R = handles.R;
    mu = handles.mu;
    dot_omega_planet = 2*pi/(handles.P_year*24*3600);
    
    % from SME SMAD equation for nodal precessoin
    a = ((-3/2*sqrt(mu)*J2*R^2*cos(i)*(1-e^2)^(-2))/dot_omega_planet)^(2/7);
    alt = (a-R)/1000;
    rmin = a*(1-e);
    
    
    if rmin<=handles.R || isreal(a)==0 % if alt doesn't work
        % throw error message
        set(handles.errormsg,'String','Inclination not possible')
        set(handles.i_edit,'String',num2str(old_i)) % put old i back in input box
    else % if alt works
        handles.a=a;
        handles.i=rad2deg(i);
        set(handles.a_edit,'String',alt);
        handles = something_changed(hObject,handles);
    end
else % if not auto mode, successful input
    handles.i=i;
    handles = something_changed(hObject,handles);
    set(handles.errormsg,'String','') % reset error message
end


guidata(hObject,handles)

function a_edit_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to scanrate_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of scanrate_edit as text
%        str2double(get(hObject,'String')) returns contents of scanrate_edit as a double
% 
% This function is called when the minimum altitude field is edited.
% This function checks the validity of the input and updates the corresponding satellite
% orbit variable, semi-major axis, in the handles structure.
% Orbit is recalcaulted and displayed.


alty = str2double(get(hObject,'String'));
old_alty = (handles.a*(1-handles.e)-handles.R)/1e3;

%check for valid input
if isinf(alty) || isnan(alty)
    set(handles.a_edit,'String',num2str(old_alty)) % put old alt back in box
    set(handles.errormsg,'String','Altitude input not valid') %error msg
else
    % semi major axis minimum if eccentric
    rmin = (alty*1E3+handles.R);
    a = rmin/(1-handles.e);
    
    %check for valid semi major axis
    if rmin<=handles.R 
        set(handles.a_edit,'String',num2str(old_alty)) % put old alt back in input box
        set(handles.errormsg,'String','Altitude not possible') % error message
    % if auto mode
    elseif handles.mode==1
%         % if change alt in auto mode, change i
        e = handles.e;
        J2 = handles.J2;
        R = handles.R;
        mu = handles.mu;
        dot_omega_planet = 2*pi/(handles.P_year*24*3600);

        expression = dot_omega_planet/(-3/2*sqrt(mu/a^3)*J2*(R/a)^2*(1-e^2)^(-2));
        radi = acos(expression);
        i = rad2deg(radi);

        if isreal(i)==0 % if i is imaginary
            set(handles.a_edit,'String',num2str(old_alty)) % put old alt back in input box
            set(handles.errormsg,'String','Altitude not possible') % error msg
        else
            set(handles.errormsg,'String','') % reset error message
            handles.i=i;
            handles.a=a;
            set(handles.i_edit,'String',i);
            handles = something_changed(hObject,handles);
        end
%         end

    % if manual mode
    else
        set(handles.errormsg,'String','') % reset error message
        handles.a = a;
    end
end

handles = something_changed(hObject,handles);
guidata(hObject,handles)

function e_edit_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to scanrate_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of scanrate_edit as text
%        str2double(get(hObject,'String')) returns contents of scanrate_edit as a double
e = str2double(get(hObject,'String'));
old_e =handles.e;

if isinf(e) || isnan(e)
    set(handles.e_edit,'String',num2str(old_e)) % put old e back in input box
    set(handles.errormsg,'String','Eccentricity input not valid')
elseif e>1 || e<0
    set(handles.e_edit,'String',num2str(old_e)) % put old e back in input box
    set(handles.errormsg,'String','Eccentricity input must be >0 and <1')
elseif (handles.a*(1-e)) <= handles.R
    set(handles.e_edit,'String',num2str(old_e)) % put old e back in input box
    set(handles.errormsg,'String','Eccentricity input not possible')
else
    set(handles.errormsg,'String','') % reset error message
    % if auto mode
    if handles.mode==1
        % if change inclination in auto mode, change a
        i = deg2rad(handles.i);
        J2 = handles.J2;
        R = handles.R;
        mu = handles.mu;
        dot_omega_planet = 2*pi/(handles.P_year*24*3600);
        original_a = handles.a;

        % from SME SMAD equation for nodal precessoin
        a = ((-3/2*sqrt(mu)*J2*R^2*cos(i)*(1-e^2)^(-2))/dot_omega_planet)^(2/7);
        alt = (a*(1-e)-R)/1000;

        if alt<0 || isreal(a)==0% if alt doesn't work
            % throw error message
            set(handles.e_edit,'String',num2str(old_e)) % put old e back in input box
            set(handles.errormsg,'String','Eccentricity not possible')
        else % if alt works
            handles.a=a;
            handles.e=e;
            set(handles.a_edit,'String',alt);
            handles = something_changed(hObject,handles);

           if handles.r_vector < handles.R
               set(handles.e_edit,'String',num2str(old_e)) % put old e back in input box
               set(handles.errormsg,'String','Eccentricity not possible')
               handles.a = original_a;
               handles = something_changed(hObject,handles);
           end

        end
    else % if manual mode
        handles.e=e;
    end
end


handles = something_changed(hObject,handles);
guidata(hObject,handles)

function w_edit_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to scanrate_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of scanrate_edit as text
%        str2double(get(hObject,'String')) returns contents of scanrate_edit as a double
w = str2double(get(hObject,'String'));
handles.w = w;
handles = something_changed(hObject,handles);
guidata(hObject,handles)

function RAAN_edit_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to scanrate_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of scanrate_edit as text
%        str2double(get(hObject,'String')) returns contents of scanrate_edit as a double
RAAN = str2double(get(hObject,'String'));
handles.RAAN = RAAN;
handles = something_changed(hObject,handles);
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
handles.earthdays = 0;
handles.righttime=0;
handles = Update_Planet_Constants(hObject,handles);
handles = Draw_Planet(hObject,handles);

% update orbit of satellite to account for change in radius of planet
handles.a = str2double(get(handles.a_edit,'String'))*1000+handles.R;
% updates sat orbit, lat/mltan graph...
handles = something_changed(hObject,handles);
% set face color left axis
set(handles.earth,'FaceColor',handles.colorz);
% update time
handles = calculate_t(hObject,handles);
% noon vector
handles = draw_noon_vectorR(hObject,handles,0);
% draw equator and poles on planet
handles = plot_equator_poles(hObject,handles,0);

if handles.mode==1
    if handles.J2 <= 0.000027
        set(handles.impossible_edit,'String','Sun synchronous orbit is not possible for this planet')
    else
        set(handles.impossible_edit','String','')
    end
else
    set(handles.impossible_edit','String','')
end

% updates handles
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
    handles.yr = 0:2*pi/handles.step:2*pi*20;
    handles.P_year = 365;
    handles.tilt = 23.5/180*pi;
elseif planet == 1
    % Mercury
    handles.mu = 2.2032e13;
    handles.J2 = 0.00006;
    handles.R = 2.439E6;
    handles.rot = 1.24e-06;
    handles.colorz = [.6 .6 .6];
    handles.name = 'mercury';
    handles.yr = 0:2*pi/handles.step:2*pi*20;
    handles.P_year = 88;
    handles.tilt = 0.03/180*pi;
elseif planet == 2
    % Venus
    handles.mu = 3.24859E14;
    handles.J2 = 0.000027;
    handles.R = 6.051E6;
    handles.rot = 2.99E-7;
    handles.colorz = [250/255, 220/255, 121/255];
    handles.name = 'venus';
    handles.yr = 0:-2*pi/handles.step:-2*pi*20;
    handles.P_year = 224.7;
    handles.tilt = 2.64/180*pi;
elseif planet == 4
    % Mars
    handles.mu = 4.282837e13;
    handles.J2 = 0.001964;
    handles.R = 3.396E6;
    handles.rot = 7.088e-05;
    handles.colorz = [.9 .3 0];
    handles.name = 'mars';
    handles.yr = 0:2*pi/handles.step:2*pi*20;
    handles.P_year = 686.93;
    handles.tilt = 25.19/180*pi;
elseif planet == 5
    % Jupiter
    handles.mu = 1.26686534E17;
    handles.J2 = 0.01475;
    handles.R = 71.492E6;
    handles.rot = 1.77E-4;
    handles.colorz = [210/255, 200/255, 150/255];
    handles.name = 'jupiter';
    handles.yr = 0:2*pi/handles.step:2*pi*20;
    handles.P_year = 11.86*365;
    handles.tilt = 3.13/180*pi;
elseif planet == 6
    % Saturn
    handles.mu = 3.7931187E16;
    handles.J2 = 0.01645;
    handles.R = 60.268E6;
    handles.rot =1.63E-4;
    handles.map = [210/255, 200/255, 150/255];
    handles.name = 'saturn';
    handles.yr = 0:2*pi/handles.step:2*pi*20;
    handles.P_year = 10755;
    handles.tilt = 26.73/180*pi;
elseif planet == 7
    % Uranus
    handles.mu = 5.794e15;
    handles.J2 = 0.012;
    handles.R = 25.559E6;
    handles.rot = -1.04E-4;
    handles.colorz = [209/255 231/255 231/255];
    handles.name = 'uranus';
    handles.yr = 0:2*pi/handles.step:2*pi*20;
    handles.P_year = 84*365;
    handles.tilt = 82.23/180*pi;
elseif planet ==8
    % Neptune
    handles.mu = 6.809e15;
    handles.J2 = 0.0004;
    handles.R = 24.764E6;
    handles.rot = 1.08E-4;
    handles.colorz = [63/255 84/255 186/255];
    handles.name = 'neptune';
    handles.yr = 0:2*pi/handles.step:2*pi*20;
    handles.P_year = 164*365;
    handles.tilt = 28.32/180*pi;
else
    %Pluto
    handles.mu = 8.71E11;
    handles.J2 = 0;
    handles.R = 1.195E6;
    handles.rot = -1.29E-5;
    handles.colorz = [.8 .7 .5];
    handles.name = 'pluto';
    handles.yr = 0:2*pi/handles.step:2*pi*20;
    handles.P_year = 248*365;
    handles.tilt = 57.47/180*pi;
end

% change default altitude based on planet
handles.a = handles.R*1.5;
set(handles.a_edit,'String',round((handles.a-handles.R)/1E3));


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

magn = sqrt(1+tan(handles.tilt)^2)*.8;
handles.light_pos = [-1/magn 0 tan(handles.tilt)/magn]; % tilt of the planet normalized
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
num_P = handles.num_P;
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

handles.ox = r.*first;
handles.oy = r.*second;
handles.oz = r.*third;
handles.r_vector = r;
% hopefully imaginary numbers aren't a problem^



% update handles w/ lat lon so that it can be graphed in outputs
handles.lat = lat;
handles.lon = lon;



guidata(hObject,handles)

function handles = yearly_sat_orbit(hObject,handles)
a = handles.a;
e = handles.e;
i = deg2rad(handles.i);
Omega = deg2rad(handles.RAAN)+225*pi/180;
J2 = handles.J2;
R = handles.R;
mu = handles.mu;
w = deg2rad(handles.w);

% year_t = 0:3600*24:3600*24*365*20;


% Mean motion
n = sqrt(mu/a^3);

P = 2*pi/n;

year_t=0:P:3600*24*365*20;

v=0;
% Eccentric anomoly
E = 2*atan( sqrt((1-e)/(1+e)) * tan(v/2) );

% Initial mean anomoly
Mo = E - e*sin(E);

% RAAN
dot_omega = 3/2*(1-e^2)^2*n*J2*(R/a)^2*cos(i);

% earth's rotation rad/s
% need to change by planet
dot_omega_earth = 2*pi/(handles.P_year*24*3600);



% initialize time dependent vectors before theyre assigned in loop
%handles.Omega_t = zeros(size(year_t));
Omega_t = zeros(size(year_t));
M_t = zeros(size(year_t));
v_t = zeros(size(year_t));



for j=1:length(year_t)
    % Mean anomoly [rad]
   M_t(j) = Mo + n*year_t(j);
   
   % True anomoly [rad]
   v_t(j) = M_t(j) + (2*e - e^3/4+e^5*e^5/96).*sin(M_t(j));
   
   % time dependent RAAN
   %handles.Omega_t(j) = Omega+(dot_omega*year_t(j))+dot_omega_earth*year_t(j);
  Omega_t(j) = Omega+(dot_omega*year_t(j))+dot_omega_earth*year_t(j);
%   Omega_t(j) = Omega+(dot_omega*year_t(j));
   
end

% orbital radius vector
r = a.*(1-e^2)./(1+e.*cos(v_t));

% rotation matrix divided into rows
first = cos(Omega_t).*cos(w+v_t) - sin(Omega_t).*sin(w+v_t).*cos(i);
second = sin(Omega_t).*cos(w+v_t) + cos(Omega_t).*sin(w+v_t).*cos(i);
third = sin(w+v_t).*sin(i);

handles.year_x = r.*first;
handles.year_y = r.*second;
handles.year_z = r.*third;


guidata(hObject,handles)

function handles = draw_orbital_plane(hObject,handles,first)
% on right axis
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

handles.t = 0:handles.dt:round(handles.num_P*handles.P);


% tail_length is one orbit length of indeces 
handles.tail_length = round(handles.P/handles.dt);

% updates handles
guidata(hObject,handles)

function eqtime = get_MLT(hObject,handles,index)
% given an orbit, figures out the local mean solar time of the ascending node
% (which is the crossing time of the equatorial plane)
% im not using handles anywhere

% noon vector components
u = get(handles.earth_quiv,'UData');
v = get(handles.earth_quiv,'VData');
n = [u v 0];

% ascending node components assuming sat starts at AN
x = handles.ox(index);
y = handles.oy(index);
an = [x y 0];
iserror=0;

try
    handles.hourangle = atan2(norm(cross(an,n)),dot(an,n)) *180/pi;
catch
    iserror=1;
    set(handles.errormsg,'String','Problem with hour angle') % reset error message
end

if iserror==0
    % convert that angle to time
    % 6 hours corresponds to 90 degrees i think
    first_hr = handles.hourangle*6/90;

    if y<handles.light_pos(2)*handles.R*1.5
       first_hr = 24-first_hr;
    end


    eqtime = first_hr;
else
    eqtime=0;
end



guidata(hObject,handles)


function handles = draw_noon_vectorR(hObject, handles,first)
% draws noon vector on right axis
x_light = handles.light_pos(1)*handles.R*1.5;
y_light = handles.light_pos(2)*handles.R*1.5;
z_light = handles.light_pos(3)*handles.R*1.5;
if first ==1
    handles.earth_quiv = quiver3(handles.earth_axes,0,0,0,x_light,y_light,...
        z_light,'Color',[0.8500 0.330 0.10]); % for right axis
else
    set(handles.earth_quiv,'UData',x_light,'VData',y_light','WData',z_light)
end

guidata(hObject,handles)

function handles = plot_equator_poles(hObject,handles,first)
% plots equator on right axis
long = -180:10:180;
lat = zeros(size(long));
alt = zeros(size(long));

planet = referenceEllipsoid(handles.name,'m');
[x,y,z] = geodetic2ecef(planet,lat,long,alt);

zpole = -handles.R*1.25:handles.R/10:handles.R*1.25;
xpole = zeros(size(zpole));
ypole = zeros(size(zpole));

hold on

if first ==1
    handles.equator = plot3(handles.earth_axes,x,y,z,'--','Color','m');
    handles.pole = plot3(handles.earth_axes,xpole,ypole,zpole,'Color','b');
else
    set(handles.equator,'XData',x,'ydata',y,'zdata',z);
    set(handles.pole,'XData',xpole,'ydata',ypole,'zdata',zpole)
end



guidata(hObject,handles)



function handles = plot_lat_MLTAN(hObject,handles,first)
% plots or updates the plot of lat/lon crossing time graph

num_pts_per_P = round(handles.P/handles.dt);
sizevar = num_pts_per_P;
% pre allocate
cross_time = zeros(sizevar);
latties = zeros(sizevar);
count=0;
step=5;
% grab every 5 latitude and longitudes 
for i=1:step:sizevar
%     if handles.lat(i) < .1
        count=count+1;
        eqtime = get_MLT(hObject,handles,i);
        cross_time(count)= eqtime;
        latties(count) = handles.lat(i);
%     end
end
indeces = cross_time~=0;
cross_time = cross_time(indeces);
latties = latties(indeces);

if first==1
    handles.linegraph = plot(handles.graph,cross_time,latties,'*');
    
else
    
    set(handles.linegraph,'XData',cross_time,'YData',latties);
end



guidata(hObject,handles)

function handles = something_changed(hObject,handles)
% this function is called when something has changed

handles = Sat_orbit(hObject,handles);

handles = plot_lat_MLTAN(hObject,handles,0);

handles = calculate_t(hObject,handles);

handles = yearly_sat_orbit(hObject,handles);

handles.earthdays =0;

handles.righttime=0;

handles.k=0; %restart loop

set(handles.J2readout,'String',strcat('Using J2 = ',num2str(handles.J2))); %display J2



%

% change light
magn = sqrt(1+tan(handles.tilt)^2)*.8;
handles.light_pos = [-1/magn 0 tan(handles.tilt)/magn]; % tilt of the planet normalized
set(handles.lighty,'Position',handles.light_pos);

guidata(hObject,handles)

function handles = auto_callback(hObject, ~, handles)

handles.mode = 1;

i = str2double(get(handles.i_edit,'String'));

if isinf(i) || isnan(i)
    set(handles.errormsg,'String','Inclination input not valid')
% if auto mode
elseif handles.mode==1
    % if successful input
    set(handles.errormsg,'String','') % reset error message
    % if change inclination in auto mode, change a
    e = handles.e;
    i = deg2rad(i);
    J2 = handles.J2;
    R = handles.R;
    mu = handles.mu;
    dot_omega_planet = 2*pi/(handles.P_year*24*3600);
    
    % from SME SMAD equation for nodal precessoin
    a = ((-3/2*sqrt(mu)*J2*R^2*cos(i)*(1-e^2)^(-2))/dot_omega_planet)^(2/7);
    alt = (a-R)/1000;
    rmin = a*(1-e);
    
    
    if rmin<=handles.R || isreal(a)==0 % if alt doesn't work
        % throw error message
        set(handles.errormsg,'String','Inclination not possible')
    else % if alt works
        handles.a=a;
        handles.i=rad2deg(i);
        set(handles.a_edit,'String',alt);
        handles = something_changed(hObject,handles);
    end
else % if not auto mode, successful input
    handles.i=i;
    handles = something_changed(hObject,handles);
    set(handles.errormsg,'String','') % reset error message
end


guidata(hObject,handles)

function handles = manual_callback(hObject, ~, handles)

handles.mode = 0;
set(handles.impossible_edit','String','')

guidata(hObject,handles)

function clothes(src,callbackdata)
%handles.k = length(handles.t)*2;
% delete(hObject);

delete(gcf)



