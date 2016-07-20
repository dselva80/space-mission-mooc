function Keplerian_Orbits
% Function KEPLERIAN_ORBITS allows for user alteration of several
% orbital parameters and animates the orbit(s) when selected by the user.
%
% Written by Patrick Voorhees 20-Jul-2016 10:49:00
%       Copyright (c) 2016 Patrick Voorhees (pwv9@cornell.edu)
% Advisors: Professors Dmitry Savransky and Daniel Selva (Cornell 
%       University, Department of Mechanical and Aerospace Engineering)
% Purpose: create a visualization of the effects of orbital parameters on
%       the shape, size, and orientation of an orbit, as well as on the 
%       velocity, position, and period of an orbiting body 
%       (MAE 4060: Introduction to Spaceflight Mechanics). 
% NOTE: two orbits may be selected to demonstrate phenomena such as
%       resonance, however each orbiting body does not effect the other, 
%       this is solved as a strictly two body problem.

% set options for figure display, open a figure and set the axes
ops = {'BackgroundColor','ForegroundColor','FontWeight','FontSize'};
H.fig = figure('Name','Keplerian_Orbits','Position',[320 60 900 600],'Color','k','Resize','off');
H.ax1 = axes('Units','Pixels','Position',[0 50 600 600],...
    'xcolor','none','ycolor','none','zcolor','none','Color','none');
uistack(H.ax1,'bottom');
H.text14 = uicontrol('Style','text','Position',[20 10 500 30],'String','',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10,'HorizontalAlignment','left');
H.text15 = uicontrol('Style','text','Position',[20 50 500 30],'String','',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10,'HorizontalAlignment','left');

% set first orbit options
H.text2 = uicontrol('Style','text','Position',[625 575 100 20],'String','FIRST ORBIT',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.text3 = uicontrol('Style','text','Position',[600 550 150 20],'String','a1 = [.723,1.523] au',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.slider1 = uicontrol('Style','slider','Position',[615 520 120 30],ops{1},'k',ops{2},'w',ops{3},...
    'Bold',ops{4},10,'Min',.723,'Max',1.523,'Value',.723);
H.edit1 = uicontrol('Style','edit','Position',[625 490 100 25],'String','.723',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.text4 = uicontrol('Style','text','Position',[600 460 150 20],'String','e1 = [0,.9]',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.slider2 = uicontrol('Style','slider','Position',[615 425 120 30],ops{1},'k',ops{2},'w',ops{3},...
    'Bold',ops{4},10,'Min',0,'Max',.9,'Value',0,'SliderStep',[1/30 1/9]);
H.edit2 = uicontrol('Style','edit','Position',[625 395 100 25],'String','0',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.text5 = uicontrol('Style','text','Position',[600 360 150 20],'String','i1 = [0,360]',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.slider3 = uicontrol('Style','slider','Position',[615 330 120 30],ops{1},'k',ops{2},'w',ops{3},...
    'Bold',ops{4},10,'Min',0,'Max',360,'Value',0,'SliderStep',[1/72 1/12]);
H.edit3 = uicontrol('Style','edit','Position',[625 300 100 25],'String','0',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.text6 = uicontrol('Style','text','Position',[600 265 150 20],'String','RAAN1 = [0,360]',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.slider4 = uicontrol('Style','slider','Position',[615 235 120 30],ops{1},'k',ops{2},'w',ops{3},...
    'Bold',ops{4},10,'Min',0,'Max',360,'Value',0,'SliderStep',[1/72 1/12]);
H.edit4 = uicontrol('Style','edit','Position',[625 205 100 25],'String','0',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.text7 = uicontrol('Style','text','Position',[600 175 150 20],'String','w1 = [0,360]',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.slider5 = uicontrol('Style','slider','Position',[615 145 120 30],ops{1},'k',ops{2},'w',ops{3},...
    'Bold',ops{4},10,'Min',0,'Max',360,'Value',0,'SliderStep',[1/72 1/12]);
H.edit5 = uicontrol('Style','edit','Position',[625 115 100 25],'String','0',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);

% set second orbit options
H.text8 = uicontrol('Style','text','Position',[775 575 100 20],'String','SECOND ORBIT',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.text9 = uicontrol('Style','text','Position',[750 550 150 20],'String','a2 = [.723,1.523] au',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.slider6 = uicontrol('Style','slider','Position',[765 520 120 30],ops{1},'k',ops{2},'w',ops{3},...
    'Bold',ops{4},10,'Min',.723,'Max',1.523,'Value',.723);
H.edit6 = uicontrol('Style','edit','Position',[775 490 100 25],'String','.723',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.text10 = uicontrol('Style','text','Position',[750 460 150 20],'String','e2 = [0,.9]',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.slider7 = uicontrol('Style','slider','Position',[765 425 120 30],ops{1},'k',ops{2},'w',ops{3},...
    'Bold',ops{4},10,'Min',0,'Max',.9,'Value',0,'SliderStep',[1/30 1/9]);
H.edit7 = uicontrol('Style','edit','Position',[775 395 100 25],'String','0',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.text11 = uicontrol('Style','text','Position',[750 360 150 20],'String','i2 = [0,360]',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.slider8 = uicontrol('Style','slider','Position',[765 330 120 30],ops{1},'k',ops{2},'w',ops{3},...
    'Bold',ops{4},10,'Min',0,'Max',360,'Value',0,'SliderStep',[1/72 1/12]);
H.edit8 = uicontrol('Style','edit','Position',[775 300 100 25],'String','0',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.text12 = uicontrol('Style','text','Position',[750 265 150 20],'String','RAAN2 = [0,360]',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.slider9 = uicontrol('Style','slider','Position',[765 235 120 30],ops{1},'k',ops{2},'w',ops{3},...
    'Bold',ops{4},10,'Min',0,'Max',360,'Value',0,'SliderStep',[1/72 1/12]);
H.edit9 = uicontrol('Style','edit','Position',[775 205 100 25],'String','0',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.text13 = uicontrol('Style','text','Position',[750 175 150 20],'String','w2 = [0,360]',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.slider10 = uicontrol('Style','slider','Position',[765 145 120 30],ops{1},'k',ops{2},'w',ops{3},...
    'Bold',ops{4},10,'Min',0,'Max',360,'Value',0,'SliderStep',[1/72 1/12]);
H.edit10 = uicontrol('Style','edit','Position',[775 115 100 25],'String','0',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);

% set animation preferences
H.orbits = uicontrol('Style','checkbox','Position',[500 60 100 25],'String','2 Orbits?',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.anim = uicontrol('Style','togglebutton','Position',[500 20 100 25],'String','Play',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.text16 = uicontrol('Style','text','Position',[612 60 126 20],'String','Manual Commands',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.pref = uicontrol('Style','popupmenu','Position',[625 20 100 25],'String',{'Rotate on','Zoom on'},...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.text17 = uicontrol('Style','text','Position',[765 60 120 20],'String','Animation Speed',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.speed = uicontrol('Style','slider','Position',[775 20 100 25],ops{1},'k',ops{2},'w',ops{3},...
    'Bold',ops{4},10,'Min',.25,'Max',3,'Value',.25);

% set background image to simulate space
Nstars = 400;
xpos = 2.9*rand(1,Nstars)-1.2;
ypos = 2.9*rand(1,Nstars)-1.2;
siz = .011*rand(1,Nstars);
th = linspace(0,2*pi,50);
hold on;
for i=1:Nstars
    % check if background stars are too small to be visible, resize, plot
    if siz(i) <= .001
        siz(i) = 10*siz(i);
    end
    fill(siz(i)*cos(th)+xpos(i),siz(i)*sin(th)+ypos(i),'w');
end
hold off;

% set second set of axes for the plotting over background image
H.ax2 = axes('Units','Pixels','Position',[0 50 600 600],...
    'xcolor','none','ycolor','none','zcolor','none','Color','none');
uistack(H.ax2,'top');
H.ax2.View = [-100,20];

% set callback functions
set(H.slider1,'Callback',{@sliderCB,H});
set(H.edit1,'Callback',{@editCB,H});
set(H.slider2,'Callback',{@sliderCB,H});
set(H.edit2,'Callback',{@editCB,H});
set(H.slider3,'Callback',{@sliderCB,H});
set(H.edit3,'Callback',{@editCB,H});
set(H.slider4,'Callback',{@sliderCB,H});
set(H.edit4,'Callback',{@editCB,H});
set(H.slider5,'Callback',{@sliderCB,H});
set(H.edit5,'Callback',{@editCB,H});
set(H.slider6,'Callback',{@sliderCB,H});
set(H.edit6,'Callback',{@editCB,H});
set(H.slider7,'Callback',{@sliderCB,H});
set(H.edit7,'Callback',{@editCB,H});
set(H.slider8,'Callback',{@sliderCB,H});
set(H.edit8,'Callback',{@editCB,H});
set(H.slider9,'Callback',{@sliderCB,H});
set(H.edit9,'Callback',{@editCB,H});
set(H.slider10,'Callback',{@sliderCB,H});
set(H.edit10,'Callback',{@editCB,H});
set(H.orbits,'Callback',{@animateOrbit,H});
set(H.anim,'CallBack',{@animateOrbit,H});
set(H.pref,'Callback',{@animateOrbit,H});
set(H.speed,'Callback',{@animateOrbit,H});

% plot initial orbit
animateOrbit(0,0,H);

function sliderCB(~,~,H)
% callback function for each slider, just resets the edit text strings to
% match the values of the sliders
set(H.edit1,'String',get(H.slider1,'Value'));
set(H.edit2,'String',get(H.slider2,'Value'));
set(H.edit3,'String',get(H.slider3,'Value'));
set(H.edit4,'String',get(H.slider4,'Value'));
set(H.edit5,'String',get(H.slider5,'Value'));
set(H.edit6,'String',get(H.slider6,'Value'));
set(H.edit7,'String',get(H.slider7,'Value'));
set(H.edit8,'String',get(H.slider8,'Value'));
set(H.edit9,'String',get(H.slider9,'Value'));
set(H.edit10,'String',get(H.slider10,'Value'));
animateOrbit(0,0,H);

function editCB(~,~,H)
% callback function for each edit text, just resets the slider values to
% match the values of the edit text strings (must check that the value is
% within the correct bounds of the slider so it does not cause an error)
set(H.slider1,'Value',checkBound(str2double(get(H.edit1,'String')),.723,1.523));
set(H.slider2,'Value',checkBound(str2double(get(H.edit2,'String')),0,.9));
set(H.slider3,'Value',checkAngle(str2double(get(H.edit3,'String')),0,360));
set(H.slider4,'Value',checkAngle(str2double(get(H.edit4,'String')),0,360));
set(H.slider5,'Value',checkAngle(str2double(get(H.edit5,'String')),0,360));
set(H.slider6,'Value',checkBound(str2double(get(H.edit6,'String')),.723,1.523));
set(H.slider7,'Value',checkBound(str2double(get(H.edit7,'String')),0,.9));
set(H.slider8,'Value',checkAngle(str2double(get(H.edit8,'String')),0,360));
set(H.slider9,'Value',checkAngle(str2double(get(H.edit9,'String')),0,360));
set(H.slider10,'Value',checkAngle(str2double(get(H.edit10,'String')),0,360));
animateOrbit(0,0,H);

function [a1,a2,e1,e2,i1,i2,R1,R2,w1,w2,sp] = read(H)
% reads in all the necessary data once for use in other functions
a1 = get(H.slider1,'Value');
e1 = get(H.slider2,'Value');
i1 = get(H.slider3,'Value')*pi/180;
R1 = get(H.slider4,'Value')*pi/180;
w1 = get(H.slider5,'Value')*pi/180;            
a2 = get(H.slider6,'Value');
e2 = get(H.slider7,'Value');
i2 = get(H.slider8,'Value')*pi/180;
R2 = get(H.slider9,'Value')*pi/180;
w2 = get(H.slider10,'Value')*pi/180;
sp = get(H.speed,'Value');

function [x,y,z,t] = calculateOrbit(a,e,i,R,w)
% This function will find the time and Cartesian coordinates of the orbits.
mu = .000295913; % au^3 / day^2    standard gravitational parameter of sun
M = linspace(0,2*pi*sqrt((a^3)/mu),100);
E = M;

% calculate eccentric anomaly to find time along the first orbit
if e > 0
    num = 0;
    error = 1;
    E = M./(1-e);
    inds = E > sqrt(6*(1-e)./e);
    E(inds) = (6*M(inds)./e).^(1/3);
    
    % Newton Raphson iteration to compute the eccentric anomaly
    while (error > eps(2*pi)) && (num <1000)
        E = E-(M-E+e.*sin(E))./(e.*cos(E)-1);
        error = max(abs(M-(E-e.*sin(E))));
        num = num+1;
    end
end
% compute time from eccentric anomaly
t = (E/sqrt((a^3)/mu)-e*sin(E/sqrt((a^3)/mu)))/sqrt(mu/((a)^3));

% Before calculating the direction cosine matrix for the orbit(s), it is
% necessary to check if the orbit is inclined or not because the right 
% ascension of the ascending node is undefined for non-inclined orbits.
if i == 0
    R = 0;
end

% define the direction cosine matrix to transform in to the inertial frame
dcm = [cos(R)*cos(w)*cos(i)-sin(R)*sin(w),-cos(R)*cos(i)*sin(w)-sin(R)*cos(w),-cos(R)*sin(i);...
    sin(R)*cos(w)*cos(i)+cos(R)*sin(w),-sin(R)*sin(w)*cos(i)+cos(R)*cos(w),-sin(R)*sin(i);...
    sin(i)*cos(w),-sin(i)*sin(w),cos(i)];

% calculate orbit using the parameters, from body frame into inertial
% note that b3 = 0 since the body frame of an orbit is located in the plane
% of that orbit
% also x and y are negative for aesthetic purposes only, to orient the
% orbit in such a way so apoapse is nearest to the user and so the body
% orbits counter clockwise as viewed from above in a prograde orbit
b1 = a*cos(2*pi*t/t(end))+a*e^2;
b2 = a*sqrt(1-e^2)*sin(2*pi*t/t(end));
x = -(dcm(1,1)*b1+dcm(1,2)*b2);
y = -(dcm(2,1)*b1+dcm(2,2)*b2);
z = dcm(3,1)*b1+dcm(3,2)*b2;

function animateOrbit(~,~,H)
% read in the data and coordinates of the orbit
[a1,a2,e1,e2,i1,i2,R1,R2,w1,w2,sp] = read(H);
[x1,y1,z1,t1] = calculateOrbit(a1,e1,i1,R1,w1);
[x2,y2,z2,t2] = calculateOrbit(a2,e2,i2,R2,w2);
mu = .000295913; %au^3 / day^2

% Calculate the ratio of periods. This will be used to ensure that the
% larger orbit has a longer period in the animation
max_a1 = get(H.slider1,'Max');
max_a2 = get(H.slider6,'Max'); 
if a1 > a2
    n = [1,(a1/a2)^1.5].*(max_a1/a1)^3;
else 
    n = [(a2/a1)^1.5,1].*(max_a2/a2)^3;
end

% save previous view
prevView = H.ax2.View;

% if two orbits are preferred by the user then plot two, else only one
if get(H.orbits,'Value') == 1
    plot3(H.ax2,x1,y1,z1,'--w',x2,y2,z2,'--w','LineWidth',1.5);
else
    plot3(H.ax2,x1,y1,z1,'--w','LineWidth',1.5);
end
hold on;
set(H.ax2,'xcolor','none','ycolor','none','zcolor','none','Color','none')

% Create a set of axes, set the axes equal, and make the axis labels invisible.
% This will enhance visual aid by providing lines of reference.
ax = -get(H.slider1,'Max')*(1.1+get(H.slider2,'Max')):get(H.slider1,'Max')*(1.1+get(H.slider2,'Max'));
plot3(H.ax2,ax,zeros(size(ax)),zeros(size(ax)),'r',zeros(size(ax)),ax,zeros(size(ax)),'g',...
    zeros(size(ax)),zeros(size(ax)),ax,'b','LineWidth',2);
axis equal;
if get(H.pref,'Value') == 1
    rotate3d on;
elseif get(H.pref,'Value') == 2
    zoom on;
end

% time index
ind1 = round(length(t1)/2); 
ind2 = ind1;

% create the data for the orbiting bodies and sun (centered at the origin)
[x,y,z] = sphere;
surf(H.ax2,max_a1*x/12,max_a1*y/12,max_a1*z/12,'EdgeColor','none','FaceColor','y');
planet1 = surf(H.ax2,max_a1*x/15+x1(ind1),max_a1*y/15+y1(ind1),max_a1*z/15+z1(ind1),...
    'EdgeColor','none','FaceColor','g');
if get(H.orbits,'Value') == 1
    planet2 = surf(H.ax2,max_a2*x/15+x2(ind2),max_a2*y/15+y2(ind2),max_a2*z/15+z2(ind2),...
        'EdgeColor','none','FaceColor','r');
end
hold off;

% apply previous view
H.ax2.View = prevView;

% reset the string of the toggle button depending on its current state
if get(H.anim,'Value') == 1
    set(H.anim,'String','Pause');
else
    set(H.anim,'String','Play');
end

% Update the data if the the user preference is set to play animation. 
% Continue animation as long as it is selected
while get(H.anim,'Value') == 1
    % update time but check that it is in the appropriate range
    ind1 = ind1 + sp*n(1);
    if ind1 > length(x1)
        ind1 = ind1 - length(x1);
    end
    ind2 = ind2 + sp*n(2);
    if ind2 > length(x2)
        ind2 = ind2 - length(x2);
    end

    % update position of the planets depending on the updated time
    set(planet1,'XData',max_a1*x/15+x1(ceil(ind1)),'YData',max_a1*y/15+y1(ceil(ind1)),...
        'ZData',max_a1*z/15+z1(ceil(ind1)));
    if get(H.orbits,'Value') == 1 && get(H.anim,'Value') == 1
        set(planet2,'XData',max_a2*x/15+x2(ceil(ind2)),'YData',max_a2*y/15+y2(ceil(ind2)),...
            'ZData',max_a2*z/15+z2(ceil(ind2)));
    end
    drawnow
    pause(eps);

    % print the distances and velocities to the figure window to enhance
    % understanding of the effects of the orbital parameters
    r1 = sqrt(x1(ceil(ind1))^2+y1(ceil(ind1))^2+z1(ceil(ind1))^2);
    v1 = sqrt(mu*(2/r1-1/a1));
    r2 = sqrt(x2(ceil(ind2))^2+y2(ceil(ind2))^2+z2(ceil(ind2))^2);
    v2 = sqrt(mu*(2/r2-1/a2));
    s1 = sprintf('r1 = %2.4f au, v1 = %2.4f au/day, T1 = %2.2f days\n',r1,v1,t1(end));
    s2 = sprintf('r2 = %2.4f au, v2 = %2.4f au/day, T2 = %2.2f days\n',r2,v2,t2(end));
    if get(H.orbits,'Value') == 1
        set(H.text15,'String',[s1 s2]);
    else
        set(H.text15,'String',s1);
    end
end

% print values
printParameters(H);

function printParameters(H)
% read in data and print to the figure after every update in the data
[a1,a2,e1,e2,i1,i2,R1,R2,w1,w2,~] = read(H);
s = sprintf('a1 = %2.3f au, e1 = %2.2f, i1 = %2.2f deg, RAAN1 = %2.2f deg, w1 = %2.2f deg\n',a1,e1,i1,R1,w1);
if get(H.orbits,'Value') == 1    
    s = [s sprintf('a2 = %2.3f au, e2 = %2.2f, i2 = %2.2f deg, RAAN2 = %2.2f deg, w2 = %2.2f deg',a2,e2,i2,R2,w2)];
end
set(H.text14,'String',s,'FontWeight','Bold','FontSize',10,'BackgroundColor','k','ForegroundColor','w');

function var = checkBound(var,low,high)
% checks to make sure that the slider value is within the appropriate range
if var < low 
    warning('variable out of range, automatically correcting to nearest in range value');
    var = low;
elseif var > high
    warning('variable out of range, automatically correcting to nearest in range value');
    var = high;
end

function var = checkAngle(var,low,high)
% checks to make sure that the angle value is within the appropriate range
while var < low
    var = var + (high - low);
end
while var > high
    var = var - (high - low);
end