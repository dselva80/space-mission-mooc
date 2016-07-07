function Keplerian_Orbits
% set options, open a figure and set the axes
ops = {'BackgroundColor','ForegroundColor','FontWeight','FontSize'};
H.fig = figure('Name','KeplerianOrbits','Position',[120 60 900 600],'Color','k');
H.ax = axes('Units','Pixels','Position',[0 50 600 600],...
    'xcolor','k','ycolor','k','zcolor','k','Color','k');
H.text14 = uicontrol('Style','text','Position',[20 10 560 30],'String','',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.text15 = uicontrol('Style','text','Position',[20 50 560 30],'String','',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);

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
    'Bold',ops{4},10,'Min',0,'Max',.9,'Value',0);
H.edit2 = uicontrol('Style','edit','Position',[625 395 100 25],'String','0',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.text5 = uicontrol('Style','text','Position',[600 360 150 20],'String','i1 = [0,360]',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.slider3 = uicontrol('Style','slider','Position',[615 330 120 30],ops{1},'k',ops{2},'w',ops{3},...
    'Bold',ops{4},10,'Min',0,'Max',360,'Value',0);
H.edit3 = uicontrol('Style','edit','Position',[625 300 100 25],'String','0',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.text6 = uicontrol('Style','text','Position',[600 265 150 20],'String','RAAN1 = [0,360]',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.slider4 = uicontrol('Style','slider','Position',[615 235 120 30],ops{1},'k',ops{2},'w',ops{3},...
    'Bold',ops{4},10,'Min',0,'Max',360,'Value',0);
H.edit4 = uicontrol('Style','edit','Position',[625 205 100 25],'String','0',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.text7 = uicontrol('Style','text','Position',[600 175 150 20],'String','w1 = [0,360]',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.slider5 = uicontrol('Style','slider','Position',[615 145 120 30],ops{1},'k',ops{2},'w',ops{3},...
    'Bold',ops{4},10,'Min',0,'Max',360,'Value',0);
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
    'Bold',ops{4},10,'Min',0,'Max',.9,'Value',0);
H.edit7 = uicontrol('Style','edit','Position',[775 395 100 25],'String','0',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.text11 = uicontrol('Style','text','Position',[750 360 150 20],'String','i2 = [0,360]',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.slider8 = uicontrol('Style','slider','Position',[765 330 120 30],ops{1},'k',ops{2},'w',ops{3},...
    'Bold',ops{4},10,'Min',0,'Max',360,'Value',0);
H.edit8 = uicontrol('Style','edit','Position',[775 300 100 25],'String','0',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.text12 = uicontrol('Style','text','Position',[750 265 150 20],'String','RAAN2 = [0,360]',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.slider9 = uicontrol('Style','slider','Position',[765 235 120 30],ops{1},'k',ops{2},'w',ops{3},...
    'Bold',ops{4},10,'Min',0,'Max',360,'Value',0);
H.edit9 = uicontrol('Style','edit','Position',[775 205 100 25],'String','0',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.text13 = uicontrol('Style','text','Position',[750 175 150 20],'String','w2 = [0,360]',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.slider10 = uicontrol('Style','slider','Position',[765 145 120 30],ops{1},'k',ops{2},'w',ops{3},...
    'Bold',ops{4},10,'Min',0,'Max',360,'Value',0);
H.edit10 = uicontrol('Style','edit','Position',[775 115 100 25],'String','0',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);

% set animation preferences
H.orbits = uicontrol('Style','checkbox','Position',[625 60 100 25],'String','2 Orbits?',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.anim = uicontrol('Style','togglebutton','Position',[625 20 100 25],'String','Play',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.pref = uicontrol('Style','popupmenu','Position',[775 60 100 25],'String',{'Rotate on','Zoom on'},...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.speed = uicontrol('Style','slider','Position',[765 20 120 25],ops{1},'k',ops{2},'w',ops{3},...
    'Bold',ops{4},10,'Min',.5,'Max',4,'Value',2);

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

% set background as space
xpos = 2.9*rand(1,400)-1.2;
ypos = 2.9*rand(1,400)-1.2;
siz = .011*rand(1,400);
th = linspace(0,2*pi,50);
hold on;
for i=1:length(siz)
    % check if background stars are too small to be visible, resize, plot
    if siz(i) <= .001
        siz(i) = 10*siz(i);
    end
    fill(siz(i)*cos(th)+xpos(i),siz(i)*sin(th)+ypos(i),[.1 .1 .1]);
end

% set background color to black and fill the stars (circles) with white
set(gca,'Color','w');
set(gca,'xcolor','w');
set(gca,'ycolor','w');
set(gca,'zcolor','w');
set(gcf,'Color','w');

% save figure as image, produce negative and set as background for figure
saveas(H.ax,'space','jpg');
neg('space.jpg');
ax = axes('units','normalized','position',[0 0 1 1]);
uistack(ax,'bottom');
im = image(imread('space.jpg'));
set(ax,'handlevisibility','off');
set(im,'alphadata',1);
hold off;

function sliderCB(~,~,H)
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
% Initialize angle and anomaly, set standard gravitational parameter of 
% the sun, and read in all of the necessary data.
mu = .000295913; % au^3 / day^2          
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
t = (E/sqrt((a^3)/mu)-e*sin(E/sqrt((a^3)/mu)))/sqrt(mu/((a)^3));

% Before calculating the direction cosine matrix for the orbit(s), it is
% necessary to check if the orbit is inclined or not because the RAAN is
% undefined for non-inclined orbits.
if i == 0
    R = 0;
end

% define the direction cosine matrix to transform in to the inertial frame
dcm = [cos(R)*cos(w)*cos(i)-sin(R)*sin(w),-cos(R)*cos(i)*sin(w)-sin(R)*cos(w),-cos(R)*sin(i);...
    sin(R)*cos(w)*cos(i)+cos(R)*sin(w),-sin(R)*sin(w)*cos(i)+cos(R)*cos(w),-sin(R)*sin(i);...
    sin(i)*cos(w),-sin(i)*sin(w),cos(i)];

% calculate first orbit using the parameters, from body frame into inertial
b1 = a*cos(2*pi*t/t(end))+a*e^2;
b2 = a*sqrt(1-e^2)*sin(2*pi*t/t(end));
x = dcm(1,1)*b1+dcm(1,2)*b2;
y = dcm(2,1)*b1+dcm(2,2)*b2;
z = dcm(3,1)*b1+dcm(3,2)*b2;

function animateOrbit(~,~,H)
% read in the data and coordinates of the orbit
[a1,a2,e1,e2,i1,i2,R1,R2,w1,w2,sp] = read(H);
[x1,y1,z1,t1] = calculateOrbit(a1,e1,i1,R1,w1);
[x2,y2,z2,t2] = calculateOrbit(a2,e2,i2,R2,w2);
mu = .000295913; %au^3 / day^2

% Calculate the ratio of periods. This will be used to ensure that the
% larger orbit has a longer period in the animation
if a1 > a2
    n = [1,a1/a2];
else 
    n = [a2/a1,1];
end

% if two orbits are preferred by the user then plot two, else only one
if get(H.orbits,'Value') == 1
    plot3(x1,y1,z1,'--w',x2,y2,z2,'--w','LineWidth',1.5);
else
    plot3(x1,y1,z1,'--w','LineWidth',1.5);
end
hold on;
% set plot color to transparent so the space background is visible
set(gca,'xcolor','None');
set(gca,'ycolor','None');
set(gca,'zcolor','None');
set(gca,'Color','None');
set(gcf,'Color','None');

% Create a set of axes, set the axes equal, and make the axis labels invisible.
% This will enhance visual aid by providing lines of reference.
ax = -get(H.slider1,'Max')*(1.1+get(H.slider2,'Max')):get(H.slider1,'Max')*(1.1+get(H.slider2,'Max'));
plot3(ax,zeros(1,length(ax)),zeros(1,length(ax)),'r','LineWidth',2);
plot3(zeros(1,length(ax)),ax,zeros(1,length(ax)),'g','LineWidth',2);
plot3(zeros(1,length(ax)),zeros(1,length(ax)),ax,'b','LineWidth',2);
axis equal;
if get(H.pref,'Value') == 1
    rotate3d on;
elseif get(H.pref,'Value') == 2
    zoom on;
end

% read in values and orbit calculations
max_a1 = get(H.slider1,'Max');
max_a2 = get(H.slider6,'Max'); 

% time index
ind1 = round(length(t1)/2); 
ind2 = round(length(t2)/2);

% create the data for the orbiting bodies and sun (centered at the origin)
[x,y,z] = sphere;
sun = fill3(max_a1*x/12,max_a1*y/12,max_a1*z/12,'y');
set(sun,'EdgeColor','y','FaceColor','y','LineWidth',3);
planet1 = fill3(max_a1*x/15+x1(ind1),max_a1*y/15+y1(ind1),max_a1*z/15+z1(ind1),'g');
set(planet1,'EdgeColor','g','FaceColor','g','LineWidth',2);
if get(H.orbits,'Value') == 1
    planet2 = fill3(max_a2*x/15+x2(ind2),max_a2*y/15+y2(ind2),max_a2*z/15+z2(ind2),'r');
    set(planet2,'EdgeColor','r','FaceColor','r','LineWidth',2);
end
hold off;

% Update the data if the the user preference is set to play animation. 
% Continue animation as long as it is selected but no longer than 10
% periods for either orbit. 
periods1 = 0; periods2 = 0;
while (get(H.anim,'Value') == 1) && (periods1 < 10) && (periods2 < 10)
    % update time but check that it is in the appropriate range
    ind1 = ind1 + sp*2*n(1);
    if ind1 > length(x1)
        ind1 = ind1 - length(x1);
        periods1 = periods1 + 1;
    end
    ind2 = ind2 + sp*2*n(2);
    if ind2 > length(x2)
        ind2 = ind2 - length(x2);
        periods2 = periods2 + 1;
    end

    % update position of the planets depending on the updated time
    set(planet1,'XData',max_a1*x/15+x1(ceil(ind1)),'YData',max_a1*y/15+y1(ceil(ind1)),'ZData',max_a1*z/15+z1(ceil(ind1)));
    if get(H.orbits,'Value') == 1 && get(H.anim,'Value') == 1
        set(planet2,'XData',max_a2*x/15+x2(ceil(ind2)),'YData',max_a2*y/15+y2(ceil(ind2)),'ZData',max_a2*z/15+z2(ceil(ind2)));
    end
    drawnow
    pause(eps);

    % print the distances and velocities to the figure window to enhance
    % understanding of the effects of the orbital parameters
    r1 = sqrt(x1(ceil(ind1))^2+y1(ceil(ind1))^2+z1(ceil(ind1))^2);
    v1 = sqrt(mu*(2/r1-1/a1));
    T1 = t1(end);
    r2 = sqrt(x2(ceil(ind2))^2+y2(ceil(ind2))^2+z2(ceil(ind2))^2);
    v2 = sqrt(mu*(2/r2-1/a2));
    T2 = t2(end);
    s1 = sprintf('r1 = %2.4f au, v1 = %2.4f au/day, T1 = %2.2f days\n',r1,v1,T1);
    s2 = sprintf('r2 = %2.4f au, v2 = %2.4f au/day, T2 = %2.2f days\n',r2,v2,T2);
    if get(H.orbits,'Value') == 1
        set(H.text15,'String',[s1 s2]);
    else
        set(H.text15,'String',s1);
    end
end
if periods1 >= 10 || periods2 >= 10
    set(H.anim,'Value',0,'String','Play');
end

% print values
printParameters(H);

function printParameters(H)
% read in data and print to the figure after every update in the data
[a1,a2,e1,e2,i1,i2,R1,R2,w1,w2,~] = read(H);
s1=sprintf('a1 = %2.3f au, e1 = %2.2f, i1 = %2.2f deg, RAAN1 = %2.2f deg, w1 = %2.2f deg\n',a1,e1,i1,R1,w1);
if get(H.orbits,'Value') == 1    
    s2= sprintf('a2 = %2.3f au, e2 = %2.2f, i2 = %2.2f deg, RAAN2 = %2.2f deg, w2 = %2.2f deg',a2,e2,i2,R2,w2);
    set(H.text14,'String',[s1 s2],'FontWeight','Bold','FontSize',10,'BackgroundColor','k','ForegroundColor','w');
else
    set(H.text14,'String',s1,'FontWeight','Bold','FontSize',10,'BackgroundColor','k','ForegroundColor','w');
end

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
    var = var + high - low;
end
while var > high
    var = var - high + low;
end

function neg(filename)
% function takes an image, reads the data and produces a negative version
I = imread(filename);
[a,b,c] = size(I);
for i = 1:a
    for j = 1:b
        for k = 1:c
            I(i,j,k) = 255-I(i,j,k);
        end
    end
end
imwrite(I,filename);