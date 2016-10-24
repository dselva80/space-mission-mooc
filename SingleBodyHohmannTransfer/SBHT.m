function SBHT

% str2double 

% set up display axes for user understanding of the angles
ops = {'BackgroundColor','ForegroundColor','FontWeight','FontSize'};
H.fig = figure('Name','Single_Body_Hohmann_Transfer','Position',...
    [10 50 1350 600],'Color','k','Resize','off');

% commands
H.animPref_text = uicontrol('Style','text','Position',...
    [600 575 150 20],'String','Animation Style',...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10);
H.planets_text = uicontrol('Style','text','Position',...
    [600 510 150 20],'String','Planet',...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10);
H.text_speed = uicontrol('Style','text','Position',...
    [600 380 150 20],'String','Animation Speed',...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10);
H.animPref = uicontrol('Style','popupmenu','Position',...
    [625 555 100 20],ops{1},'no',ops{2},'w',ops{3},...
    'Bold',ops{4},10,'Value',1,'String',{'Automatic','Manual'});
H.planets = uicontrol('Style','popupmenu','Position',[625 490 100 20],...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10,'Value',3,...
    'String',{'Mercury','Venus','Earth','Mars','Jupiter',...
    'Saturn','Uranus','Neptune','Moon'});
H.print_planet = uicontrol('Style','text','Position',[600 400 150 80],...
    'String',[sprintf('\nmu = 398600.00 (km^3)/(s^2) \n'),...
    sprintf('      R = 6378.10 km')]...
    ,ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10);
H.parameters = uicontrol('Style','text','Position',...
    [700 50 200 200],'String','',ops{1},'no',...
    ops{2},'w',ops{3},'Bold',ops{4},10);
H.speed = uicontrol('Style','slider','Position',...
    [625 350 100 25],'Min',.5,'Max',2,'Value',1.4,...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10);
H.animate = uicontrol('Style','togglebutton','Position',...
    [625 310 100 25],'Value',0,'String','Play',...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10);

% first orbit parameters
H.text_orb1 = uicontrol('Style','text','Position',...
    [750 575 150 20],'String','-- Initial Orbit --',...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10);
H.text_altp1 = uicontrol('Style','text','Position',...
    [750 550 150 20],'String','alt at periapse (km)',...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10);
H.text_alta1 = uicontrol('Style','text','Position',...
    [750 500 150 20],'String','alt at apoapse (km)',...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10);
H.text_i1 = uicontrol('Style','text','Position',...
    [750 450 150 20],'String','inclination = [0,180]',...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10);
H.text_R1 = uicontrol('Style','text','Position',...
    [750 400 150 20],'String','RAAN = [0,360]',...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10);
H.text_w1 = uicontrol('Style','text','Position',...
    [760 335 130 40],'String','argument of periapse = [0,360]',...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10);
H.altp1 = uicontrol('Style','edit','Position',...
    [775 530 100 20],'String','200',...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10);
H.alta1 = uicontrol('Style','edit','Position',...
    [775 480 100 20],'String','5000',...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10);
H.i1 = uicontrol('Style','edit','Position',...
    [775 430 100 20],'String','30',...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10);
H.R1 = uicontrol('Style','edit','Position',...
    [775 380 100 20],'String','10',...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10);
H.w1 = uicontrol('Style','edit','Position',...
    [775 315 100 20],'String','10',...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10);

% second orbit parameters
H.text_orb2 = uicontrol('Style','text','Position',...
    [890 575 150 20],'String','-- Final Orbit --',...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10);
H.text_altp2 = uicontrol('Style','text','Position',...
    [890 550 150 20],'String','alt at periapse (km)',...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10);
H.text_alta2 = uicontrol('Style','text','Position',...
[890 500 150 20],'String','alt at apoapse (km)',...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10);
H.text_i2 = uicontrol('Style','text','Position',...
    [890 450 150 20],'String','inclination = [0,180]',...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10);
H.text_R2 = uicontrol('Style','text','Position',...
    [890 400 150 20],'String','RAAN = [0,360]',...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10);
H.text_w2 = uicontrol('Style','text','Position',...
    [900 335 130 40],'String','argument of periapse = [0,360]',...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10);
H.altp2 = uicontrol('Style','edit','Position',...
    [915 530 100 20],'String','200',...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10);
H.alta2 = uicontrol('Style','edit','Position',...
    [915 480 100 20],'String','200',...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10);
H.i2 = uicontrol('Style','edit','Position',...
    [915 430 100 20],'String','0',...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10);
H.R2 = uicontrol('Style','edit','Position',...
    [915 380 100 20],'String','0',...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10);
H.w2 = uicontrol('Style','edit','Position',...
    [915 315 100 20],'String','0',...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10);

% manual commands (all visible off unless manual animation style selected)
% first burn
H.text_burn1 = uicontrol('Style','text','Position',...
    [1040 575 150 20],'String','-- First Burn --',...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10); 
H.text_nu1 = uicontrol('Style','text','Position',...
    [1050 530 130 40],'String','true anomaly in orbit = [0,360]',...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10);
H.text_dv1 = uicontrol('Style','text','Position',...
    [1040 485 150 20],'String','delta v (km/s)',...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10);
H.text_imp1 = uicontrol('Style','text','Position',...
    [1050 415 130 40],'String','thrust impulse angle = [0,360]',...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10);
H.text_th1 = uicontrol('Style','text','Position',...
    [1040 345 150 40],'String','inclination change angle = [0,360]',...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10);
H.nu1 = uicontrol('Style','edit','Position',...
    [1065 510 100 20],'String','45',...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10);
H.dv1 = uicontrol('Style','text','Position',...
    [1065 465 100 20],'String','.4',...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10);
H.imp1 = uicontrol('Style','text','Position',...
    [1065 395 100 20],'String','0',...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10);
H.th1 = uicontrol('Style','text','Position',...
    [1065 325 100 20],'String','0',...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10);
% second burn
H.text_burn2 = uicontrol('Style','text','Position',...
    [1180 575 150 20],'String','-- Second Burn --',...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10); 
H.text_nu2 = uicontrol('Style','text','Position',...
    [1190 530 130 40],'String','true anomaly in orbit = [0,360]',...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10);
H.text_dv2 = uicontrol('Style','text','Position',...
    [1180 485 150 20],'String','delta v (km/s)',...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10);
H.text_imp2 = uicontrol('Style','text','Position',...
    [1190 415 130 40],'String','thrust impulse angle = [0,360]',...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10);
H.text_th2 = uicontrol('Style','text','Position',...
    [1180 345 150 40],'String','inclination change angle = [0,360]',...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10);
H.nu2 = uicontrol('Style','edit','Position',...
    [1205 510 100 20],'String','0',...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10);
H.dv2 = uicontrol('Style','text','Position',...
    [1205 465 100 20],'String','0',...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10);
H.imp2 = uicontrol('Style','text','Position',...
    [1205 395 100 20],'String','0',...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10);
H.th2 = uicontrol('Style','text','Position',...
    [1205 325 100 20],'String','0',...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10);

% set axes to display the thrust angles to the user
H.disp = axes('Units','Pixels','Position',[900 -50 450 450],...
    'xcolor','no','ycolor','no','zcolor','no','Color','no');
set(H.disp,'View',[-42,35]);
display(H);

% set axes to plot orbits and animate
H.ax = axes('Units','Pixels','Position',[-50 0 700 700],...
    'xcolor','no','ycolor','no','zcolor','no','Color','no');
H.ax.View = [-100,20];

% set callback functions
set(H.planets,'Callback',{@CB,H})
set(H.animPref,'Callback',{@CB,H})
set(H.speed,'Callback',{@CB,H})
set(H.animate,'Callback',{@CB,H})
set(H.altp1,'Callback',{@CB,H})
set(H.alta1,'Callback',{@CB,H})
set(H.i1,'Callback',{@CB,H})
set(H.R1,'Callback',{@CB,H})
set(H.w1,'Callback',{@CB,H})
set(H.nu1,'Callback',{@CB,H})
set(H.altp2,'Callback',{@CB,H})
set(H.alta2,'Callback',{@CB,H})
set(H.i2,'Callback',{@CB,H})
set(H.R2,'Callback',{@CB,H})
set(H.w2,'Callback',{@CB,H})
set(H.nu2,'Callback',{@CB,H})
set(H.dv1,'Callback',{@CB,H})
set(H.th1,'Callback',{@CB,H})
set(H.imp1,'Callback',{@CB,H})
set(H.dv2,'Callback',{@CB,H})
set(H.th2,'Callback',{@CB,H})
set(H.imp2,'Callback',{@CB,H})
end

function CB(~,~,H)
% set the correct styles depending on the selected animation style
if get(H.animPref,'Value') == 1
    style1 = 'edit';
    style2 = 'text';
else
    style1 = 'text';
    style2 = 'edit';
end
set(H.altp2,'Style',style1)
set(H.alta2,'Style',style1)
set(H.i2,'Style',style1)
set(H.R2,'Style',style1)
set(H.w2,'Style',style1)
set(H.dv1,'Style',style2)
set(H.th1,'Style',style2)
set(H.imp1,'Style',style2)
set(H.dv2,'Style',style2)
set(H.th2,'Style',style2)
set(H.imp2,'Style',style2)

if str2double(get(H.altp1,'String')) > str2double(get(H.alta1,'String')) 
    tmp = sprintf('%2.0f',str2double(get(H.altp1,'String')));
    set(H.altp1,'String',get(H.alta1,'String'))
    set(H.alta1,'String',tmp)
end
if str2double(get(H.altp2,'String')) > str2double(get(H.alta2,'String')) 
    tmp = sprintf('%2.0f',str2double(get(H.altp2,'String')));
    set(H.altp2,'String',get(H.alta2,'String'))
    set(H.alta2,'String',tmp)
end

% update the thrust angle display
display(H)

% continue to animation/calculation
if get(H.animPref,'Value') == 1
    calcAuto(H)
elseif get(H.animPref,'Value') == 2
    calcManual(H)
end
end

function [af,ef,If,Rf,wf,nuf,h] = calcs(a0,e0,I0,R0,w0,nu0,dv,imp,th,mu)

% direction cosine matrix to transform between frames
dcm1 = [cos(w0)*cos(I0)*cos(R0)-sin(w0)*sin(R0),...
    cos(w0)*cos(I0)*sin(R0)+sin(w0)*cos(R0),cos(w0)*sin(I0);...
    -sin(w0)*cos(I0)*cos(R0)-cos(w0)*sin(R0),...
    -sin(w0)*cos(I0)*sin(R0)+cos(w0)*cos(R0),-sin(w0)*sin(I0);...
    -sin(I0)*cos(R0),-sin(I0)*sin(R0),cos(I0)];
dcm1 = dcm1';

% all calculations initally set in orbital plane and then rotated into
% inertial frame using dcm1
% position - scalar and vector
r_mag = a0*(1-e0^2)/(1+e0*cos(nu0));
r = r_mag*[cos(nu0);sin(nu0);0];

% tangent of the orbit at the point of the burn specified by nu1
dydx = -sqrt(1-e0^2)*(r(1)+a0*e0)/sqrt(a0^2-(r(1)+a0*e0)^2);
if nu0 > pi
    dydx = -dydx;
end

% rotate radius into perifocal frame
r = dcm1*r;

% inital velocity in the direction of the derivative
v0 = -sqrt(mu*(2/r_mag-1/a0))*[cos(atan(dydx));sin(atan(dydx));0];
if nu0 > pi
    v0 = -v0;
end
v0 = dcm1*v0;

% vectorize the delta v and add to initial velocity vector 
dv = dv*[cos(nu0+pi/2-imp)*cos(th);sin(nu0+pi/2-imp)*cos(th);sin(th)]; 
dv = dcm1*dv;
vf = [v0(1)+dv(1);v0(2)+dv(2);v0(3)+dv(3)];
vf_mag = norm(vf);

% unit vectors
j = [0;1;0];
k = [0;0;1];

% specific angular momentum and magnitude
h = cross(r,vf);
h_mag = norm(h);

% specific orbital energy
energy = .5*vf_mag^2-mu/r_mag;

% nodal vector and magnitude
n_vec = cross(k,h)/h_mag;
n_mag = norm(n_vec);

% calculate the (second) orbital elements
af = -mu/2/energy;
ecc2 = cross(vf,h)/mu-r/r_mag;
ef = norm(ecc2);
nuf = acos(dot(ecc2,r)/r_mag/ef);
if dot(r,vf) < 0
    nuf = 2*pi - nuf;
end
If = acos(dot(h,k)/h_mag);
if If ~= 0 
    wf = acos(dot(n_vec,ecc2)/ef/n_mag);
    if ecc2(3) < 0
        wf = 2*pi - wf;
    end
    wf = wf-pi/2;
    Rf = acos(dot(j,n_vec)/n_mag);
    if n_vec(2) < 0 
        Rf = 2*pi - Rf;
    end
    Rf = Rf-pi;
end
if If == 0
    wf = atan2(ecc2(2),ecc2(1));
end
end

function calcAuto(H)
%[a1,e1,i1,R1,w1,nu1,a2,e2,i2,R2,w2,nu2,sp,~,~,~,~,~,~,mu,R,col] = parameters(H);
end

function calcManual(H)
% read in all data
[a1,e1,i1,R1,w1,nu1,~,~,~,~,~,nu2,sp,dv1,dv2,...
    imp1,imp2,th1,th2,mu,Radius,col] = parameters(H);

% calculate orbital parameters depending on user inputs
[a1,e1,i1,R1,w1,nu1,h1] = calcs(a1,e1,i1,R1,w1,nu1,0,imp1,th1,mu);
[aT,eT,iT,RT,wT,nuT,hT] = calcs(a1,e1,i1,R1,w1,nu1,dv1,imp1,th1,mu);
[a2,e2,i2,R2,w2,nu3,h2] = calcs(aT,eT,iT,RT,wT,nuT,dv2,imp2,th2,mu);

% update the final and initial orbit texts for the user to view
s = [sprintf('a1 = %2.0f km  e1 = %2.4f  \n',a1,e1),...
    sprintf('a2 = %2.0f km  e2 = %2.4f  ',a2,e2)];
set(H.parameters,'String',s);

% reset text values for user display
altp2 = a2*(1-e2)-Radius;
alta2 = a2*(1+e2)-Radius;
if altp2 < 160
    warning('altitude extremely low or negative')
end
set(H.altp2,'String',sprintf('%1.0f',altp2));
set(H.alta2,'String',sprintf('%1.0f',alta2));
set(H.i2,'String',sprintf('%1.0f',i2*180/pi));
set(H.R2,'String',sprintf('%1.0f',R2*180/pi));
set(H.w2,'String',sprintf('%1.0f',w2*180/pi));

% calculate the orbits
[r1,t1,t_inf1] = orbit(a1,e1,i1,R1,w1,h1,mu);
x1 = r1(1,:);
y1 = r1(2,:);
z1 = r1(3,:);
[rT,tT,t_infT] = orbit(aT,eT,iT,RT,wT,hT,mu);
xT = rT(1,:);
yT = rT(2,:);
zT = rT(3,:);
%zT = zT + rad1(3,1);
[r2,t2,t_inf2] = orbit(a2,e2,i2,R2,w2,h2,mu);
x2 = r2(1,:);
y2 = r2(2,:);
z2 = r2(3,:);
%z2 = z2 + rad1(3,1) + rad2(3,1);

% find index for animation
ind1 = calcIndex(t1,t_inf1,nu1,e1);
ind2 = calcIndex(tT,t_infT,nuT,eT);
ind3 = calcIndex(tT,t_infT,nu2,eT);
ind4 = calcIndex(t2,t_inf2,nu3,e2);

% redefine the orbits to meet index dimensions
r1 = [r1,r1(:,1:ind1)];
if eT < 1
    rT = [rT(:,ind2:end),rT(:,1:ind3)];
else
    rT = rT(:,ind2:ind3);
end
if e2 < 1 
    r2 = [r2(:,ind4:end),r2];
else
    r2 = r2(:,ind4:end);
end
path = [r1(1,:),rT(1,:),r2(1,:);...
    r1(2,:),rT(2,:),r2(2,:);...
    r1(3,:),rT(3,:),r2(3,:)];

% if e1 == 0 && e2 == 0
%     path = [x1,x1,x2;y1,y1,y2;z1,z1,z2];
% end

% create axes and axis limits
axlim = real(1.1*max(max(a2*(1+e2)),max(a1*(1+e1),aT*(1+eT))));
ax = linspace(-axlim,axlim,10);

% create value to correctly animate resonance and correct the time
if a1 > a2
    n = [1,(a1/aT)^1.5,(a1/a2)^1.5].*(axlim/a1)^3;
else 
    n = [(a2/a1)^1.5,(a2/aT)^1.5,1].*(axlim/a2)^3;
end

% save previous view
prevView = H.ax.View;

% plot orbits
plot3(H.ax,ax,zeros(size(ax)),zeros(size(ax)),'r',...
    zeros(size(ax)),ax,zeros(size(ax)),'g',...
    zeros(size(ax)),zeros(size(ax)),ax,'b','LineWidth',2);
hold(H.ax,'on');
plot3(H.ax,x1,y1,z1,'w',xT,yT,zT,'y',x2,y2,z2,'m');
%plot3(H.ax,path(1,:),path(2,:),path(3,:),'w');
set(H.ax,'Color','no','xcolor','no','ycolor','no','zcolor','no')

% create the objects for the surfaces
ind = 1;
[x,y,z] = sphere;
surf(H.ax,axlim*x/12,axlim*y/12,axlim*z/12,...
'EdgeColor','none','FaceColor',col);
sat = surf(H.ax,axlim*x/15+x1(ind),axlim*y/15+y1(ind),...
    axlim*z/15+z1(ind),'EdgeColor','none','FaceColor','c');
hold(H.ax,'off');

% apply previous view
H.ax2.View = prevView;

% reset the string of the toggle button depending on its current state
if get(H.animate,'Value') == 1
    set(H.animate,'String','Pause');
else
    set(H.animate,'String','Play');
end
    
% animate if selected 
while get(H.animate,'Value') == 1 && ceil(ind) < length(path)
    if ind > ind1+length(x1)+ind2+length(xT)*2+ind3
        ind = ind + sp*n(3);
    elseif ind > ind1 && ind < ind1+length(x1)+ind2+length(xT)*2+ind3
        ind = ind + sp*n(2);
    elseif ind < ind1+length(x1)
        ind = ind + sp*n(1);
    end
    if ceil(ind) > length(path) && e2 < 1
        ind = ind - length(x2);
    end
    set(sat,'XData',axlim*x/15+path(1,ceil(ind)),...
        'YData',axlim*y/15+path(2,ceil(ind)),...
        'ZData',axlim*z/15+path(3,ceil(ind)));
    drawnow;
    pause(eps);
end
end

function ind = calcIndex(t,t_inf,nu,e)
ind = 1;
if e == 0
    ind = round(nu*length(t)/2/pi);
elseif e < 1
    while ind < length(t) && nu > t_inf*(t(ind)-t(end/2))/t(end)
        ind = ind + 1;
    end
end
ind = checkAngle(ind+length(t)/2,1,length(t));
end

function [a1,e1,i1,R1,w1,nu1,a2,e2,i2,R2,w2,nu2,...
    sp,dv1,dv2,imp1,imp2,th1,th2,mu,R,col] = parameters(H)
% read in the data to retrieve the planet the user selected and define the
% standard gravitational parameter and color accordingly
if get(H.planets,'Value') == 1
    mu = 22032;         % standard gravitational parameter
    col = [1 .5 .5];    % color of planet to display (optional)
    R = 2439.7;         % radius of planet
elseif get(H.planets,'Value') == 2
    mu = 324860; 
    col = [.9 .6 .1];
    R = 6051.8;
elseif get(H.planets,'Value') == 3
    mu = 398600; 
    col = [0 1 0];
    R = 6378.1;
elseif get(H.planets,'Value') == 4
    mu = 42828; 
    col = [1 0 0];
    R = 3396.2;
elseif get(H.planets,'Value') == 5
    mu = 1.26687*10^8; 
    col = [1 .5 .3];
    R = 71492;
elseif get(H.planets,'Value') == 6
    mu = 3.79311*10^7; 
    col = [1 .8 0];
    R = 60268;
elseif get(H.planets,'Value') == 7
    mu = 5.79394*10^6; 
    col = [.7 .9 .9];
    R = 25559;
elseif get(H.planets,'Value') == 8
    mu = 6.83653*10^6; 
    col = [.7 .9 .9];
    R = 24764;
elseif get(H.planets,'Value') == 9
    mu = 4904.87; 
    col = [1 1 1];
    R = 1738.1;
end

% update the planet parameters for user
set(H.print_planet,'HorizontalAlignment','Center','String',...
    sprintf('\nmu = %2.2f (km^3)/(s^2)\n      R = %2.2f km',mu,R))

% check boundaries
set(H.altp1,'String',sprintf('%2.0f',...
    checkBound(str2double(get(H.altp1,'String')),200,inf)));
set(H.alta1,'String',sprintf('%2.0f',...
    checkBound(str2double(get(H.alta1,'String')),200,inf)));
set(H.altp2,'String',sprintf('%2.0f',...
    checkBound(str2double(get(H.altp2,'String')),200,inf)));
set(H.alta2,'String',sprintf('%2.0f',...
    checkBound(str2double(get(H.alta2,'String')),200,inf)));

% calculate all outputs for future use
a1 = (2*R+str2double(get(H.altp1,'String'))+...
    str2double(get(H.alta1,'String')))/2;
e1 = 1-(R+str2double(get(H.altp1,'String')))/a1;
i1 = checkAngle(str2double(get(H.i1,'String')),0,360)*pi/180;
R1 = checkAngle(str2double(get(H.R1,'String')),0,360)*pi/180;
w1 = checkAngle(str2double(get(H.w1,'String')),0,360)*pi/180;
nu1 = checkAngle(str2double(get(H.nu1,'String')),0,360)*pi/180;
a2 = (2*R+str2double(get(H.altp2,'String'))+...
    str2double(get(H.alta2,'String')))/2;
e2 = 1-(R+str2double(get(H.altp2,'String')))/a1;
i2 = checkAngle(str2double(get(H.i2,'String')),0,360)*pi/180;
R2 = checkAngle(str2double(get(H.R2,'String')),0,360)*pi/180;
w2 = checkAngle(str2double(get(H.w2,'String')),0,360)*pi/180;
nu2 = checkAngle(str2double(get(H.nu2,'String')),0,360)*pi/180;
sp = get(H.speed,'Value');
dv1 = str2double(get(H.dv1,'String'));
th1 = checkAngle(str2double(get(H.th1,'String')),0,360)*pi/180;
imp1 = checkAngle(str2double(get(H.imp1,'String')),0,360)*pi/180;
dv2 = str2double(get(H.dv2,'String'));
th2 = checkAngle(str2double(get(H.th2,'String')),0,360)*pi/180;
imp2 = checkAngle(str2double(get(H.imp2,'String')),0,360)*pi/180;

% check that the values are in an appropriate range
set(H.i1,'String',sprintf('%1.0f',i1*180/pi));
set(H.R1,'String',sprintf('%1.0f',R1*180/pi));
set(H.w1,'String',sprintf('%1.0f',w1*180/pi));
set(H.nu1,'String',sprintf('%1.0f',nu1*180/pi));
set(H.th1,'String',sprintf('%1.0f',th1*180/pi));
set(H.imp1,'String',sprintf('%1.0f',imp1*180/pi));
set(H.i2,'String',sprintf('%1.0f',i2*180/pi));
set(H.R2,'String',sprintf('%1.0f',R2*180/pi));
set(H.w2,'String',sprintf('%1.0f',w2*180/pi));
set(H.nu2,'String',sprintf('%1.0f',nu2*180/pi));
set(H.th2,'String',sprintf('%1.0f',th2*180/pi));
set(H.imp2,'String',sprintf('%1.0f',imp2*180/pi));
end

function [r,t,t_inf] = orbit(a,e,I,R,w,h,mu)
% calculate the orbital shape by calculating eccentric anomaly and
% solving for time 

% circular orbits
if e == 0
    % mean motion
    n = sqrt(mu/(a^3));
    
    t_inf = 2*pi;
    E = linspace(0,t_inf,360);
    t = (E-e*sin(E))/n;

    % calculate coordinates of the orbit
    nu = t_inf*(t-t(end/2))/t(end);
    r_mag = a*(1-e^2)./(1+e*cos(nu));
    r = [r_mag.*cos(nu);r_mag.*sin(nu);zeros(1,length(t))];
    
% elliptical orbits
elseif e < 1
    % mean motion
    n = sqrt(mu/(a^3));
    
    t_inf = 2*pi;
    M = linspace(0,t_inf,360);
    num = 0;
    error = 1;
    E = M./(1-e);
    inds = E > sqrt(6*(1-e)./e);
    E(inds) = (6*M(inds)./e).^(1/3);

    % Newton Raphson iteration
    while (error > eps(2*pi)) && (num < 1000)
        E = E-(M-E+e.*sin(E))./(e.*cos(E)-1);
        error = max(abs(M-(E-e.*sin(E))));
        num = num+1;
    end
    t = (E-e*sin(E))/n;
    
    % calculate coordinates of the orbit
    r_mag = a*(1-e*cos(E));
    nu = atan2(sin(E)*sqrt(1-e^2)./(1+e*cos(E)),(cos(E)-e)./(1+e*cos(E)));
    r = [r_mag.*cos(nu);r_mag.*sin(nu);zeros(1,length(t))];
    
% parabolic orbits
elseif e == 1
    % mean motion
    p = dot(h,h)/mu;
    n = 2*sqrt(mu/(p^3));

    t_inf = pi/4;
    M = linspace(-t_inf,t_inf,360);
    D = M;
    num = 0;
    error = 1;
    
    % Newton Raphson iteration
    while (error > eps(2*pi)) && (num < 1000)
        D = D-(M-D-(D.^3)/3)./(-1-D^2);
        error = max(abs(M-D-(D.^3)/3));
        num = num+1;
    end
    t = (D+((D)^3)/3)/n;
    
    % calculate coordinates of the orbit
    r_mag = .5*p*(1+D.^2);
    nu = atan2(p*D./r_mag,(p-r_mag)./r_mag);
    r = [a*cos(nu)-a*e;a*sqrt(1-e^2)*sin(nu);zeros(1,length(t))];
    
% hyperbolic orbits
elseif e > 1
    % mean motion
    n = sqrt(-mu/(a^3));

    t_inf = atan(sqrt(e^2-1));
    M = linspace(-t_inf,t_inf,360);
    F = M;
    num = 0;
    error = 1;
    
    % Newton Raphson iteration
    while (error > eps(2*pi)) && (num < 1000)
        F = F-(M+F-e.*sinh(F))./(1-e.*cosh(F));
        error = max(abs(M-(e.*sinh(F)-F)));
        num = num+1;
    end
    t = (e*sinh(F)-F)/n;
    
    % calculate coordinates of the orbit
    nu = atan2((-sinh(F)*sqrt(e^2-1)./(1-e*cosh(F))),(cosh(F)-e)./(1-e*cosh(F)));
    r_mag = a*(1-e*cosh(F));
    r = [r_mag.*cos(nu);r_mag.*sin(nu);zeros(1,length(t))];
end

% check if orbit is inclined, if not then set RAAN=0 (undefined for i=0)
if I == 0 
    R = 0;
end

% transform the ellipse depending on the angles
% direction cosine matrix to transform between frames
r = DCM323(r,I,R,w);
end

function var = checkBound(var,low,high)
% checks to make sure that the slider value is within the appropriate range
if var < low 
    var = low;
elseif var > high
    var = high;
end
end

function var = checkAngle(var,low,high)
% checks to make sure that the angle value is within the appropriate range
for i = 1:length(var)
    while var(i) < low
        var(i) = var(i) + (high - low);
    end
    while var(i) > high
        var(i) = var(i) - (high - low);
    end
end
end

function display(H)
% read parameters to display the thrust elements to the user
[a,e,i,R,w,nu,~,~,~,~,~,~,~,dv,~,imp,~,inc,~,mu,PR,col] = parameters(H); 

% coordinates for orbit and sphere, and set axis limit
[r,~] = orbit(a,e,i,R,w,0,mu);
x = r(1,:);
y = r(2,:);
z = r(3,:);
[xs,ys,zs] = sphere;
axis_lim = linspace(-1.2*a*(1+e),1.2*a*(1+e),3);

% save previous view
prevView = H.disp.View;

% predefine string for the legend and redefine if any of the angles = 0
str = {'orbital radius','delta v','true anomaly','inclination angle',...
    'impulse angle','reference direction'};

% radius as scalar and vector
r = a*(1-e^2)/(1+e*cos(nu))*[cos(nu);sin(nu);0];
r = DCM323(r,i,R,w);
quiver3(H.disp,0,0,0,r(1),r(2),r(3),'-r',...
    'LineWidth',2,'AutoScale','off','MaxHeadSize',.5)
hold(H.disp,'on')

% delta v vector (not to scale, only to demonstrate the direction)
if dv ~= 0
    dv = .5*a*[-sin(nu-imp)*cos(inc);cos(nu-imp)*cos(inc);sin(inc)];
    dv = DCM323(dv,i,R,w);
    dv = [dv(1,1),dv(2,1),dv(3,1)];
    quiver3(H.disp,r(1),r(2),r(3),dv(1),dv(2),dv(3),'-y',...
        'LineWidth',2,'AutoScaleFactor',1.5,'MaxHeadSize',.5)
else
    str = {str{1},str{3:end}};
end

% reference line for first axis of perifocal frame (eccentricity vector)
e_line = [axis_lim(1),axis_lim(3);0,0;0,0];
e_line = DCM323(e_line,i,R,w);

% reference line for second axis of perifocal frame
q_line = [0,0;axis_lim(1),axis_lim(3);0,0];
q_line = DCM323(q_line,i,R,w);

% reference line for third axis of perifocal frame
h_line = [0,0;0,0;axis_lim(1),axis_lim(3)];
h_line = DCM323(h_line,i,R,w);

% reference line for thrust impulse angle
perp = [0,-a*cos(nu-pi/2)/1.5;0,-a*sin(nu-pi/2)/1.5;0,0];
perp = DCM323(perp,i,R,w);
perp = [perp(1,:)+r(1);perp(2,:)+r(2);perp(3,:)+r(3)];

% projection of the delta v vector in the orbital plane
if imp == 0 && inc == 0
    dv_proj = [0,0;0,0;0,0];
else
    dv_proj = DCM323([0,.5*a;0,0;0,0],0,pi/2-imp,0);
    dv_proj = [dv_proj(1,:)+a*(1-e^2)/(1+e*cos(nu));...
        dv_proj(2,:);zeros(size(dv_proj(3,:)))];
    dv_proj = DCM323(dv_proj,i,R,w+nu);
end

% add arrows -- draw arcs and add arrow at the end
th = linspace(0,2*pi,180);
% the first index to check for redefining the legend
n = length(str)-3;

% true anomaly arrow
if nu <= pi/45
    str = {str{1:n-1},str{n+1:n+3}};
    nu_curve = [0;0;0];
else
    nu_curve = .5*a*(1-e)*[cos(th(1:floor(nu*180/pi/2)));...
        sin(th(1:floor(nu*180/pi/2)));...
        zeros(1,floor(nu*180/pi/2))];
    nu_curve = DCM323(nu_curve,i,R,w);
    nu_arr = [nu_curve(1,end)-nu_curve(1,end-1),...
        nu_curve(2,end)-nu_curve(2,end-1),...
        nu_curve(3,end)-nu_curve(3,end-1)];
    nu_ind = length(nu_curve)-1;
    quiver3(H.disp,nu_curve(1,nu_ind),nu_curve(2,nu_ind),...
        nu_curve(3,nu_ind),nu_arr(1),nu_arr(2),nu_arr(3),'-b',...
        'LineWidth',2,'AutoScale','off','MaxHeadSize',14);
    n = n + 1;
end

% inclination change angle arrow
if inc == 0
    str = {str{1:n-1},str{n+1:n+2}};
    inc_curve = [0;0;0];
else
    inc_curve = .35*a*[cos(th(1:floor(inc*180/pi/2)));...
        zeros(1,floor(inc*180/pi/2));...
        sin(th(1:floor(inc*180/pi/2)))];
    inc_curve = DCM323(inc_curve,0,pi/2-imp,0);
    inc_curve = [inc_curve(1,:)+a*(1-e^2)/(1+e*cos(nu));...
        inc_curve(2,:);inc_curve(3,:)];
    inc_curve = DCM323(inc_curve,i,R,w+nu);
    inc_arr = [inc_curve(1,end)-inc_curve(1,end-1),...
        inc_curve(2,end)-inc_curve(2,end-1),...
        inc_curve(3,end)-inc_curve(3,end-1)];
    inc_ind = length(inc_curve)-1;
    quiver3(H.disp,inc_curve(1,inc_ind),inc_curve(2,inc_ind),...
        inc_curve(3,inc_ind),inc_arr(1),inc_arr(2),inc_arr(3),'-g',...
        'LineWidth',2,'AutoScale','off','MaxHeadSize',.5);
    n = n + 1;
end

% impulse angle arrow
if imp == 0
    str = {str{1:n-1},str{n+1}};
    imp_curve = [0;0;0];
else
    imp_curve = .35*a*[-sin(th(1:floor(imp*180/pi/2)));...
        cos(th(1:floor(imp*180/pi/2)));...
        zeros(1,floor(imp*180/pi/2))];
    imp_curve = DCM323(imp_curve,0,-imp+(3/180)*pi,0);
    imp_curve = [imp_curve(1,:)+a*(1-e^2)/(1+e*cos(nu));...
        imp_curve(2,:);imp_curve(3,:)];
    imp_curve = DCM323(imp_curve,i,R,w+nu);
    imp_arr = [imp_curve(1,1)-imp_curve(1,2);...
        imp_curve(2,1)-imp_curve(2,2);...
        imp_curve(3,1)-imp_curve(3,2)];
    quiver3(H.disp,imp_curve(1,1),imp_curve(2,1),imp_curve(3,1),...
        imp_arr(1),imp_arr(2),imp_arr(3),'-m',...
        'LineWidth',2,'AutoScale','off','MaxHeadSize',.5);
end

% plot the curves
plot3(H.disp,nu_curve(1,:),nu_curve(2,:),nu_curve(3,:),'b',...
    inc_curve(1,:),inc_curve(2,:),inc_curve(3,:),'g',...
    imp_curve(1,:),imp_curve(2,:),imp_curve(3,:),'m','LineWidth',2);

% plot perifocal frame
plot3(H.disp,e_line(1,:),e_line(2,:),e_line(3,:),'--w',...
    q_line(1,:),q_line(2,:),q_line(3,:),'--w',...
    h_line(1,:),h_line(2,:),h_line(3,:),'--w')

% plot the axes, horizontal and vertical reference lines for the thrust,
% the in plane projection of the delta v vector, and the orbit itself
plot3(H.disp,axis_lim,zeros(size(axis_lim)),zeros(size(axis_lim)),':w',...
    zeros(size(axis_lim)),axis_lim,zeros(size(axis_lim)),':w',...
    zeros(size(axis_lim)),zeros(size(axis_lim)),axis_lim,':w',...
    perp(1,:),perp(2,:),perp(3,:),'--w',...
    dv_proj(1,:),dv_proj(2,:),dv_proj(3,:),'--w',...
    x,y,z,'w');

% plot the central and orbiting bodies
surf(H.disp,3*(PR^(2/3))*xs,3*(PR^(2/3))*ys,3*(PR^(2/3))*zs,...
    'EdgeColor','no','FaceColor',col);
surf(H.disp,2*(PR^(2/3))*xs+r(1),2*(PR^(2/3))*ys+r(2),...
    2*(PR^(2/3))*zs+r(3),'EdgeColor','no','FaceColor','c');

% set axis background, limits, and legend, and allow rotation
set(H.disp,'Color','no','xcolor','no','ycolor','no','zcolor','no');
xlim([-1.2*a*(1+e) 1.2*a*(1+e)]);
axis(H.disp,'equal')
legend(H.disp,str,'Location','EastOutside',...
    'TextColor','w','Color','k','Box','on');
hold(H.disp,'off')

% hold the new view
H.disp.View = prevView;
end

function OUT = DCM323(IN,I,R,w)
% this will be the direction cosine matrix for a typical orbital euler
% angle set (3-2-3), with the variable IN being a 3D column vector or a 
% 3xN matrix where N may be any positive integer
if I == 0
    R = 0;
end
dcm = [cos(w)*cos(I)*cos(R)-sin(w)*sin(R),...
    cos(w)*cos(I)*sin(R)+sin(w)*cos(R),cos(w)*sin(I);...
    -sin(w)*cos(I)*cos(R)-cos(w)*sin(R),...
    -sin(w)*cos(I)*sin(R)+cos(w)*cos(R),-sin(w)*sin(I);...
    -sin(I)*cos(R),-sin(I)*sin(R),cos(I)];
dcm = dcm';
OUT(1,:) = IN(1,:)*dcm(1,1)+IN(2,:)*dcm(1,2)+IN(3,:)*dcm(1,3);
OUT(2,:) = IN(1,:)*dcm(2,1)+IN(2,:)*dcm(2,2)+IN(3,:)*dcm(2,3);
OUT(3,:) = IN(1,:)*dcm(3,1)+IN(2,:)*dcm(3,2)+IN(3,:)*dcm(3,3);
end