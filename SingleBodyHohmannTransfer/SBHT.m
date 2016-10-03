function SBHT
% set up display axes for user understanding of the angles
ops = {'BackgroundColor','ForegroundColor','FontWeight','FontSize'};
H.fig = figure('Name','Single_Body_Hohmann_Transfer','Position',...
    [10 50 1350 600],'Color','k','Resize','off');

% commands
H.animPref_text = uicontrol('Style','text','Position',[600 575 150 20],'String','Animation Style',...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10);
H.planets_text = uicontrol('Style','text','Position',[600 510 150 20],'String','Planet',...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10);
H.text_speed = uicontrol('Style','text','Position',[600 380 150 20],'String','Animation Speed',...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10);
H.animPref = uicontrol('Style','popupmenu','Position',[625 555 100 20],ops{1},'no',ops{2},'w',ops{3},...
    'Bold',ops{4},10,'Value',1,'String',{'Automatic','Manual'});
H.planets = uicontrol('Style','popupmenu','Position',[625 490 100 20],ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10,...
    'Value',3,'String',{'Mercury','Venus','Earth','Mars','Jupiter','Saturn','Uranus','Neptune','Moon'});
H.print_planet = uicontrol('Style','text','Position',[600 400 150 80],'String',...
    sprintf('\nmu = 398600.00 (km^3)/(s^2) \n      R = 6378.10 km'),ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10);
H.parameters = uicontrol('Style','text','Position',[700 50 200 200],'String','',...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10);
H.speed = uicontrol('Style','slider','Position',[625 350 100 25],'Min',.25,'Max',4,'Value',2,...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10);
H.animate = uicontrol('Style','togglebutton','Position',[625 310 100 25],'Value',0,...
    'String','Play',ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10);

% first orbit parameters
H.text_orb1 = uicontrol('Style','text','Position',[750 575 150 20],'String','-- Initial Orbit --',...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10);
H.text_altp1 = uicontrol('Style','text','Position',[750 550 150 20],'String','alt at periapse (km)',...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10);
H.text_alta1 = uicontrol('Style','text','Position',[750 500 150 20],'String','alt at apoapse (km)',...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10);
H.text_i1 = uicontrol('Style','text','Position',[750 450 150 20],'String','inclination = [0,180]',...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10);
H.text_R1 = uicontrol('Style','text','Position',[750 400 150 20],'String','RAAN = [0,360]',...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10);
H.text_w1 = uicontrol('Style','text','Position',[760 335 130 40],'String','argument of periapse = [0,360]',...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10);
H.altp1 = uicontrol('Style','edit','Position',[775 530 100 20],'String','200',...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10);
H.alta1 = uicontrol('Style','edit','Position',[775 480 100 20],'String','200',...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10);
H.i1 = uicontrol('Style','edit','Position',[775 430 100 20],'String','0',...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10);
H.R1 = uicontrol('Style','edit','Position',[775 380 100 20],'String','0',...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10);
H.w1 = uicontrol('Style','edit','Position',[775 315 100 20],'String','0',...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10);

% second orbit parameters
H.text_orb2 = uicontrol('Style','text','Position',[890 575 150 20],'String','-- Final Orbit --',...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10);
H.text_altp2 = uicontrol('Style','text','Position',[890 550 150 20],'String','alt at periapse (km)',...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10);
H.text_alta2 = uicontrol('Style','text','Position',[890 500 150 20],'String','alt at apoapse (km)',...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10);
H.text_i2 = uicontrol('Style','text','Position',[890 450 150 20],'String','inclination = [0,180]',...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10);
H.text_R2 = uicontrol('Style','text','Position',[890 400 150 20],'String','RAAN = [0,360]',...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10);
H.text_w2 = uicontrol('Style','text','Position',[900 335 130 40],'String','argument of periapse = [0,360]',...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10);
H.altp2 = uicontrol('Style','edit','Position',[915 530 100 20],'String','200',...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10);
H.alta2 = uicontrol('Style','edit','Position',[915 480 100 20],'String','200',...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10);
H.i2 = uicontrol('Style','edit','Position',[915 430 100 20],'String','0',...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10);
H.R2 = uicontrol('Style','edit','Position',[915 380 100 20],'String','0',...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10);
H.w2 = uicontrol('Style','edit','Position',[915 315 100 20],'String','0',...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10);

% manual commands (all visible off unless manual animation style selected)
% first burn
H.text_burn1 = uicontrol('Style','text','Position',[1040 575 150 20],'String','-- First Burn --',...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10); 
H.text_nu1 = uicontrol('Style','text','Position',[1050 530 130 40],'String','true anomaly in orbit = [0,360]',...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10);
H.text_dv1 = uicontrol('Style','text','Position',[1040 485 150 20],'String','delta v (km/s)',...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10);
H.text_imp1 = uicontrol('Style','text','Position',[1050 415 130 40],'String','thrust impulse angle = [0,360]',...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10);
H.text_th1 = uicontrol('Style','text','Position',[1040 345 150 40],'String','inclination change angle = [0,360]',...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10);
H.nu1 = uicontrol('Style','edit','Position',[1065 510 100 20],'String','0',...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10);
H.dv1 = uicontrol('Style','text','Position',[1065 465 100 20],'String','0',...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10);
H.imp1 = uicontrol('Style','text','Position',[1065 395 100 20],'String','0',...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10);
H.th1 = uicontrol('Style','text','Position',[1065 325 100 20],'String','0',...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10);
% second burn
H.text_burn2 = uicontrol('Style','text','Position',[1180 575 150 20],'String','-- Second Burn --',...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10); 
H.text_nu2 = uicontrol('Style','text','Position',[1190 530 130 40],'String','true anomaly in orbit = [0,360]',...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10);
H.text_dv2 = uicontrol('Style','text','Position',[1180 485 150 20],'String','delta v (km/s)',...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10);
H.text_imp2 = uicontrol('Style','text','Position',[1190 415 130 40],'String','thrust impulse angle = [0,360]',...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10);
H.text_th2 = uicontrol('Style','text','Position',[1180 345 150 40],'String','inclination change angle = [0,360]',...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10);
H.nu2 = uicontrol('Style','edit','Position',[1205 510 100 20],'String','0',...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10);
H.dv2 = uicontrol('Style','text','Position',[1205 465 100 20],'String','0',...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10);
H.imp2 = uicontrol('Style','text','Position',[1205 395 100 20],'String','0',...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10);
H.th2 = uicontrol('Style','text','Position',[1205 325 100 20],'String','0',...
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

function [a2,e2,I2,R2,w2,nu2] = calcs(a1,e1,I1,R1,w1,nu1,dv,imp,th,mu)
a2 = a1; e2 = e1; I2 = I1; R2 = R1; w2 = w1; nu2 = nu1;
if dv ~= 0
    % direction cosine matrix to transform between frames
    dcm1 = [cos(w1)*cos(I1)*cos(R1)-sin(w1)*sin(R1),...
        cos(w1)*cos(I1)*sin(R1)+sin(w1)*cos(R1),cos(w1)*sin(I1);...
        -sin(w1)*cos(I1)*cos(R1)-cos(w1)*sin(R1),...
        -sin(w1)*cos(I1)*sin(R1)+cos(w1)*cos(R1),-sin(w1)*sin(I1);...
        -sin(I1)*cos(R1),-sin(I1)*sin(R1),cos(I1)];
    dcm1 = dcm1';
    
    % all calculations initally set in orbital plane and then rotated into
    % inertial frame using dcm1
    % position - scalar and vector
    r_mag = a1*(1-e1^2)/(1+e1*cos(nu1));
    r = r_mag*[cos(nu1);sin(nu1);0];
    % tangent of the orbit at the point of the burn specified by nu1
    dydx = -sqrt(1-e1^2)*(r(1)+a1*e1)/sqrt(a1^2-(r(1)+a1*e1)^2);
    if nu1 > pi
        dydx = -dydx;
    end
    % rotate radius into perifocal frame
    r = dcm1*r;
    % inital velocity in the direction of the derivative
    v0 = -sqrt(mu*(2/r_mag-1/a1))*[cos(atan(dydx));sin(atan(dydx));0];
    if nu1 > pi
        v0 = -v0;
    end
    v0 = dcm1*v0;
    % vectorize the delta v and add to initial velocity vector 
    dv = dv*[cos(nu1+pi/2-imp)*cos(th);sin(nu1+pi/2-imp)*cos(th);sin(th)]; 
    dv = dcm1*dv;
    vf = [v0(1)+dv(1);v0(2)+dv(2);v0(3)+dv(3)];
    vf_mag = sqrt(dot(vf,vf));
	% unit vector
    j = [0;1;0];
    k = [0;0;1];
    % specific angular momentum and magnitude
    h = cross(r,vf);
    h_mag = norm(h);
    % specific orbital energy
    energy = .5*vf_mag^2-mu/r_mag;
    % nodal vector and magnitude
    n_vec = cross(k,h);
    n_mag = norm(n_vec);
    % calculate the (second) orbital elements
    a2 = -mu/2/energy;
    ecc2 = cross(vf,h)/mu-r/r_mag;
    e2 = norm(ecc2);
    nu2 = acos(dot(ecc2,r)/r_mag/e2);
    if dot(r,vf) < 0
        nu2 = 2*pi - nu2;
    end
    I2 = acos(dot(h,k)/h_mag);
    if I2 ~= 0 
        w2 = acos(dot(n_vec,ecc2)/e2/n_mag);
        if ecc2(3) < 0
            w2 = 2*pi - w2;
        end
        w2 = w2-pi/2;
        R2 = acos(dot(j,n_vec)/n_mag);
        if n_vec(2) < 0 
            R2 = 2*pi - R2;
        end
        R2 = R2-pi;
    end
    if I2 == 0
        w2 = atan2(ecc2(2),ecc2(1));
    end
end
end

function calcAuto(H)
%[a1,e1,i1,R1,w1,nu1,a2,e2,i2,R2,w2,nu2,sp,~,~,~,~,~,~,mu,R,col] = parameters(H);
end

function calcManual(H)
% read in all data
[a1,e1,i1,R1,w1,nu1,~,~,~,~,~,nu2,sp,dv1,dv2,imp1,imp2,th1,th2,mu,Radius,col] = parameters(H);

% calculate orbital parameters depending on user inputs
[aT,eT,iT,RT,wT,nuT] = calcs(a1,e1,i1,R1,w1,nu1,dv1,imp1,th1,mu);
[a2,e2,i2,R2,w2,nu3] = calcs(aT,eT,iT,RT,wT,nuT,dv2,imp2,th2,mu);

%fprintf('%2.2f  %2.2f  %2.2f  %2.2f  %2.2f  %2.2f\n',aT,eT,iT,RT,wT,nuT)
%fprintf('%2.2f  %2.2f  %2.2f  %2.2f  %2.2f  %2.2f\n',a2,e2,i2,R2,w2,nu3)

% reset text values for user display
altp2 = a2*(1-e2)-Radius;
alta2 = a2*(1+e2)-Radius;
set(H.altp2,'String',sprintf('%1.0f',altp2));
set(H.alta2,'String',sprintf('%1.0f',alta2));
set(H.i2,'String',sprintf('%1.0f',i2*180/pi));
set(H.R2,'String',sprintf('%1.0f',R2*180/pi));
set(H.w2,'String',sprintf('%1.0f',w2*180/pi));

% calculate the orbits
[x1,y1,z1,t1] = ellipse(a1,e1,i1,R1,w1,mu);
[xT,yT,zT,tT] = ellipse(aT,eT,iT,RT,wT,mu);
%zT = zT + rad1(3,1);
[x2,y2,z2,t2] = ellipse(a2,e2,i2,R2,w2,mu);
%z2 = z2 + rad1(3,1) + rad2(3,1);

% find index for animation
%{
r1 = a1*(1-e1^2)/(1+e1*cos(nu1));
ind1 = checkAngle(findIndex(x1,y1,z1,r1,nu1,e1),1,360);
r2 = aT*(1-eT^2)/(1+eT*cos(nuT));
ind2 = checkAngle(findIndex(xT,yT,zT,r2,nuT,eT),1,360);
r3 = aT*(1-eT^2)/(1+eT*cos(nu2));
ind3 = checkAngle(findIndex(xT,yT,zT,r3,nu2,eT),1,360);
r4 = a2*(1-e2^2)/(1+e2*cos(nu3));
ind4 = checkAngle(findIndex(x2,y2,z2,r4,nu3,e2),1,360);
%}
ind1 = calcIndex(a1,e1,nu1,mu,t1);
ind2 = calcIndex(aT,eT,nuT,mu,tT);
ind3 = calcIndex(aT,eT,nu2,mu,tT);
ind4 = calcIndex(a2,e2,nu3,mu,t2);
%{
% calculate the path of the satellite
if ind1 < 90 
    x1 = [x1,x1(1:ind1)]; y1 = [y1,y1(1:ind1)]; z1 = [z1,z1(1:ind1)];
else
    x1 = x1(1:ind1); y1 = y1(1:ind1); z1 = z1(1:ind1);
end
%}
if (ind2 < ind3) && (abs(ind3 - ind2) > 45)
    path = [x1,xT(ind2:ind3),x2(ind4:end),x2;...
        y1,yT(ind2:ind3),y2(ind4:end),y2;...
        z1,zT(ind2:ind3),z2(ind4:end),z2];
elseif (ind3 <= ind2) || (abs(ind3 - ind2) < 45)
    path = [x1,xT(ind2:end),xT(1:ind3),x2(ind4:end),x2;...
        y1,yT(ind2:end),yT(1:ind3),y2(ind4:end),y2;...
        z1,zT(ind2:end),zT(1:ind3),z2(ind4:end),z2];
end

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
plot3(H.ax,ax,zeros(size(ax)),zeros(size(ax)),'r',zeros(size(ax)),ax,zeros(size(ax)),'g',...
    zeros(size(ax)),zeros(size(ax)),ax,'b','LineWidth',2);
hold(H.ax,'on');
plot3(H.ax,x1,y1,z1,'w',xT,yT,zT,'y',x2,y2,z2,'c');
%plot3(H.ax,path(1,:),path(2,:),path(3,:),'w');
set(H.ax,'Color','no','xcolor','no','ycolor','no','zcolor','no')

% create the objects for the 
ind = 1;
[x,y,z] = sphere;
surf(H.ax,axlim*x/12,axlim*y/12,axlim*z/12,'EdgeColor','none','FaceColor',col);
sat = surf(H.ax,axlim*x/15+x1(ind),axlim*y/15+y1(ind),axlim*z/15+z1(ind),...
    'EdgeColor','none','FaceColor','g');
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
while get(H.animate,'Value') == 1
    if ind > ind1+(ind3-ind2)
        ind = ind + sp*n(3);
    elseif ind > ind1 && ind < ind1+(ind3-ind2)
        ind = ind + sp*n(2);
    elseif ind < ind1
        ind = ind + sp*n(1);
    end
    if ceil(ind) > length(path)
        ind = ind - length(x2);
    end
    set(sat,'XData',axlim*x/15+path(1,ceil(ind)),'YData',axlim*y/15+path(2,ceil(ind)),...
        'ZData',axlim*z/15+path(3,ceil(ind)));
    drawnow;
    pause(eps);
end
end

function ind = calcIndex(a,e,nu,mu,t)
E = acos((e+cos(nu))/(1+e*cos(nu)));
E = atan2(sin(E),cos(E));
time = (E-e*sin(E))/sqrt(mu/(a^3));
ind = 1;
while t(ind) < time && ind <= length(t)
    ind = ind + 1;
end
end

function [a1,e1,i1,R1,w1,nu1,a2,e2,i2,R2,w2,nu2,sp,dv1,dv2,imp1,imp2,th1,th2,mu,R,col] = parameters(H)
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
set(H.print_planet,'String',sprintf('\nmu = %2.2f (km^3)/(s^2)\n      R = %2.2f km',mu,R),'HorizontalAlignment','Center')

% check boundaries
set(H.altp1,'String',sprintf('%2.0f',checkBound(str2double(get(H.altp1,'String')),200,inf)));
set(H.alta1,'String',sprintf('%2.0f',checkBound(str2double(get(H.alta1,'String')),200,inf)));
set(H.altp2,'String',sprintf('%2.0f',checkBound(str2double(get(H.altp2,'String')),200,inf)));
set(H.alta2,'String',sprintf('%2.0f',checkBound(str2double(get(H.alta2,'String')),200,inf)));

% calculate all outputs for future use
a1 = (2*R+str2double(get(H.altp1,'String'))+str2double(get(H.alta1,'String')))/2;
e1 = 1-(R+str2double(get(H.altp1,'String')))/a1;
i1 = checkAngle(str2double(get(H.i1,'String')),0,360)*pi/180;
R1 = checkAngle(str2double(get(H.R1,'String')),0,360)*pi/180;
w1 = checkAngle(str2double(get(H.w1,'String')),0,360)*pi/180;
nu1 = checkAngle(str2double(get(H.nu1,'String')),0,360)*pi/180;
a2 = (2*R+str2double(get(H.altp2,'String'))+str2double(get(H.alta2,'String')))/2;
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
s = sprintf('a1 = %2.0f km  e1 = %2.0f  \na2 = %2.0f km  e2 = %2.0f  ',a1,e1,a2,e2);
set(H.parameters,'String',s);
end

function [x,y,z,t] = ellipse(a,e,I,R,w,mu)
% calculate the orbital ellipse by calculating eccentric anomaly and
% solving for time 
if e < .01 
    e = 0;
end
M = linspace(0,2*pi*sqrt((a^3)/mu),360);
E = M;
if e > 0
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
end
t = (E/sqrt((a^3)/mu)-e*sin(E/sqrt((a^3)/mu)))/sqrt(mu/(a^3));
path = ([a*cos(2*pi*t/t(end))-a*e;a*sqrt(1-e^2)*sin(2*pi*t/t(end));zeros(1,length(t))]);

% check if orbit is inclined, if not then set RAAN=0 (undefined for i=0)
if I == 0 
    R = 0;
end

% transform the ellipse depending on the angles
% direction cosine matrix to transform between frames
path = DCM323_vec(path,I,R,w);
x = path(1,:);
y = path(2,:);
z = path(3,:);
%x = [path(1,end/2:end),path(1,1:end/2)];
%y = [path(2,end/2:end),path(2,1:end/2)];
%z = [path(3,end/2:end),path(3,1:end/2)];
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
while var < low
    var = var + (high - low);
end
while var > high
    var = var - (high - low);
end
end

function display(H)
% read parameters to display the thrust elements to the user
[a,e,i,R,w,nu,~,~,~,~,~,~,~,dv,~,imp,~,inc,~,mu,PR,col] = parameters(H); 

% coordinates for orbit and sphere, and set axis limit
[x,y,z,~] = ellipse(a,e,i,R,w,mu);
[xs,ys,zs] = sphere;
axis_lim = linspace(-1.2*a*(1+e),1.2*a*(1+e),3);
%i = -i;
%R = -R;
%w = -w;

% save previous view
prevView = H.disp.View;

% predefine string for the legend and redefine if any of the angles = 0
str = {'orbital radius','delta v','true anomaly','inclination angle',...
    'impulse angle','reference direction'};

% radius as scalar and vector
r = a*(1-e^2)/(1+e*cos(nu))*[cos(nu);sin(nu);0];
r = DCM323_vec(r,i,R,w);
quiver3(H.disp,0,0,0,r(1),r(2),r(3),'-r','LineWidth',2,'AutoScale','off','MaxHeadSize',.5)
hold(H.disp,'on')

% delta v vector (not to scale, only to demonstrate the direction)
if dv ~= 0
    dv = .5*a*[-sin(nu-imp)*cos(inc);cos(nu-imp)*cos(inc);sin(inc)];
    dv = DCM323_vec(dv,i,R,w);
    dv = [dv(1,1),dv(2,1),dv(3,1)];
    quiver3(H.disp,r(1),r(2),r(3),dv(1),dv(2),dv(3),'-y','LineWidth',2,'AutoScaleFactor',1.5,'MaxHeadSize',.5)
else
    str = {str{1},str{3:end}};
end

% reference line for periapse
arg_line = [axis_lim(1),axis_lim(3);0,0;0,0];
arg_line = DCM323_vec(arg_line,i,R,w);

% reference line for thrust inclination angle
vert = [0,0;0,0;min(axis_lim),max(axis_lim)];
vert = DCM323_vec(vert,i,R,w);
vert = [vert(1,:)+r(1);vert(2,:)+r(2);vert(3,:)+r(3)];

% reference line for thrust impulse angle
perp = [0,-a*cos(nu-pi/2)/1.5;0,-a*sin(nu-pi/2)/1.5;0,0];
perp = DCM323_vec(perp,i,R,w);
perp = [perp(1,:)+r(1);perp(2,:)+r(2);perp(3,:)+r(3)];

% projection of the delta v vector in the orbital plane
if imp == 0 && inc == 0
    dv_proj = [0,0;0,0;0,0];
else
    dv_proj = DCM323_vec([0,.5*a;0,0;0,0],0,pi/2-imp,0);
    dv_proj = [dv_proj(1,:)+a*(1-e^2)/(1+e*cos(nu));dv_proj(2,:);zeros(size(dv_proj(3,:)))];
    dv_proj = DCM323_vec(dv_proj,i,R,w+nu);
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
    nu_curve = .5*a*(1-e)*[cos(th(1:floor(nu*180/pi/2)));sin(th(1:floor(nu*180/pi/2)));zeros(1,floor(nu*180/pi/2))];
    nu_curve = DCM323_vec(nu_curve,i,R,w);
    nu_arr = [nu_curve(1,end)-nu_curve(1,end-1),nu_curve(2,end)-nu_curve(2,end-1),...
        nu_curve(3,end)-nu_curve(3,end-1)];
    nu_ind = length(nu_curve)-1;
    quiver3(H.disp,nu_curve(1,nu_ind),nu_curve(2,nu_ind),nu_curve(3,nu_ind),...
        nu_arr(1),nu_arr(2),nu_arr(3),'-b','LineWidth',2,'AutoScale','off','MaxHeadSize',14);
    n = n + 1;
end

% inclination change angle arrow
if inc == 0
    str = {str{1:n-1},str{n+1:n+2}};
    inc_curve = [0;0;0];
else
    inc_curve = .35*a*[cos(th(1:floor(inc*180/pi/2)));zeros(1,floor(inc*180/pi/2));sin(th(1:floor(inc*180/pi/2)))];
    inc_curve = DCM323_vec(inc_curve,0,pi/2-imp,0);
    inc_curve = [inc_curve(1,:)+a*(1-e^2)/(1+e*cos(nu));inc_curve(2,:);inc_curve(3,:)];
    inc_curve = DCM323_vec(inc_curve,i,R,w+nu);
    inc_arr = [inc_curve(1,end)-inc_curve(1,end-1),inc_curve(2,end)-inc_curve(2,end-1),...
        inc_curve(3,end)-inc_curve(3,end-1)];
    inc_ind = length(inc_curve)-1;
    quiver3(H.disp,inc_curve(1,inc_ind),inc_curve(2,inc_ind),inc_curve(3,inc_ind),...
        inc_arr(1),inc_arr(2),inc_arr(3),'-g','LineWidth',2,'AutoScale','off','MaxHeadSize',.5);
    n = n + 1;
end

% impulse angle arrow
if imp == 0
    str = {str{1:n-1},str{n+1}};
    imp_curve = [0;0;0];
else
    imp_curve = .35*a*[-sin(th(1:floor(imp*180/pi/2)));cos(th(1:floor(imp*180/pi/2)));zeros(1,floor(imp*180/pi/2))];
    imp_curve = DCM323_vec(imp_curve,0,-imp+(3/180)*pi,0);
    imp_curve = [imp_curve(1,:)+a*(1-e^2)/(1+e*cos(nu));imp_curve(2,:);imp_curve(3,:)];
    imp_curve = DCM323_vec(imp_curve,i,R,w+nu);
    imp_arr = [imp_curve(1,1)-imp_curve(1,2);imp_curve(2,1)-imp_curve(2,2);...
        imp_curve(3,1)-imp_curve(3,2)];
    quiver3(H.disp,imp_curve(1,1),imp_curve(2,1),imp_curve(3,1),...
        imp_arr(1),imp_arr(2),imp_arr(3),'-m','LineWidth',2,'AutoScale','off','MaxHeadSize',.5);
end

% plot the curves and reference direction
plot3(H.disp,arg_line(1,:),arg_line(2,:),arg_line(3,:),'--w',...
    nu_curve(1,:),nu_curve(2,:),nu_curve(3,:),'b',...
    inc_curve(1,:),inc_curve(2,:),inc_curve(3,:),'g',...
    imp_curve(1,:),imp_curve(2,:),imp_curve(3,:),'m','LineWidth',2);

% plot the axes, horizontal and vertical reference lines for the thrust,
% the in plane projection of the delta v vector, and the orbit itself
plot3(H.disp,axis_lim,zeros(size(axis_lim)),zeros(size(axis_lim)),':w',...
    zeros(size(axis_lim)),axis_lim,zeros(size(axis_lim)),':w',...
    zeros(size(axis_lim)),zeros(size(axis_lim)),axis_lim,':w',...
    perp(1,:),perp(2,:),perp(3,:),'--w',...
    vert(1,:),vert(2,:),vert(3,:),'--w',...
    dv_proj(1,:),dv_proj(2,:),dv_proj(3,:),'--w',...
    x,y,z,'w');

% plot the central and orbiting bodies
surf(H.disp,3*(PR^(2/3))*xs,3*(PR^(2/3))*ys,3*(PR^(2/3))*zs,'EdgeColor','no','FaceColor',col);
surf(H.disp,2*(PR^(2/3))*xs+r(1),2*(PR^(2/3))*ys+r(2),2*(PR^(2/3))*zs+r(3),'EdgeColor','no','FaceColor','c');

% set axis background, limits, and legend, and allow rotation
set(H.disp,'Color','no','xcolor','no','ycolor','no','zcolor','no');
xlim([-1.2*a*(1+e) 1.2*a*(1+e)]);
axis(H.disp,'equal')
legend(H.disp,str,'Location','EastOutside','TextColor','w','Color','k','Box','on');
hold(H.disp,'off')

% hold the new view
H.disp.View = prevView;
end

function OUT = DCM323_vec(IN,I,R,w)
% this will be the direction cosine matrix for a typical orbital euler
% angle set (3-2-3), with the variable IN being a 3D column vector or a 
% 3xN matrix where N may be any positive integer
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