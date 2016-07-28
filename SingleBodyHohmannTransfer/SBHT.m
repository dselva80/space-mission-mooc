function SBHT
% set up display axes for user understanding of the angles
ops = {'BackgroundColor','ForegroundColor','FontWeight','FontSize'};
H.fig = figure('Name','Single_Body_Hohmann_Transfer','Position',[10 50 1350 600],'Color','k');
H.disp = axes('Units','Pixels','Position',[820 -150 600 600],...
    'xcolor','no','ycolor','no','zcolor','no','Color','no');
display(H.disp);
set(H.disp,'View',[-42,35]);

% set up the axes to plot orbits and animate
H.ax = axes('Units','Pixels','Position',[-50 0 700 700],...
    'xcolor','no','ycolor','no','zcolor','no','Color','no');
H.ax.View = [-100,20];

% commands
H.animPref_text = uicontrol('Style','text','Position',[600 575 150 20],'String','Animation Style',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.planets_text = uicontrol('Style','text','Position',[600 510 150 20],'String','Planet',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.text_speed = uicontrol('Style','text','Position',[600 340 150 20],'String','Animation Speed',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.animPref = uicontrol('Style','popupmenu','Position',[625 555 100 20],ops{1},'k',ops{2},'w',ops{3},...
    'Bold',ops{4},10,'Value',1,'String',{'Automatic','Manual'});
H.planets = uicontrol('Style','popupmenu','Position',[625 490 100 20],ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10,...
    'Value',3,'String',{'Mercury','Venus','Earth','Mars','Jupiter','Saturn','Uranus','Neptune','Moon'});
H.print_planet = uicontrol('Style','text','Position',[600 400 150 80],'String',...
    sprintf('\nmu = 398600.00 (km^3)/(s^2) \n      R = 6378.10 km'),ops{1},'none',ops{2},'w',ops{3},'Bold',ops{4},10);
H.parameters = uicontrol('Style','text','Position',[500 50 300 200],'String','',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.pref = uicontrol('Style','popupmenu','Position',[625 370 100 25],'String',{'Rotate on','Zoom on'},...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.speed = uicontrol('Style','slider','Position',[625 310 100 25],'Min',.25,'Max',4,'Value',2,...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.animate = uicontrol('Style','togglebutton','Position',[625 270 100 25],'Value',0,...
    'String','Play',ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);

% first orbit parameters
H.text_orb1 = uicontrol('Style','text','Position',[750 575 150 20],'String','-- Initial Orbit --',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.text_altp1 = uicontrol('Style','text','Position',[750 550 150 20],'String','alt at periapse (km)',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.text_alta1 = uicontrol('Style','text','Position',[750 500 150 20],'String','alt at apoapse (km)',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.text_i1 = uicontrol('Style','text','Position',[750 450 150 20],'String','inclination = [0,180]',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.text_R1 = uicontrol('Style','text','Position',[750 400 150 20],'String','RAAN = [0,360]',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.text_w1 = uicontrol('Style','text','Position',[760 335 130 40],'String','argument of periapse = [0,360]',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.altp1 = uicontrol('Style','edit','Position',[775 530 100 20],'String','200',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.alta1 = uicontrol('Style','edit','Position',[775 480 100 20],'String','200',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.i1 = uicontrol('Style','edit','Position',[775 430 100 20],'String','0',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.R1 = uicontrol('Style','edit','Position',[775 380 100 20],'String','0',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.w1 = uicontrol('Style','edit','Position',[775 315 100 20],'String','0',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);

% second orbit parameters
H.text_orb2 = uicontrol('Style','text','Position',[890 575 150 20],'String','-- Final Orbit --',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.text_altp2 = uicontrol('Style','text','Position',[890 550 150 20],'String','alt at periapse (km)',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.text_alta2 = uicontrol('Style','text','Position',[890 500 150 20],'String','alt at apoapse (km)',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.text_i2 = uicontrol('Style','text','Position',[890 450 150 20],'String','inclination = [0,180]',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.text_R2 = uicontrol('Style','text','Position',[890 400 150 20],'String','RAAN = [0,360]',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.text_w2 = uicontrol('Style','text','Position',[900 335 130 40],'String','argument of periapse = [0,360]',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.altp2 = uicontrol('Style','edit','Position',[915 530 100 20],'String','200',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.alta2 = uicontrol('Style','edit','Position',[915 480 100 20],'String','200',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.i2 = uicontrol('Style','edit','Position',[915 430 100 20],'String','0',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.R2 = uicontrol('Style','edit','Position',[915 380 100 20],'String','0',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.w2 = uicontrol('Style','edit','Position',[915 315 100 20],'String','0',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);

% manual commands (all visible off unless manual animation style selected)
% first burn
H.text_burn1 = uicontrol('Style','text','Position',[1040 575 150 20],'String','-- First Burn --',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10); 
H.text_nu1 = uicontrol('Style','text','Position',[1050 530 130 40],'String','true anomaly in orbit = [0,360]',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.text_dv1 = uicontrol('Style','text','Position',[1040 485 150 20],'String','delta v (km/s)',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.text_imp1 = uicontrol('Style','text','Position',[1050 415 130 40],'String','thrust impulse angle = [0,360]',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.text_th1 = uicontrol('Style','text','Position',[1040 345 150 40],'String','inclination change angle = [0,360]',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.nu1 = uicontrol('Style','edit','Position',[1065 510 100 20],'String','0',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.dv1 = uicontrol('Style','text','Position',[1065 465 100 20],'String','0',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.imp1 = uicontrol('Style','text','Position',[1065 395 100 20],'String','0',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.th1 = uicontrol('Style','text','Position',[1065 325 100 20],'String','0',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
% second burn
H.text_burn2 = uicontrol('Style','text','Position',[1180 575 150 20],'String','-- Second Burn --',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10); 
H.text_nu2 = uicontrol('Style','text','Position',[1190 530 130 40],'String','true anomaly in orbit = [0,360]',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.text_dv2 = uicontrol('Style','text','Position',[1180 485 150 20],'String','delta v (km/s)',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.text_imp2 = uicontrol('Style','text','Position',[1190 415 130 40],'String','thrust impulse angle = [0,360]',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.text_th2 = uicontrol('Style','text','Position',[1180 345 150 40],'String','inclination change angle = [0,360]',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.nu2 = uicontrol('Style','edit','Position',[1205 510 100 20],'String','0',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.dv2 = uicontrol('Style','text','Position',[1205 465 100 20],'String','0',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.imp2 = uicontrol('Style','text','Position',[1205 395 100 20],'String','0',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.th2 = uicontrol('Style','text','Position',[1205 325 100 20],'String','0',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);

% set callbacks
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
% set the correct visibilities depending on the selected animation style
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

calcManual(H)
end

function [a2,e2,i2,R2,w2,nu2] = calcs(a1,e1,i1,R1,w1,nu1,dv,imp,th,mu)
a2 = a1; e2 = e1; i2 = i1; R2 = R1; w2 = w1; nu2 = nu1;
if dv ~= 0
    % orbital radius at first burn
    r1 = a1*(1-e1^2)/(1+e1*cos(nu1));

    % use position derivative to define the intial velocity as vector
    dydx = sqrt(1-e1^2)*(r1*cos(nu1)/a1+e1)/sqrt(1-(r1*cos(nu1)/a1+e1)^2);
    v0 = [sqrt(mu*(2/r1-1/a1))*cos(atan(dydx)),sqrt(mu*(2/r1-1/a1))*sin(atan(dydx))];
    [v0(1),v0(2)] = checkVelocity(v0(1),v0(2),nu1);
    
    % define initial position as vector
    r0 = [r1*cos(abs(nu1)),r1*sin(abs(nu1))];
    if nu1 > 0 && nu1 < pi
        r0(2) = -r0(2);
    end
    if nu1 > 3*pi/2 && nu1 < pi/2
        r0(1) = -r0(1);
    end
    
    % define delta v as vector
    dv_vec = [dv*sin(nu1-imp)*sin(th),dv*cos(nu1-imp)*sin(th),dv*cos(th)];
    
    % define final velocity as vector
    beta = acos((r0(1)*v0(1)+r0(2)*v0(2))/sqrt(mu*(2/r1-1/a1))/r1);
    vf = [v0(1)+dv_vec(1),v0(2)+dv_vec(2),dv_vec(3)];
    vf2 = sqrt(dv^2+mu*(2/r1-1/a1)-2*sqrt(mu*(2/r1-1/a1))*dv*cos(pi/2-imp+beta));
    fprintf('%2.2f   %2.2f   %2.2f   %2.2f\n',vf2,vf(1),vf(2),vf(3))
    
    % calculate semi major axis
    a2 = 1./(2/r1-sqrt(vf(1)^2+vf(2)^2+vf(3)^2)./mu);
    
    % before calculating eccentricity, make sure the signs are correct
    [vf(1),vf(2)] = checkVelocity(vf(1),vf(2),nu1);

    % calculate eccentricity
    e2 = sqrt(1-((r1*cos(nu1)*vf(2)-r1*sin(nu1)*vf(1))^2)/mu/a2);
    if ~isreal(e2)
        e2 = 0;
    end
    
    % calculate true anomaly in second orbit
    nu2 = acos((a2*(1-e2^2)/r1-1)/e2);
    if ~isreal(nu2)
        nu2 = 0;
    end
    
    % find RAAN for the orbit
    if nu1 > 180 %&& nu2 > 180
        R2 = -R1 + nu1 - nu2;
    else
        R2 = R1 + nu1 - nu2;
    end

    % find argument of periapse for the orbit
    w2 = w1;
    i2 = i1;
end
end

function calcManual(H)
% read in all data
[a1,e1,i1,R1,w1,nu1,~,~,~,~,~,nu2,sp,dv1,dv2,imp1,imp2,th1,th2,mu,R,col] = read(H);

% calculate orbital parameters depending on user inputs
[aT,eT,iT,RT,wT,nuT] = findOrb(a1,e1,i1,R1,w1,nu1,dv1,imp1,th1,mu);
%[aT,eT,iT,RT,wT,nuT] = calcs(a1,e1,i1,R1,w1,nu1,dv1,imp1,th1,mu);
[a2,e2,i2,R2,w2,nu3] = findOrb(aT,eT,iT,RT,wT,nuT,dv2,imp2,th2,mu);

% reset text values for user display
altp2 = a2*(1-e2)-R;
alta2 = 2*(a2-R)-(a2*(1-e2)-R);
set(H.altp2,'String',sprintf('%2.0f',altp2));
set(H.alta2,'String',sprintf('%2.0f',alta2));
set(H.i2,'String',sprintf('%2.0f',i2));
set(H.R2,'String',sprintf('%2.0f',R2));
set(H.w2,'String',sprintf('%2.0f',w2));

% calculate the orbits
[x1,y1,z1] = ellipse(a1,e1,i1,R1,w1,mu);
[xT,yT,zT] = ellipse(aT,eT,iT,RT+pi,wT,mu);
[x2,y2,z2] = ellipse(a2,e2,i2,R2+pi,w2,mu);

% find index for animation
r1 = a1*(1-e1^2)/(1+e1*cos(nu1));
ind1 = checkAngle(findIndex(x1,y1,z1,r1,nu1),1,360);
r2 = aT*(1-eT^2)/(1+eT*cos(nuT));
ind2 = checkAngle(findIndex(xT,yT,zT,r2,nuT)+180,1,360);
r3 = aT*(1-eT^2)/(1+eT*cos(nu2));
ind3 = checkAngle(findIndex(xT,yT,zT,r3,nu2)+180,1,360);
r4 = a2*(1-e2^2)/(1+e2*cos(nu3));
ind4 = checkAngle(findIndex(x2,y2,z2,r4,nu3)+180,1,360);

%fprintf('%2.2f  %2.2f  %2.2f  %2.2f\n',ind1,ind2,ind3,ind4);

% calculate the path of the satellite
if ind1 < 90 
    x1 = [x1,x1(1,ind1)]; y1 = [y1,y1(1:ind1)]; z1 = [z1,z1(1:ind1)];
else
    x1 = x1(1:ind1); y1 = y1(1:ind1); z1 = z1(1:ind1);
end
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
axlim = 1.1*max(max(a2*(1+e2)),max(a1*(1+e1),aT*(1+eT)));
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
plot3(ax,zeros(size(ax)),zeros(size(ax)),'r',zeros(size(ax)),ax,zeros(size(ax)),'g',...
    zeros(size(ax)),zeros(size(ax)),ax,'b','LineWidth',2);
hold on;
%plot3(x1,y1,z1,'w',xT,yT,zT,'y',x2,y2,z2,'c');
plot3(path(1,:),path(2,:),path(3,:),'g');
set(gca,'Color','no','xcolor','no','ycolor','no','zcolor','no')
if get(H.pref,'Value') == 1
    rotate3d on;
elseif get(H.pref,'Value') == 2
    zoom on;
end

% create the objects for the 
ind = 1;
[x,y,z] = sphere;
surf(axlim*x/12,axlim*y/12,axlim*z/12,'EdgeColor','none','FaceColor',col);
sat = surf(axlim*x/15+x1(ind),axlim*y/15+y1(ind),axlim*z/15+z1(ind),...
    'EdgeColor','none','FaceColor','g');
hold off;

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

function i = findIndex(x,y,z,r,nu)
% find index for animation
r2 = sqrt(x.^2+y.^2+z.^2);
if nu > 180 
    i = 1;
    while r2(i) > r
        i = i + 1;
    end
else
    i = 181;
    while r2(i) < r
        i = i + 1;
        if i > length(x)
            i = i - length(x);
        end
    end
end
end

function [a2,e2,i2,R2,w2,nu2] = findOrb(a1,e1,i1,R1,w1,nu1,dv,imp,th,mu)
if dv == 0 
    a2 = a1; e2 = e1; i2 = i1; R2 = R1; w2 = w1; nu2 = nu1;
else
    % orbital radius at first burn
    r1 = a1*(1-e1^2)/(1+e1*cos(nu1));

    % new magnitude of velocity after first burn
    %vf = sqrt(mu*(2/r1-1/a1))+dv;
    dydx = sqrt(1-e1^2)*(r1*cos(nu1)/a1+e1)/sqrt(1-(r1*cos(nu1)/a1+e1)^2);
    v0x = sqrt(mu*(2/r1-1/a1))*cos(atan(dydx));
    v0y = sqrt(mu*(2/r1-1/a1))*sin(atan(dydx));
    [v0x,v0y] = checkVelocity(v0x,v0y,nu1);
    r0x = r1*cos(abs(nu1));
    r0y = r1*sin(abs(nu1));
    if nu1 > 0 && nu1 < pi
        r0y = -r0y;
    end
    if nu1 > 3*pi/2 && nu1 < pi/2
        r0x = -r0x;
    end
    beta = acos((r0x*v0x+r0y*v0y)/sqrt(mu*(2/r1-1/a1))/r1);
    vf = sqrt(dv^2+mu*(2/r1-1/a1)-2*sqrt(mu*(2/r1-1/a1))*dv*cos(pi/2-imp+beta));
    
    if nu1 < pi
        dvx = -dv*cos(pi/2+nu1-imp);
        dvy = -dv*sin(pi/2+nu1-imp);
    else
        dvx = -dv*cos(-pi/2+nu1-imp);
        dvy = -dv*sin(-pi/2+nu1-imp);
    end
    ang = acos((dvx*v0x+dvy*v0y)/dv/sqrt(mu*(2/r1-1/a1)));
    
    % semi major axis of transfer ellipse
    a2 = 1/(2/r1-(vf^2)/mu);

    % find the component velocities to calculate eccentricity
    vfx = vf*cos(atan(sqrt(1-e1^2)*(r1*cos(nu1)/a1+e1)/sqrt(1-(r1*cos(nu1)/a1+e1)^2)));
    vfy = vf*sin(atan(sqrt(1-e1^2)*(r1*cos(nu1)/a1+e1)/sqrt(1-(r1*cos(nu1)/a1+e1)^2)));

    % before calculating eccentricity, make sure the signs are correct
    [vfx,vfy] = checkVelocity(vfx,vfy,nu1);

    % calculate eccentricity
    e2 = sqrt(1-((r1*cos(nu1)*vfy-r1*sin(nu1)*vfx)^2)/mu/a2);
    if ~isreal(e2)
        e2 = 0;
    end

    % find the true anomaly of the second orbit at the first burn
    nu2 = acos((a2*(1-e2^2)/r1-1)/e2);
    % alternate form :nuT = acos((aT*(1-eT^2)*(new_vel*new_vel/mu+1/aT)/2-1)/eT);
    if ~isreal(nu2)
        nu2 = 0;
    end

    % find RAAN for the orbit
    if nu1 > 180 %&& nu2 > 180
        R2 = -R1 + nu1 - nu2;
    else
        R2 = R1 + nu1 - nu2;
    end

    % find argument of periapse for the orbit
    w2 = w1;
    i2 = i1;
end
end

function [a1,e1,i1,R1,w1,nu1,a2,e2,i2,R2,w2,nu2,sp,dv1,dv2,imp1,imp2,th1,th2,mu,R,col] = read(H)
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

% check that the automatic values are in an appropriate range
set(H.i1,'String',sprintf('%3.0f',i1*180/pi));
set(H.R1,'String',sprintf('%3.0f',R1*180/pi));
set(H.w1,'String',sprintf('%3.0f',w1*180/pi));
set(H.nu1,'String',sprintf('%3.0f',nu1*180/pi));
set(H.i2,'String',sprintf('%3.0f',i2*180/pi));
set(H.R2,'String',sprintf('%3.0f',R2*180/pi));
set(H.w2,'String',sprintf('%3.0f',w2*180/pi));
set(H.nu2,'String',sprintf('%3.0f',nu2*180/pi));

% check that the manual values are in an appropriate range
set(H.th1,'String',sprintf('%3.0f',th1*180/pi));
set(H.imp1,'String',sprintf('%3.0f',imp1*180/pi));
set(H.th2,'String',sprintf('%3.0f',th2*180/pi));
set(H.imp2,'String',sprintf('%3.0f',imp2*180/pi));
end

function [x,y,z] = ellipse(a,e,i,R,w,mu)
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
path = ([a*cos(2*pi*t/t(end))+a*e;a*sqrt(1-e^2)*sin(2*pi*t/t(end))]);

% transform the ellipse depending on the angles
dcm = [cos(R)*cos(w)*cos(i)-sin(R)*sin(w),-cos(R)*cos(i)*sin(w)-sin(R)*cos(w),-cos(R)*sin(i);...
    sin(R)*cos(w)*cos(i)+cos(R)*sin(w),-sin(R)*sin(w)*cos(i)+cos(R)*cos(w),-sin(R)*sin(i);...
    sin(i)*cos(w),-sin(i)*sin(w),cos(i)];
x = -(dcm(1,1)*path(1,:)+dcm(1,2)*path(2,:));
y = -(dcm(2,1)*path(1,:)+dcm(2,2)*path(2,:));
z = dcm(3,1)*path(1,:)+dcm(3,2)*path(2,:);
x = [x(end/2:end),x(1:end/2)];
y = [y(end/2:end),y(1:end/2)];
z = [z(end/2:end),z(1:end/2)];
end

function [vx,vy] = checkVelocity(vx,vy,nu)
if nu >= 0 && nu < pi/2
    if vx < 0
        vx = -vx;
    end
    if vy > 0 
        vy = -vy;
    end
elseif nu >= pi/2 && nu < pi
    if vx < 0
        vx = -vx;
    end
    if vy < 0 
        vy = -vy;
    end
elseif nu >= pi && nu < 3*pi/2
    if vx > 0
        vx = -vx;
    end
    if vy < 0 
        vy = -vy;
    end
elseif nu >= 3*pi/2
    if vx > 0
        vx = -vx;
    end
    if vy > 0 
        vy = -vy;
    end
end
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

function display(ax)
% create random ellitical orbit
[x,y,z] = ellipse(100,.5,0,0,0,1);
% calculate radius to display (true anomaly = 120)
r = [0,-50;0,86.6;0,0];
% vertical line through the position in the orbit
vert = [-50,-50;86.6,86.6;-100,100];
% horizaontal line perpendicular to the radius vector
perp = [-200,25;0,130;0,0];
% delta v vector
dv_arrow = [-50,-111.2;86.6,122;0,70.7];
dv_arrow_proj = [-50,-111.2;86.6,122;0,0];
dv_arrow_vert = [-111.2,-111.2;122,122;0,70.7];
dv_arrow_head = [-93.9,-111.2,-101.2;122,122,104.7;56.6,70.7,56.6];
% angle vectors (the arrowheads are all at 30 degree angles)
th = linspace(0,2*pi,360);
nu = [30*cos(th(1:120));30*sin(th(1:120));zeros(size(th(1:120)))];
nu_arrow_head = [-5,-15,-10;26,26,34.7;0,0,0];
% dcm1 = [cos(60),-sin(60),0;sin(60),cos(60),0;0,0,1]
imp = [50*cos(th(1:60))*cos(5*pi/6)-50*sin(th(1:60))*sin(5*pi/6)-50;...
    50*cos(th(1:60))*sin(5*pi/6)+50*sin(th(1:60))*cos(5*pi/6)+86.6;zeros(size(th(1:60)))];
imp_arrow_head = [-102,-93.3,-93.3;106.6,111.6,101.6;0,0,0];
% dcm2 = [1,0,0;0,cos(45),-sin(45);0,sin(45),cos(45)]
% dcm3 = [cos(60),-sin(60),0;sin(60),cos(60),0;0,0,1]
inc = [(50*cos(th(1:45))*cos(pi/4)-50*sin(th(1:45))*sin(pi/4))*-sin(pi/3)-50;...
    (50*cos(th(1:45))*cos(pi/4)-50*sin(th(1:45))*sin(pi/4))*cos(pi/3)+86.6;...
    50*cos(th(1:45))*sin(pi/4)+50*sin(th(1:45))*cos(pi/4)];
inc_arrow_head = [-78.1,-80.6,-72;100.4,104.7,99.7;44.2,35.4,35.4];

% x and y axis limits
limitx = linspace(-225,75,10);
limity = linspace(-100,100,10);
plot3(ax,dv_arrow(1,:),dv_arrow(2,:),dv_arrow(3,:),'y',...
    imp(1,:),imp(2,:),imp(3,:),'g',...
    inc(1,:),inc(2,:),inc(3,:),'b',...
    nu(1,:),nu(2,:),nu(3,:),'m',...
    r(1,:),r(2,:),r(3,:),'r',...
    x,y,z,'w',...
    vert(1,:),vert(2,:),vert(3,:),'--w',...
    perp(1,:),perp(2,:),perp(3,:),'--w',...
    dv_arrow_proj(1,:),dv_arrow_proj(2,:),dv_arrow_proj(3,:),'--w',...
    dv_arrow_vert(1,:),dv_arrow_vert(2,:),dv_arrow_vert(3,:),'--w',...
    limitx,zeros(size(limitx)),zeros(size(limitx)),':w',...
    zeros(size(limity)),limity,zeros(size(limity)),':w',...
    zeros(size(limity)),zeros(size(limity)),limity,':w');
hold on;
triang(dv_arrow_head(1,:),dv_arrow_head(2,:),dv_arrow_head(3,:),'y');
triang(imp_arrow_head(1,:),imp_arrow_head(2,:),imp_arrow_head(3,:),'g');
triang(inc_arrow_head(1,:),inc_arrow_head(2,:),inc_arrow_head(3,:),'b');
triang(nu_arrow_head(1,:),nu_arrow_head(2,:),nu_arrow_head(3,:),'m');
set(gca,'Color','no','xcolor','no','ycolor','no','zcolor','no');
[xs,ys,zs] = sphere;
surf(ax,13*xs,13*ys,13*zs,'EdgeColor','no','FaceColor','y');
surf(ax,5*xs-50,5*ys+86.6,5*zs,'EdgeColor','no','FaceColor','c');
axis equal;
hold off;
legend({'delta v','impulse angle','inclination angle','true anomaly','orbital radius'},...
    'Location',[.8033 .1 .0948 .1433],'TextColor','w','Color','k','Box','on');
end

function triang(x,y,z,col)
tri = delaunay(x,y);
trisurf(tri,x,y,z,'EdgeColor','no','FaceColor',col);
end