function SBHT
ops = {'BackgroundColor','ForegroundColor','FontWeight','FontSize'};
H.fig = figure('Name','Single_Body_Hohmann_Transfer','Position',[10 50 1350 600],'Color','k');
H.ax = axes('Units','Pixels','Position',[0 50 600 600],...
    'xcolor','none','ycolor','none','zcolor','none','Color','none');
H.ax.View = [-100,20];

% commands
H.planets_text = uicontrol('Style','text','Position',[750 575 200 20],'String','Select Planet',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.animPref_text = uicontrol('Style','text','Position',[600 575 150 20],'String','Select Animation Style',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.planets = uicontrol('Style','popupmenu','Position',[800 550 100 20],ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10,...
    'Value',1,'String',{'Mercury','Venus','Earth','Mars','Jupiter','Saturn','Uranus','Neptune','Moon'});
H.animPref = uicontrol('Style','popupmenu','Position',[625 550 100 20],ops{1},'k',ops{2},'w',ops{3},...
    'Bold',ops{4},10,'Value',1,'String',{'Automatic','Manual'});
H.print_planet = uicontrol('Style','text','Position',[925 555 400 20],'String',...
    'Mercury:    mu = 22032.00 (km^3)/(s^2)    R = 2439.70 km',ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.parameters = uicontrol('Style','text','Position',[600 50 700 60],'String','',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.pref = uicontrol('Style','popupmenu','Position',[750 200 100 25],'String',{'Rotate on','Zoom on'},...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.text_speed = uicontrol('Style','text','Position',[850 230 150 20],'String','Animation Speed',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.speed = uicontrol('Style','slider','Position',[875 200 100 25],'Min',.25,'Max',4,'Value',2,...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.animate = uicontrol('Style','togglebutton','Position',[1000 200 100 25],'Value',0,...
    'String','Play',ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);

% first orbit parameters
H.text_orb1 = uicontrol('Style','text','Position',[600 510 300 20],'String','-- First Orbit --',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.text_altp1 = uicontrol('Style','text','Position',[600 485 150 20],'String','alt at periapse (km)',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.text_alta1 = uicontrol('Style','text','Position',[750 485 150 20],'String','alt at apoapse (km)',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.text_i1 = uicontrol('Style','text','Position',[600 430 150 20],'String','inclination = [0,180]',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.text_R1 = uicontrol('Style','text','Position',[750 430 150 20],'String','RAAN = [0,360]',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.text_w1 = uicontrol('Style','text','Position',[610 355 130 40],'String','argument of periapse = [0,360]',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.text_nu1 = uicontrol('Style','text','Position',[750 355 150 40],'String','true anomaly at escape = [0,360]',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.altp1 = uicontrol('Style','edit','Position',[625 460 100 20],'String','200',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.alta1 = uicontrol('Style','edit','Position',[775 460 100 20],'String','200',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.i1 = uicontrol('Style','edit','Position',[625 405 100 20],'String','0',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.R1 = uicontrol('Style','edit','Position',[775 405 100 20],'String','0',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.w1 = uicontrol('Style','edit','Position',[625 335 100 20],'String','0',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.nu1 = uicontrol('Style','edit','Position',[775 335 100 20],'String','0',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);

% second orbit parameters
H.text_orb2 = uicontrol('Style','text','Position',[950 510 300 20],'String','-- Second Orbit --',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.text_altp2 = uicontrol('Style','text','Position',[950 485 150 20],'String','alt at periapse (km)',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.text_alta2 = uicontrol('Style','text','Position',[1100 485 150 20],'String','alt at apoapse (km)',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.text_i2 = uicontrol('Style','text','Position',[950 430 150 20],'String','inclination = [0,180]',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.text_R2 = uicontrol('Style','text','Position',[1100 430 150 20],'String','RAAN = [0,360]',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.text_w2 = uicontrol('Style','text','Position',[960 355 130 40],'String','argument of periapse = [0,360]',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.text_nu2 = uicontrol('Style','text','Position',[1100 355 150 40],'String','true anomaly at capture = [0,360]',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.altp2 = uicontrol('Style','edit','Position',[975 460 100 20],'String','200',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.alta2 = uicontrol('Style','edit','Position',[1125 460 100 20],'String','200',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.i2 = uicontrol('Style','edit','Position',[975 405 100 20],'String','0',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.R2 = uicontrol('Style','edit','Position',[1125 405 100 20],'String','0',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.w2 = uicontrol('Style','edit','Position',[975 335 100 20],'String','0',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);
H.nu2 = uicontrol('Style','edit','Position',[1125 335 100 20],'String','0',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10);

% manual commands (all visible off unless manual animation style selected)
% first burn
H.text_burn1 = uicontrol('Style','text','Position',[950 510 150 20],'String','-- First Burn --',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10,'Visible','off'); 
H.text_dv1 = uicontrol('Style','text','Position',[950 485 150 20],'String','delta v (km/s)',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10,'Visible','off');
H.text_th1 = uicontrol('Style','text','Position',[950 430 150 20],'String','theta = [0,360]',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10,'Visible','off');
H.text_imp1 = uicontrol('Style','text','Position',[960 355 130 40],'String','thrust impulse angle = [0,360]',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10,'Visible','off');
H.dv1 = uicontrol('Style','edit','Position',[975 460 100 20],'String','0',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10,'Visible','off');
H.th1 = uicontrol('Style','edit','Position',[975 405 100 20],'String','0',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10,'Visible','off');
H.imp1 = uicontrol('Style','edit','Position',[975 335 100 20],'String','0',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10,'Visible','off');
% second burn
H.text_burn2 = uicontrol('Style','text','Position',[1100 510 150 20],'String','-- Second Burn --',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10,'Visible','off'); 
H.text_dv2 = uicontrol('Style','text','Position',[1100 485 150 20],'String','delta v (km/s)',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10,'Visible','off');
H.text_th2 = uicontrol('Style','text','Position',[1100 430 150 20],'String','theta = [0,360]',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10,'Visible','off');
H.text_imp2 = uicontrol('Style','text','Position',[1110 355 130 40],'String','thrust impulse angle = [0,360]',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10,'Visible','off');
H.dv2 = uicontrol('Style','edit','Position',[1125 460 100 20],'String','0',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10,'Visible','off');
H.th2 = uicontrol('Style','edit','Position',[1125 405 100 20],'String','0',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10,'Visible','off');
H.imp2 = uicontrol('Style','edit','Position',[1125 335 100 20],'String','0',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10,'Visible','off');

% add clarifications for plane change
H.text_percent = uicontrol('Style','text','Position',[950 290 150 40],'String','percent delta v for plane change',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10,'Visible','off');
H.percent = uicontrol('Style','edit','Position',[975 270 100 20],'String','0',...
    ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10,'Visible','off');
H.text_PC = uicontrol('Style','text','Position',[1100 300 150 20],'String',...
    'plane change burn',ops{1},'k',ops{2},'w',ops{3},'Bold',ops{4},10,'Visible','off');
H.PC = uicontrol('Style','popupmenu','Position',[1125 275 100 20],ops{1},'k',ops{2},'w',ops{3},...
    'Bold',ops{4},10,'Value',1,'String',{'First Burn','Second Burn'},'Visible','off');

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
set(H.percent,'Callback',{@CB,H})
set(H.PC,'Callback',{@CB,H})

function CB(~,~,H)
% set the correct visibilities depending on the selected animation style
if get(H.animPref,'Value') == 1
    vis1 = 'on';
    vis2 = 'off';
else
    vis1 = 'off';
    vis2 = 'on';
end
set(H.text_orb2,'Visible',vis1)
set(H.text_altp2,'Visible',vis1)
set(H.text_alta2,'Visible',vis1)
set(H.text_i2,'Visible',vis1)
set(H.text_R2,'Visible',vis1)
set(H.text_w2,'Visible',vis1)
%set(H.text_nu2,'Visible',vis1)
set(H.altp2,'Visible',vis1)
set(H.alta2,'Visible',vis1)
set(H.i2,'Visible',vis1)
set(H.R2,'Visible',vis1)
set(H.w2,'Visible',vis1)
%set(H.nu2,'Visible',vis1)
set(H.text_burn1,'Visible',vis2)
set(H.text_burn2,'Visible',vis2)
set(H.text_dv1,'Visible',vis2)
set(H.text_th1,'Visible',vis2)
set(H.text_imp1,'Visible',vis2)
set(H.text_dv2,'Visible',vis2)
set(H.text_th2,'Visible',vis2)
set(H.text_imp2,'Visible',vis2)
set(H.dv1,'Visible',vis2)
set(H.th1,'Visible',vis2)
set(H.imp1,'Visible',vis2)
set(H.dv2,'Visible',vis2)
set(H.th2,'Visible',vis2)
set(H.imp2,'Visible',vis2)
set(H.text_percent,'Visible',vis2)
set(H.percent,'Visible',vis2)
set(H.text_PC,'Visible',vis2)
set(H.PC,'Visible',vis2)

% change position of true anomaly of capture text/edits
if get(H.animPref,'Value') == 2
    set(H.text_nu2,'Position',[675 285 150 40])
    set(H.nu2,'Position',[700 265 100 20])
end

calcManual(H)

function calcManual(H)
[a1,e1,i1,R1,w1,nu1,~,~,~,~,~,nu2,sp,dv1,dv2,th1,th2,imp1,imp2,percent,mu,~,~] = read(H);
if get(H.PC,'Value') == 1
    d = dv1;
    dv1 = d*(1-percent/100);
    dvi1 = d*percent/100;
    dvi2 = 0;
else
    d = dv2;
    dv2 = d*(1-percent/100);
    dvi1 = 0;
    dvi2 = d*percent/100;
end

% calculate orbital parameters depending on user inputs
[aT,eT,iT,RT,wT,nuT] = findOrb(a1,e1,i1,R1,w1,nu1,dv1,dvi1,mu);
[a2,e2,i2,R2,w2,nu3] = findOrb(aT,eT,iT,RT,wT,nuT,dv2,dvi2,mu);

% calculate the orbits
[x1,y1,z1] = ellipse(a1,e1,i1,R1,w1,mu);
[xT,yT,zT] = ellipse(aT,eT,iT,RT,wT,mu);
[x2,y2,z2] = ellipse(a2,e2,i2,R2,w2,mu);

% find index for animation
r1 = a1*(1-e1^2)/(1+e1*cos(nu1));
ind1 = findIndex(x1,y1,z1,r1,nu1);
r2 = aT*(1-eT^2)/(1+eT*cos(nuT));
ind2 = findIndex(xT,yT,zT,r2,nuT);
r3 = aT*(1-e1^2)/(1+e1*cos(nu2));
ind3 = findIndex(xT,yT,zT,r3,nu2);
r4 = a2*(1-e2^2)/(1+e2*cos(nu3));
ind4 = findIndex(x2,y2,z2,r4,nu3);

% calculate the path of the satellite
%nu1 = ceil(nu1*180/pi)+1;
%nuT = ceil(nuT*180/pi)+1;
%nu2 = ceil(nu2*180/pi)+1;
%nu3 = ceil(nu3*180/pi)+1;
%fprintf('%2.2f  %2.2f  %2.2f  %2.2f\n',nu1,nuT,nu2,nu3)
path = [x1(1:ind1),xT(ind2:ind3),x2(ind4:end),x2;...
    y1(1:ind1),yT(ind2:ind3),y2(ind4:end),y2;...
    z1(1:ind1),zT(ind2:ind3),z2(ind4:end),z2];

% create axes and axis limits
axlim = 1.1*max(max(a2*(1+e2)),max(a1*(1+e1),aT*(1+eT)));
ax = linspace(-axlim,axlim,10);

% save previous view
prevView = H.ax.View;

% plot orbits
plot3(ax,zeros(size(ax)),zeros(size(ax)),'r',zeros(size(ax)),ax,zeros(size(ax)),'g',...
    zeros(size(ax)),zeros(size(ax)),ax,'b','LineWidth',2);
hold on;
plot3(x1,y1,z1,'w',xT,yT,zT,'y',x2,y2,z2,'c');
%plot3(path(1,:),path(2,:),path(3,:),'w')
set(gca,'Color','no','xcolor','no','ycolor','no','zcolor','no')
if get(H.pref,'Value') == 1
    rotate3d on;
elseif get(H.pref,'Value') == 2
    zoom on;
end
hold off;

% apply previous view
H.ax2.View = prevView;


% iteration method to find e -- did not work
%{
if aT == r1
    eT = 0;
else 
    eT = .6;
    num = 0;
    error = 1;
    %val = ((e1^2-1)*cos(nu1)/(1+e1*cos(nu1))-e1)*sqrt(1-e^2)/sqrt(1-((e1^2-1)*cos(nu1)/(1+e1*cos(nu1))-a1*e1)^2);
    vT = acos((aT*(1-eT^2)*(1+e1*cos(nu1))/a1/(1-e1^2)-1)/eT);
    %disp(vT)
    while (abs(error) > eps(2*pi)) && (num < 1000)
        %eT = eT-(val+((1-eT^2)^1.5)*cos(vT)/(1+e1*cos(vT))/sqrt(1-((1-eT^2)*cos(vT)/(1+eT*cos(vT)))^2));
        %der = 2*(2*eT*cos(vT)/(1+eT*cos(vT))+(eT^2-1)*((cos(vT))^2)/(1+eT*cos(vT))^2-1)*((eT^2-1)*cos(vT)/(1+eT*cos(vT))-eT)...
         %   +2*(-sqrt(1-eT^2)*sin(vT)*cos(vT)/((1+eT*cos(vT))^2)+eT*sin(vT)/-sqrt(1-eT^2)/(1+eT*cos(vT)))*...
          %  (-sqrt(1-eT^2)*sin(vT)/(1+eT*cos(vT)));
        % the derivative assumes that d/d_eT(vT) is negligible, derivative of eq
        deriv = 2*((1-eT^2)/(1+eT*cos(vT)))*(-2*eT/(1+eT*cos(vT))+(1-eT^2)*sin(vT)/((1+eT*cos(vT))^2))+...
            2*(1-3*eT^2)*cos(vT)/(1+eT*cos(vT))+2*(eT-eT^3)*sin(vT)*cos(vT)/((1+eT*cos(vT))^2)+2*eT;
        % equation to find eT
        error = (((eT^2-1)*cos(vT)/(1+eT*cos(vT))-eT)^2+(sqrt(1-eT^2)*sin(vT)/(1+eT*cos(vT)))^2-1);
        eT = eT-error/deriv;
        disp(eT)
        vT = (aT*(1-eT^2)*(1+e1*cos(nu1))/a1/(1-e1^2)-1)/eT;
        %disp(vT)
        %error = abs(error);
        %error = abs(val-((1-eT^2)^1.5)*vT/(1+e1*vT)/sqrt(1-((1-eT^2)*vT/(1+eT*vT))^2));
        disp(error)
        num = num+1;
    end
end
if eT < 0 
    eT = 0;
end
vT = acos((aT*(1-eT^2)*(1+e1*cos(nu1))/a1/(1-e1^2)-1)/eT);
%}


function ind = findIndex(x,y,z,r,nu)
% find index for animation
ind = 1;
r2 = sqrt(x(ind)^2+y(ind)^2+z(ind)^2);
if nu > 180 
    ind = 181;
    while r2 > r
        r2 = sqrt(x(ind)^2+y(ind)^2+z(ind)^2);
        ind = ind +1;
    end
else
    ind = 1;
    while r2 < r
        r2 = sqrt(x(ind)^2+y(ind)^2+z(ind)^2);
        ind = ind +1;
    end
end

function [a2,e2,i2,R2,w2,nu2] = findOrb(a1,e1,i1,R1,w1,nu1,dv1,dvi,mu)
% orbital radius at first burn
r1 = a1*(1-e1^2)/(1+e1*cos(nu1));

% new magnitude of velocity after first burn
new_vel = sqrt(mu*(2/r1-1/a1))+dv1;

% semi major axis of transfer ellipse
a2 = 1/(2/r1-(new_vel^2)/mu);

% find the component velocities to calculate eccentricity
vx = new_vel*cos(atan(sqrt(1-e1^2)*(r1*cos(nu1)/a1+e1)/sqrt(1-(r1*cos(nu1)/a1+e1)^2)));
vy = new_vel*sin(atan(sqrt(1-e1^2)*(r1*cos(nu1)/a1+e1)/sqrt(1-(r1*cos(nu1)/a1+e1)^2)));

% before calculating eccentricity, make sure the signs are correct
if nu1 >= 0 && nu1 < pi/2
    if vx < 0
        vx = -vx;
    end
    if vy > 0 
        vy = -vy;
    end
elseif nu1 >= pi/2 && nu1 < pi
    if vx < 0
        vx = -vx;
    end
    if vy < 0 
        vy = -vy;
    end
elseif nu1 >= pi && nu1 < 3*pi/2
    if vx > 0
        vx = -vx;
    end
    if vy < 0 
        vy = -vy;
    end
elseif nu1 >= 3*pi/2
    if vx > 0
        vx = -vx;
    end
    if vy > 0 
        vy = -vy;
    end
end

% calculate eccentricity
e2 = sqrt(1-((r1*cos(nu1)*vy-r1*sin(nu1)*vx)^2)/mu/a2);
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
    R2 = -(R1 - nu1 + nu2);
else
    R2 = R1 - nu2 + nu1;
end

% find argument of periapse for the orbit
w2 = w1;
i2 = i1;

function [a1,e1,i1,R1,w1,nu1,a2,e2,i2,R2,w2,nu2,sp,dv1,dv2,th1,th2,imp1,imp2,percent,mu,R,col] = read(H)
% read in the data to retrieve the planet the user selected and define the
% standard gravitational parameter and color accordingly
if get(H.planets,'Value') == 1
    mu = 22032;         % standard gravitational parameter
    str = 'Mercury';    % string of planet to print
    col = [1 .5 .5];    % color of planet to display (optional)
    R = 2439.7;         % radius of planet
elseif get(H.planets,'Value') == 2
    mu = 324860; 
    str = 'Venus';
    col = [.9 .6 .1];
    R = 6051.8;
elseif get(H.planets,'Value') == 3
    mu = 398600; 
    str = 'Earth';
    col = [0 1 0];
    R = 6378.1;
elseif get(H.planets,'Value') == 4
    mu = 42828; 
    str = 'Mars';
    col = [1 0 0];
    R = 3396.2;
elseif get(H.planets,'Value') == 5
    mu = 1.26687*10^8; 
    str = 'Jupiter';
    col = [1 .5 .3];
    R = 71492;
elseif get(H.planets,'Value') == 6
    mu = 3.79311*10^7; 
    str = 'Saturn';
    col = [1 .8 0];
    R = 60268;
elseif get(H.planets,'Value') == 7
    mu = 5.79394*10^6; 
    str = 'Uranus';
    col = [.7 .9 .9];
    R = 25559;
elseif get(H.planets,'Value') == 8
    mu = 6.83653*10^6; 
    str = 'Neptune';
    col = [.7 .9 .9];
    R = 24764;
elseif get(H.planets,'Value') == 9
    mu = 4904.87; 
    str = 'the Moon';
    col = [1 1 1];
    R = 1738.1;
end

% update the planet parameters for user
s = [str,sprintf(':    mu = %2.2f (km^3)/(s^2)    R = %2.2f km',mu,R)];
set(H.print_planet,'String',s)

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
percent = checkBound(str2double(get(H.percent,'String')),0,100);

% check that the automatic values are in an appropriate range
set(H.i1,'String',sprintf('%1.0f',i1*180/pi));
set(H.R1,'String',sprintf('%1.0f',R1*180/pi));
set(H.w1,'String',sprintf('%1.0f',w1*180/pi));
set(H.nu1,'String',sprintf('%1.0f',nu1*180/pi));
set(H.i2,'String',sprintf('%1.0f',i2*180/pi));
set(H.R2,'String',sprintf('%1.0f',R2*180/pi));
set(H.w2,'String',sprintf('%1.0f',w2*180/pi));
set(H.nu2,'String',sprintf('%1.0f',nu2*180/pi));

% check that the manual values are in an appropriate range
set(H.th1,'String',sprintf('%1.0f',th1*180/pi));
set(H.imp1,'String',sprintf('%1.0f',imp1*180/pi));
set(H.th2,'String',sprintf('%1.0f',th2*180/pi));
set(H.imp2,'String',sprintf('%1.0f',imp2*180/pi));
set(H.percent,'String',sprintf('%1.0f',percent));

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

function var = checkBound(var,low,high)
% checks to make sure that the slider value is within the appropriate range
if var < low 
    var = low;
elseif var > high
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