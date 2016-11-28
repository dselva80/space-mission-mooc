function Single_Body_Hohmann_Transfer

% clear display
clear 
close all
clc

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
    [625 350 100 25],'Min',1,'Max',10,'Value',7,...
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
    [1040 345 150 40],'String','burn inclination = [0,360]',...
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
    [1180 530 150 40],'String','true anomaly in orbit = [0,360]',...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10);
H.text_dv2 = uicontrol('Style','text','Position',...
    [1180 485 150 20],'String','delta v (km/s)',...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10);
H.text_imp2 = uicontrol('Style','text','Position',...
    [1190 415 130 40],'String','thrust impulse angle = [0,360]',...
    ops{1},'no',ops{2},'w',ops{3},'Bold',ops{4},10);
H.text_th2 = uicontrol('Style','text','Position',...
    [1180 345 150 40],'String','burn inclination = [0,360]',...
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
set(H.planets,'Callback',{@MainCallBack,H})
set(H.animPref,'Callback',{@MainCallBack,H})
set(H.speed,'Callback',{@MainCallBack,H})
set(H.animate,'Callback',{@MainCallBack,H})
set(H.altp1,'Callback',{@MainCallBack,H})
set(H.alta1,'Callback',{@MainCallBack,H})
set(H.i1,'Callback',{@MainCallBack,H})
set(H.R1,'Callback',{@MainCallBack,H})
set(H.w1,'Callback',{@MainCallBack,H})
set(H.nu1,'Callback',{@MainCallBack,H})
set(H.altp2,'Callback',{@MainCallBack,H})
set(H.alta2,'Callback',{@MainCallBack,H})
set(H.i2,'Callback',{@MainCallBack,H})
set(H.R2,'Callback',{@MainCallBack,H})
set(H.w2,'Callback',{@MainCallBack,H})
set(H.nu2,'Callback',{@MainCallBack,H})
set(H.dv1,'Callback',{@MainCallBack,H})
set(H.th1,'Callback',{@MainCallBack,H})
set(H.imp1,'Callback',{@MainCallBack,H})
set(H.dv2,'Callback',{@MainCallBack,H})
set(H.th2,'Callback',{@MainCallBack,H})
set(H.imp2,'Callback',{@MainCallBack,H})
end

function MainCallBack(~,~,H)
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

% if altitude at periapse is higher than that at apoapse then switch
if str2double(get(H.altp1,'String')) > str2double(get(H.alta1,'String')) 
    tmp = get(H.altp1,'String');
    set(H.altp1,'String',get(H.alta1,'String'))
    set(H.alta1,'String',tmp)
end
if str2double(get(H.altp2,'String')) > abs(str2double(get(H.alta2,'String'))) 
    tmp = get(H.altp2,'String');
    set(H.altp2,'String',get(H.alta2,'String'))
    set(H.alta2,'String',tmp)
end

% update the thrust/angle display
display(H)

% continue animation/calculation
if get(H.animPref,'Value') == 1
    calculate_Automatic(H)
elseif get(H.animPref,'Value') == 2
    calculate_Manual(H)
end
end

function calculate_Automatic(H)
end

function calculate_Manual(H)
% read in all data
[ai,ei,Ii,Ri,wi,nu1,~,~,~,~,~,nu2,sp,dv1,dv2,...
    imp1,imp2,th1,th2,mu,PR,col] = params(H);

% define the eccentricity as vector for initial orbit 
ei = [ei; 0; 0];
ei = DCM323(ei,Ii,Ri,wi);

% guard against imaginary numbers (discontinuity at nu ~ pi)
if (nu1 > pi*(1-.0090/pi) && nu1 < pi*(1+.0090/pi)) ||...
        (nu1 > -.0090/pi && nu1 < .0090/pi)
    nu1 = nu1+pi*.5/180;
end
if (nu2 > pi*(1-.0090/pi) && nu2 < pi*(1+.0090/pi)) ||...
        (nu2 > -.0090/pi && nu2 < .0090/pi)
    nu2 = nu2+pi*.5/180;
end

% calculate all the orbital elements
[~,~,ai,ei,Ii,Ri,wi,~,h1] = calcs(ai,ei,Ii,Ri,wi,nu1,0,0,0,mu);
[~,~,at,et,It,Rt,wt,nut,ht] = calcs(ai,ei,Ii,Ri,wi,nu1,dv1,imp1,th1,mu);
[~,~,af,ef,If,Rf,wf,nu3,h2] = calcs(at,et,It,Rt,wt,nu2,dv2,imp2,th2,mu);

% update the final and initial orbit texts for the user to view
s = [sprintf('a1 = %2.0f km  e1 = %2.4f  \n',ai,ei),...
    sprintf('a2 = %2.0f km  e2 = %2.4f  ',af,ef)];
set(H.parameters,'String',s);

% reset text values for user display
altp2 = af*(1-ef)-PR;
alta2 = af*(1+ef)-PR;
if altp2 < 160
    warning('altitude extremely low or negative')
end
set(H.altp2,'String',num2str(altp2));
set(H.alta2,'String',num2str(alta2));
set(H.i2,'String',num2str(If*180/pi));
set(H.R2,'String',num2str(Rf*180/pi));
set(H.w2,'String',num2str(wf*180/pi));
if -mu/2*af >= 0
    set(H.alta2,'Visible','off');
else
    set(H.alta2,'Visible','on');
end

% calculate initial period for comparison
P1 = 2*pi*sqrt((ai^3)/mu);

% calculate the orbits
[RA,~,nuA] = orbit(ai,ei,Ii,Ri,wi,h1,mu,max(ai,at),P1);
[RB,~,nuB] = orbit(at,et,It,Rt,wt,ht,mu,max(ai,at),P1);
[RC,~,nuC] = orbit(af,ef,If,Rf,wf,h2,mu,max(ai,at),P1);

% set limits on the second burn location if the transfer orbit is open
if et < 1
    TAL = sprintf('true anomaly in orbit = [0,360]');
else
    TAL = sprintf('true anomaly in orbit = [%2.4f,%2.4f]',...
        nuB(1)*180/pi,nuB(end)*180/pi);
end
set(H.text_nu2,'String',TAL);

% create axes and axis limits
axlim = 1.1*max(max(sqrt(RA(1,:).^2+RA(2,:).^2+RA(3,:).^2)),...
    max(max(sqrt(RB(1,:).^2+RB(2,:).^2+RB(3,:).^2)),...
    max(sqrt(RC(1,:).^2+RC(2,:).^2+RC(3,:).^2))));
ax = linspace(-axlim,axlim,10);

% save previous view
prevView = H.ax.View;

% plot orbits
p3 = plot3(H.ax,RC(1,:),RC(2,:),RC(3,:),'--w');
hold(H.ax,'on')
p2 = plot3(H.ax,RB(1,:),RB(2,:),RB(3,:),'--m');
p1 = plot3(H.ax,RA(1,:),RA(2,:),RA(3,:),'y');
legend({'final orbit','intermediate orbit','initial orbit'},...
    'Location','South','TextColor','w','color','k')
%plot3(H.ax,path(1,:),path(2,:),path(3,:),'w');
plot3(H.ax,ax,zeros(size(ax)),zeros(size(ax)),'r',...
    zeros(size(ax)),ax,zeros(size(ax)),'g',...
    zeros(size(ax)),zeros(size(ax)),ax,'b','LineWidth',2);
set(H.ax,'Color','no','xcolor','no','ycolor','no','zcolor','no')

% find index for animation
if nu1 >= pi
    [~,i] = min(nuA);
    nuA = [nuA(i:end),nuA(1:i-1)];
    ind1 = calcIndex(nuA,nu1);
    ind1 = length(RA)/2+ind1;
else
    ind1 = calcIndex(nuA,nu1);
end
if nut >= pi
    [~,i] = min(nuB);
    nuB = [nuB(i:end),nuB(1:i-1)];
    ind2 = calcIndex(nuB,nut);
    ind2 = length(RB)/2+ind2;
else
    ind2 = calcIndex(nuB,nut);
end
if nu2 >= pi 
    [~,i] = min(nuB);
    nuB = [nuB(i:end),nuB(1:i-1)];
    ind3 = calcIndex(nuB,nu2);
    ind3 = length(RB)/2+ind3;
else
    ind3 = calcIndex(nuB,nu2);
end
if nu3 >= pi && ef >= 1
    nu3 = nu3 - 2*pi;
    ind4 = calcIndex(nuC,nu3);
    ind4 = length(RC)-ind4;
elseif nu3 >= pi
    nu3 = -nu3 + 2*pi;
    ind4 = calcIndex(nuC,nu3);
    ind4 = length(RC)-ind4;
else
    ind4 = calcIndex(nuC,nu3);
end

% round and check index of inds to protect against error/warning messages
ind1 = round(ckAngle(ind1,1,length(nuA)));
ind2 = round(ckAngle(ind2,1,length(nuB)));
if (nu2 <= pi && nu1 >= pi) 
    ind3 = round(ckAngle(ind3-length(nuB)/2,1,length(nuB)));
else
    ind3 = round(ckAngle(ind3,1,length(nuB)));
end
ind4 = round(ckAngle(ind4,1,length(nuC)));

% calculate path
if ind1 > length(RA)/30
    RA = RA(:,1:ind1);
else
    RA = [RA,RA(:,1:ind1)];
end
if ind3 > ind2+length(RB)/20
    RB = RB(:,ind2:ind3);
else
    RB = [RB(:,ind2:end),RB(:,1:ind3)];
end
if ef < 1
    path = [RA,RB,RC(:,ind4:end),RC];
else
    path = [RA,RB,RC(:,ind4:end)];
end

% create the objects for the surfaces
ind = 1;
[x,y,z] = sphere;
surf(H.ax,PR*x/2,PR*y/2,PR*z/2,'EdgeColor','none','FaceColor',col);
sat = surf(H.ax,PR*x/8+path(1,ind),PR*y/8+path(2,ind),...
    PR*z/8+path(3,ind),'EdgeColor','none','FaceColor','c');
hold(H.ax,'off');

% apply previous view
H.ax2.View = prevView;

% reset the string of the toggle button depending on its current state
if get(H.animate,'Value') == 1
    set(H.animate,'String','Pause');
else
    set(H.animate,'String','Play');
end

% initialize variable to count for burns
count = 0;

% animate if selected 
while get(H.animate,'Value') == 1 && ceil(ind) < length(path)
    % increase the index for animate update
    ind = ind + sp;
    
    % check for the location of each burn to pause animation
    if ind > length(RA) && count == 0
        count = 1;
        pause(.7)
        set(p1,'LineStyle','--')
        set(p2,'LineStyle','-')
    elseif ind > length(RA)+length(RB) && count == 1
        count = 2;
        pause(.7)
        set(p2,'LineStyle','--')
        set(p3,'LineStyle','-')
    end
    
    % check if index is out of bounds and final orbit is closed, then loop,
    % else set index to last value and end animation
    if ceil(ind) >= length(path)-1 && ef < 1
        ind = ind - length(RC);
    elseif ceil(ind) >= length(path)
        ind = length(path);
    end
    
    % update data and draw
    set(sat,'XData',PR*x/10+path(1,ceil(ind)),...
        'YData',PR*y/10+path(2,ceil(ind)),...
        'ZData',PR*z/10+path(3,ceil(ind)));
    drawnow;
    pause(eps);
end

end

function [r,t,nu] = orbit(a,e,I,R,w,h,mu,a0,P0)
% make eccentricity vector scalar
e = norm(e);

% closed orbits
if e < 1 && e > 0
    % mean motion
    n = sqrt(mu/a^3);
    
    % period
    P = 2*pi/n;
    
    % mean anomaly
    nu_inf = 2*pi;
    M = linspace(0,nu_inf,round(P*360/P0));
    
    % Newton Raphson iteration
    num = 0;
    error = 1;
    E = M./(1-e);
    if e > 0
        inds = E > sqrt(6*(1-e)./e);
        E(inds) = (6*M(inds)./e).^(1/3);
        while (error > eps(2*pi)) && (num < 1000)
            E = E-(M-E+e.*sin(E))./(e.*cos(E)-1);
            error = max(abs(M-(E-e.*sin(E))));
            num = num+1;
        end
    end
    
    % time, true anomaly, and position array
    t = (E-e*sin(E))/n;
    nu = atan2(sin(E)*sqrt(1-e^2)./(1-e*cos(E)),(cos(E)-e)./(1-e*cos(E)));
    r = a*(1-e^2)./(1+e*cos(nu));
    r = [r.*cos(nu);r.*sin(nu);zeros(size(nu))];
    
% parabolic orbit
elseif e == 1
    % semi parameter
    p = dot(h,h)/mu;
    
    % mean motion
    n = 2*sqrt(mu/p^3);
    
    % parameterize endpoints
    nu_inf = acos(dot(h,h)/3/a0/mu-1);
    D_inf = tan(nu_inf/2);
    M_inf = D_inf+(1/3)*D_inf^3; 
    
    % calculate the time of flight
    P = 2*M_inf/n;
    
    % mean anomaly
    M = linspace(-M_inf,M_inf,round(P*360/P0));
    
    % Newton Raphson iteration
    D = M;
    num = 0;
    error = 1;
    while (error > eps(2*pi)) && (num < 1000)
        D = D-(M-D-(D.^3)/3)./(-1-D.^2);
        error = max(abs(M-D-(D.^3)/3));
        num = num+1;
    end
    
    % time and position array
    t = 2*(D+(1/3)*D.^3)/n;
    r = dot(h,h)*(1+D.^2)/2/mu;
    nu = atan2(D*p./r,(p-r)./r);
    r = [r.*cos(nu);r.*sin(nu);zeros(size(nu))];
    
% hyperbolic orbit
elseif e > 1
    % mean motion
    n = 2*sqrt(-mu/a^3);
    
    % parameterize endpoints
    nu_inf = acos((dot(h,h)/3/a0/mu-1)/norm(e));
    F_inf = 2*atanh(tan(nu_inf/2)*sqrt((e-1)/(1+e)));
    M_inf = e*sinh(F_inf)-F_inf;
    
    % calculate the time of flight
    P = 2*M_inf/n;
    
    % mean anomaly
    M = linspace(-M_inf,M_inf,round(P*360/P0));
    
    % Newton Raphson iteration
    F = M;
    num = 0;
    error = 1;
    while (error > eps(2*pi)) && (num < 1000)
        F = F+(M-e*sinh(F)+F)./(e*cosh(F)-1);
        error = max(abs(M-e*sinh(F)+F));
        num = num+1;
    end
    
    % time and position array
    t = (e*sinh(F)-F)/n;
    nu = atan2(-sinh(F)*sqrt(e^2-1)./(1-e*cosh(F)),...
        (cosh(F)-e)./(1-e*cosh(F)));
    r = a*(1-e^2)./(1+e*cos(nu));
    r = [r.*cos(nu);r.*sin(nu);zeros(size(nu))];
end

% rotate the position vector in 3D
r = DCM323(r,I,R,w);

end

function [r,vf,af,ef,If,Rf,wf,nuf,h] = calcs(a0,e0,I0,R0,w0,nu0,dv,imp,th,mu)
% predefine RAAN in the case of zero inclination, in which it is undefined
Rf = 0;

% direction cosine matrix to transform between frames
dcm1 = [cos(w0)*cos(I0)*cos(R0)-sin(w0)*sin(R0),...
    cos(w0)*cos(I0)*sin(R0)+sin(w0)*cos(R0),cos(w0)*sin(I0);...
    -sin(w0)*cos(I0)*cos(R0)-cos(w0)*sin(R0),...
    -sin(w0)*cos(I0)*sin(R0)+cos(w0)*cos(R0),-sin(w0)*sin(I0);...
    -sin(I0)*cos(R0),-sin(I0)*sin(R0),cos(I0)];
dcm1 = transpose(dcm1);

% all calculations initally set in orbital plane and then rotated into
% inertial frame using dcm1
% position - scalar and vector
e0 = norm(e0);
r_mag = a0*(1-e0^2)/(1+e0*cos(nu0));
r = r_mag*[cos(nu0);sin(nu0);0];

% tangent of the orbit at the point of the burn specified by nu0
if e0 < 1
    dydx = -sqrt(1-e0^2)*(r(1)+a0*e0)/sqrt(a0^2-(r(1)+a0*e0)^2);
else
    dydx = -sqrt(e0^2-1)*(r(1)+a0*e0)/sqrt((r(1)+a0*e0)^2-a0^2);
end
if nu0 >= pi
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
j = [0; 1; 0];
k = [0; 0; 1];

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
eccf = cross(vf,h)/mu-r/r_mag;
ef = norm(eccf);
nuf = acos(dot(eccf,r)/r_mag/ef);
if dot(r,vf) < 0
    nuf = 2*pi - nuf;
end
If = acos(dot(h,k)/h_mag);
if If ~= 0 
    wf = acos(dot(n_vec,eccf)/ef/n_mag);
    if eccf(3) < 0
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
    wf = atan2(eccf(2),eccf(1));
end
end

function ind = calcIndex(nu_orbit,nu_burn)
ind = 1;
nu_burn = ckAngle(nu_burn,-pi,pi);
while ind < length(nu_orbit) && ckAngle(nu_orbit(ind),-pi,pi) < nu_burn
    ind = ind + 1;
end
end

function [a1,e1,i1,R1,w1,nu1,a2,e2,i2,R2,w2,nu2,...
    sp,dv1,dv2,imp1,imp2,th1,th2,mu,R,col] = params(H)
% read in the data to retrieve the planet the user selected and define the
% standard gravitational parameter and color accordingly
if get(H.planets,'Value') == 1
    mu = 22032;         % standard gravitational parameter
    col = [1 .5 .5];    % color of planet to display 
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
set(H.altp1,'String',num2str(Bound(str2double(get(H.altp1,'String')),200,inf)));
set(H.alta1,'String',num2str(Bound(str2double(get(H.alta1,'String')),200,inf)));
set(H.altp2,'String',num2str(Bound(str2double(get(H.altp2,'String')),200,inf)));
set(H.alta2,'String',num2str(Bound(str2double(get(H.alta2,'String')),200,inf)));

% calculate all outputs for future use
a1 = (2*R+str2double(get(H.altp1,'String'))+...
    str2double(get(H.alta1,'String')))/2;
e1 = 1-(R+str2double(get(H.altp1,'String')))/a1;
i1 = ckAngle(str2double(get(H.i1,'String')),0,360)*pi/180;
R1 = ckAngle(str2double(get(H.R1,'String')),0,360)*pi/180;
w1 = ckAngle(str2double(get(H.w1,'String')),0,360)*pi/180;
nu1 = ckAngle(str2double(get(H.nu1,'String')),0,360)*pi/180;
a2 = (2*R+str2double(get(H.altp2,'String'))+...
    str2double(get(H.alta2,'String')))/2;
e2 = 1-(R+str2double(get(H.altp2,'String')))/a2;
i2 = ckAngle(str2double(get(H.i2,'String')),0,360)*pi/180;
R2 = ckAngle(str2double(get(H.R2,'String')),0,360)*pi/180;
w2 = ckAngle(str2double(get(H.w2,'String')),0,360)*pi/180;
nu2 = ckAngle(str2double(get(H.nu2,'String')),0,360)*pi/180;
sp = get(H.speed,'Value');
dv1 = abs(str2double(get(H.dv1,'String')));
th1 = ckAngle(str2double(get(H.th1,'String')),0,360)*pi/180;
imp1 = ckAngle(str2double(get(H.imp1,'String')),0,360)*pi/180;
dv2 = abs(str2double(get(H.dv2,'String')));
th2 = ckAngle(str2double(get(H.th2,'String')),0,360)*pi/180;
imp2 = ckAngle(str2double(get(H.imp2,'String')),0,360)*pi/180;

% check that the values are in an appropriate range
set(H.i1,'String',num2str(i1*180/pi))
set(H.R1,'String',num2str(R1*180/pi))
set(H.w1,'String',num2str(w1*180/pi))
set(H.dv1,'String',num2str(dv1))
set(H.nu1,'String',num2str(nu1*180/pi))
set(H.th1,'String',num2str(th1*180/pi))
set(H.imp1,'String',num2str(imp1*180/pi))
set(H.i2,'String',num2str(i2*180/pi))
set(H.R2,'String',num2str(R2*180/pi))
set(H.w2,'String',num2str(w2*180/pi))
set(H.dv2,'String',num2str(dv2))
set(H.nu2,'String',num2str(nu2*180/pi))
set(H.th2,'String',num2str(th2*180/pi))
set(H.imp2,'String',num2str(imp2*180/pi))
end

function var = Bound(var,low,high)
% checks to make sure that the slider value is within the appropriate range
if var < low 
    var = low;
elseif var > high
    var = high;
end
end

function var = ckAngle(var,low,high)
% checks to make sure that the angle value is within the appropriate range
for i = 1:length(var)
    while var(i) < low
        var(i) = var(i) + (high - low);
    end
    while var(i) >= high
        var(i) = var(i) - (high - low);
    end
end
end

function display(H)
% read parameters to display the thrust elements to the user
[a,e,i,R,w,nu,~,~,~,~,~,~,~,~,~,imp,~,inc,~,mu,PR,col] = params(H); 

% coordinates for orbit and sphere, and set axis limit
[r,~,~] = orbit(a,e,i,R,w,0,mu,a,2*pi*sqrt((a^3)/mu));
x = r(1,:);
y = r(2,:);
z = r(3,:);
[xs,ys,zs] = sphere;
axis_lim = linspace(-1.2*a*(1+e),1.2*a*(1+e),3);

% save previous view
prevView = H.disp.View;

% predefine string for the legend and redefine if any of the angles = 0
str = {'orbital radius','delta v','true anomaly','inclination angle',...
    'impulse angle'};

% radius as scalar and vector
r = a*(1-e^2)/(1+e*cos(nu))*[cos(nu);sin(nu);0];
r = DCM323(r,i,R,w);
quiver3(H.disp,0,0,0,r(1),r(2),r(3),'-r',...
    'LineWidth',2,'AutoScale','off','MaxHeadSize',.5)
hold(H.disp,'on')

% delta v vector (not to scale, only to demonstrate the direction)
dv = .5*a*[-sin(nu-imp)*cos(inc);cos(nu-imp)*cos(inc);sin(inc)];
dv = DCM323(dv,i,R,w);
dv = [dv(1,1),dv(2,1),dv(3,1)];
quiver3(H.disp,r(1),r(2),r(3),dv(1),dv(2),dv(3),'-y',...
    'LineWidth',2,'AutoScaleFactor',1.5,'MaxHeadSize',.5)

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
if inc <= pi/45
    str = {str{1:n-1},str{n+1:n+2}};
    inc_curve = [0;0;0];
else
    inc_curve = .35*a*[zeros(1,floor(inc*180/pi/2));...
        cos(th(1:floor(inc*180/pi/2)));...
        sin(th(1:floor(inc*180/pi/2)))];
    inc_curve = DCM323(inc_curve,0,0,-imp);
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
if imp <= pi/45
    str = {str{1:n-1},str{n+1}};
    imp_curve = [0;0;0];
else
    imp_curve = .35*a*[sin(th(1:floor(imp*180/pi/2)));...
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

% hold the new view
H.disp.View = prevView;

hold(H.disp,'off')
end

function OUT = DCM323(IN,I,R,w)
% this will be the direction cosine matrix for a typical orbital euler
% angle set (3-2-3), with the variable IN being a 3D column vector or a 
% 3xN matrix where N may be any positive integer
if I == 0
    R = 0;
end
rot1 = [cos(R), sin(R), 0; -sin(R), cos(R), 0; 0, 0, 1];
rot2 = [cos(I), 0, sin(I); 0, 1, 0; -sin(I), 0, cos(I)];
rot3 = [cos(w), sin(w), 0; -sin(w), cos(w), 0; 0, 0, 1];
dcm = rot3*rot2*rot1;
dcm = transpose(dcm);
OUT(1,:) = IN(1,:)*dcm(1,1)+IN(2,:)*dcm(1,2)+IN(3,:)*dcm(1,3);
OUT(2,:) = IN(1,:)*dcm(2,1)+IN(2,:)*dcm(2,2)+IN(3,:)*dcm(2,3);
OUT(3,:) = IN(1,:)*dcm(3,1)+IN(2,:)*dcm(3,2)+IN(3,:)*dcm(3,3);
end