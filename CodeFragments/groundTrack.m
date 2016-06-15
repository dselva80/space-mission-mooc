function w = groundTrack(I,P,f)
%groundTrack - Draw ground track of inclined, circular orbit over one day
%      
%      INPUTS
%      I    Oribtal inclination (degrees)
%      P    Orbital period (sidereal days)
%      
%      OUTPUTS
%      w    Handle to map object in figure 1

if ~exist('I','var')
    I = 30*pi/180; %inclination
else
    I = I*pi/180;
end
if ~exist('P','var')
    P = 1; %period sidereal days
end
if ~exist('f','var')
    f = 1; %period sidereal days
end

Omega = 2*pi; %earth rotation rate (rad/day)
t = linspace(0,1,1000); %time (days)
M = 2*pi*t./P; % mean anomaly

I = pi/2-I;
phi = asin(cos(I).*sin(M));
lambda = atan2(cos(Omega.*t).*sin(I).*sin(M) - cos(M).*sin(Omega.*t),...
    (sin(Omega.*t).*sin(I).*sin(M) + cos(M).*cos(Omega.*t)));

figure(f)
clf
set(f,'Position',[66    94   800   500],'Color',[1,1,1]);

w = worldmap('World');
cdat = load('coast');
plotm(cdat.lat, cdat.long)
hold on
plotm(phi*180/pi,lambda*180/pi,'.r','LineWidth',2)
hold off