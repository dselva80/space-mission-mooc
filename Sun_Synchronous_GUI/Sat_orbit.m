function [ lat,lon,z, P,vel ] = Sat_orbit( a,e,i,Omega,v,w,num_P,gtrack,dt,planet )
% Allegra Moran 6/17/16
%   Groundtrack takes in orbital parameters and returns two corresponding
%   vectors, latitude and longitude. Plot on a map to represent
%   groundtrack of an orbit.

% Equations are adapted from 
% "Satellites Orbits and Missions" by Michael Capderou

% INPUTS

% a = semimajor axis [m]
% e = eccentricity [0 to 1]
% i = inclinations [deg]
% Omega = RAAN = right ascension of the ascending node [deg]
% gtrack is a 0 or 1. if 1, we take into account the rotation of the earth
%  (for the case of computing lat and long for a groundtrack). 
% dt is time step in seconds- controls velocity of the plot
% planet represents position to sun
% 1 - Mercury
% 2 - Venus
% 3 - Earth
% 4 - Mars
% 5 - Jupiter
% 6 - Saturn
% 7 - Uranus
% 8 - Neptune
% 9 - Pluto

% OUTPUTS 

% latitude is in degrees -180 to 180
% longitude is in degrees from -180 to 180
% ^ meant to be used with map feature
% z is the altitude in Earth radii
% P is period of orbit in seconds

% convert inputs to radians
i = deg2rad(i);
Omega = deg2rad(Omega);
w = deg2rad(w);
v = deg2rad(v);

% CONSTANTS

% Earth grav constant [m^3/s^2]

if planet == 3
    %earth
    mu = 3.986E14;
    J2 = 0.00108263;
    R = 6.371E6;
    % rotation of the Earth [rad/s]
    rot = 7.2921e-5;
elseif planet == 1
    % Mercury
    mu = 2.2032e13;
    J2 = 0.00006;
    R = 2.439E6;
    rot = 1.24e-06;
elseif planet == 2
    % Venus
    mu = 3.24859E14;
    J2 = 0.000027;
    R = 6.051E6;
    rot = 2.99E-7;
elseif planet == 4
    % Mars
    mu = 4.282837e13;
    J2 = 0.001964;
    R = 3.396E6;
    rot = 7.088e-05;
elseif planet == 5
    % Jupiter
    mu = 1.26686534E17;
    J2 = 0.01475;
    R = 71.492E6;
    rot = 1.77E-4;
elseif planet == 6
    % Saturn
    mu = 3.7931187E16;
    J2 = 0.01645;
    R = 60.268E6;
    rot =1.63E-4;
elseif planet == 7
    % Uranus
    mu = 5.794e15;
    J2 = 0.012;
    R = 25.559E6;
    rot = -1.04E-4;
elseif planet ==8
    % Neptune
    mu = 6.809e15;
    J2 = 0.0004;
    R = 24.764E6;
    rot = 1.08E-4;
else
    %Pluto
    mu = 8.71E11;
    J2 = 0;
    R = 1.195E6;
    rot = -1.29E-5;
end


% ORBITAL VARIABLES CALCULATIONS

% Mean motion
n = sqrt(mu/a^3);

% Eccentric anomoly
E = 2*atan( sqrt((1-e)/(1+e)) * tan(v/2) );

% Initial mean anomoly
Mo = E - e*sin(E);

% period
P = 2*pi/n;

% time vector
tf = round(P*num_P);
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
    if gtrack == 1
        Omega_t(j) = Omega+(dot_omega-rot)*t(j);
    else
        Omega_t(j) = Omega+(dot_omega*t(j));
    end
   
end

% orbital radius vector
r = a.*(1-e^2)./(1+e.*cos(v_t));

vel = sqrt( mu*( 2./r - 1/a ));

% rotation matrix divided into rows
first = cos(Omega_t).*cos(w+v_t) - sin(Omega_t).*sin(w+v_t).*cos(i);
second = sin(Omega_t).*cos(w+v_t) + cos(Omega_t).*sin(w+v_t).*cos(i);
third = sin(w+v_t).*sin(i);


% conversion to latitude longitude from rotation matrix
lat = rad2deg(asin(third));
longi = rad2deg( acos(first./cos(deg2rad(lat))) );
lon = longi .* sign(second);

% returns 0 altitude for groundtrack. returns altidude in earth radii for an
% orbit
if gtrack ==0
    z=(r-R)./R;
else
    z = zeros(size(lat));
end

end