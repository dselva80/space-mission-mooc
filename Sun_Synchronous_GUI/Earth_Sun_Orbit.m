function [ rotation_matrix] = Earth_Sun_Orbit( dt )
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

% OUTPUTS 

% latitude is in degrees -180 to 180
% longitude is in degrees from -180 to 180
% ^ meant to be used with map feature
% z is the altitude in Earth radii
% P is period of orbit in seconds

% Inputs in degrees
i = 7.155;
% Omega = 174.9;
w = 288.1;
v = 0;
a = 149.6E9;
e = 0.0167;

num_P=1;

% CONSTANTS

% Earth grav constant [m^3/s^2]
% mu = 3.986E14;

%Sun grav constant
mu = 1.327E20;

% rotation of the Earth [rad/s]
%rotE = 1.9910e-07;



%oblateness of the Earth
% J2 = 0.00108263;


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
%dot_omega = 3/2*(1-e^2)^2*n*J2*(Re/a)^2*cos(i);

% initialize time dependent vectors before theyre assigned in loop
M_t = zeros(size(t));
v_t = zeros(size(t));
Omega_t = zeros(size(t));



for j=1:length(t)
    
   
   % Mean anomoly [rad]
   M_t(j) = Mo + n*t(j);
   
   % True anomoly [rad]
   v_t(j) = M_t(j) + (2*e - e^3/4+e^5*e^5/96).*sin(M_t(j));
   
   
end

% orbital radius vector
r = a.*(1-e^2)./(1+e.*cos(v_t));

vel = sqrt( mu*( 2./r - 1/a ));

% rotation matrix divided into rows
first = cos(Omega_t).*cos(w+v_t) - sin(Omega_t).*sin(w+v_t).*cos(i);
second = sin(Omega_t).*cos(w+v_t) + cos(Omega_t).*sin(w+v_t).*cos(i);
third = sin(w+v_t).*sin(i);

rotation_matrix = [first;second;third];


% conversion to latitude longitude from rotation matrix
% lat = rad2deg(asin(third));
% longi = rad2deg( acos(first./cos(deg2rad(lat))) );
% lon = longi .* sign(second);



end

