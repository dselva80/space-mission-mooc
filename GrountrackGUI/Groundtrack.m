function [ latitude,longitude ] = Groundtrack( a,e,i,v,w,num_P)
% Allegra Moran 6/7/16
%   Groundtrack takes in orbital parameters and returns two corresponding
%   vectors, latitude and ongitude. Plotted on a map to represent
%   groundtrack of an orbit.

% Equations are adapted from 
% "Satellites Orbits and Missions" by Michael Capderou

% INPUTS

% a = semimajor axis [m]
% e = eccentricity [-1 to 1]
% i = inclinations [deg]
% Omega = RAAN = right ascension of the ascending node [deg]

% OUTPUTS 
% latitude is in degrees -180 to 180
% longitude is in degrees from -180 to 180
% ^ meant to be used with map feature

% convert input to radians

i = deg2rad(i);
v = deg2rad(v);
w = deg2rad(w);

% CONSTANTS

% Earth grav constant [m^3/s^2]
mu = 3.986E14;

% rotation of the Earth [rad/s]
rotE = 7.2921158553e-5;

% radius of the Earth [m]
Re = 6.378E6; 

%oblateness of the Earth
J2 = 0.00108263;


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
dt = 60;
tf = round(P*num_P);
t = 0:dt:tf ;

% RAAN
dot_omega = 3/2*(1-e^2)^2*n*J2*(Re/a)^2*cos(i);


M = zeros(size(t));
v_t = zeros(size(t));
Omega_t = zeros(size(t));

for j=1:length(t)
    
   % Mean anomoly [rad]
   M(j) = Mo + n*t(j);
   
   % True anomoly [rad]
   v_t(j) = M(j) + (2*e - e^3/4+e^5*e^5/96)*sin(M(j));
   
   % time dependent RAAN
   Omega_t(j) = (dot_omega-rotE)*t(j);
   
end


% rotation matrix divided into rows
first = cos(Omega_t).*cos(w+v_t) - sin(Omega_t).*sin(w+v_t).*cos(i);
second = sin(Omega_t).*cos(w+v_t) + cos(Omega_t).*sin(w+v_t).*cos(i);
third = sin(w+v_t).*sin(i);


% converted to arctan needlessly :(
latitude = rad2deg(asin(third));
% arccos x = 1 / ((ATN x)^2 + 1)^.5
% longitude = rad2deg(acos(first./cos(deg2rad(latitude))));
x = first./cos(deg2rad(latitude));
longi = rad2deg(pi/2 - atan(x./sqrt(1-x.^2)));
longitude = longi .* sign(second);

end

