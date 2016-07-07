function [ circle_lat,circle_lon,lambda,x1 ] = swath( alt, eta, latg, long)
% Allegra Moran 6/20/16

% equations are from section 8.3.1 in SME: The New SMAD
% assumes spherical earth. 

%INPUTS
% alt is altidude above the earth's surface in earth radii
% eta is the nadir angle of the spacecraft [degrees]
% latg is a single latitude in degrees
% long is a single longitude in degrees

%OUTPUTS
% circle_lat is an array of latitudes [degrees]
% circle_lon is the corresponding array of longitude point [degrees]
% plotted together, they draw a circle on a map
% x and y are perpendicular local distances

earth = referenceEllipsoid('earth','m');

etar = deg2rad(eta);

Re = 6.378e6;

sin_rho = Re /(Re+alt*Re);

if abs( sin(etar)/sin_rho ) > 1
    
    lambda = 90 - rad2deg(asin(sin_rho));
else

    epsilon = rad2deg(acos(sin(etar)./ sin_rho));

    lambda = 90 - epsilon - eta;
end

[circle_lat,circle_lon] = scircle1(latg,long,lambda);

% hopefully this is the index 90 degrees from the first circle_lat
half = round(length(circle_lon)/2);

[arc,az] = distance(circle_lat(1),circle_lon(1),circle_lat(half),circle_lon(half));

x1 = nm2km(arc*60);
end

