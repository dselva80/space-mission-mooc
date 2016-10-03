
function [xsurf, ysurf, zsurf] = conical_animation(r,h,latg,long,lat,lon)
% This function takes 2 LAT/LON/ALT locations and draws a cone.
% Meant to be used as a satellite sensor FOV visual

% INPUTS
% lat [deg] = latitude where r=R of circular swath
% lon [deg] = corresponding^ longitude at radius of circular swath
% h [altitude m] = altiude of satellite in orbit
% lat0 [deg] = latitude where r=0 (groundtrack underneath sat) of circular swath
% lon0 [deg] = correaponding^ longitude where r=0 (groundtrack underneath sat) of circular swath
% h0 [altitudem] = altitude of groundtrackundersat (should be zero for this
% program)


    earth = referenceEllipsoid('earth');
    
    [circlat,circlon] = scircle1(latg,long,r);

    [xh,yh,zh] = geodetic2ecef(earth,lat,lon,h);

    [xcirc,ycirc,zcirc] = geodetic2ecef(earth,circlat,circlon,0);

    xsurf = zeros(2,length(xcirc));
    xsurf(1,:) = xcirc; 
    xsurf(2,:) = xh;

    ysurf = zeros(2,length(ycirc));
    ysurf(1,:) = ycirc;
    ysurf(2,:) = yh;

    zsurf = zeros(2,length(zcirc));
    zsurf(1,:) = zcirc;
    zsurf(2,:) = zh;

end

