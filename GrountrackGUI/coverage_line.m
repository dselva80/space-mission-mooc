function [ lats, lons] = coverage_line( latg,long,lambda,circle_lats,circle_lons )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
earth = referenceEllipsoid('earth');

[xEast,yNorth,zUp] = geodetic2enu(lambda,lambda,0,0,0,0,earth);

r = sqrt(xEast^2+yNorth^2);

len = 2*(length(latg)-1);

lats=zeros(1,len);
lons=zeros(1,len);


for n =1: length(long)-1
    %angle to latitude -> slope of groundtrack
    theta = atand((latg(n+1) - latg(n)) / (long(n+1)-long(n)));
    
    % finding local x y z
    y1 = r*cosd(theta);
    x1 = r*sind(theta);
    
    [xECEF1, yECEF1, zECEF1] = enu2ecef(x1, -y1, 0, latg(n), long(n), 0, earth);

    [xECEF2, yECEF2, zECEF2] = enu2ecef(-x1, y1, 0, latg(n), long(n), 0, earth);
    
    [lat1,lon1] = ecef2geodetic(earth,xECEF1,yECEF1,zECEF1);
    
    lats(n)=lat1;
    lons(n)=lon1;

    [lat2,lon2] = ecef2geodetic(earth,xECEF2,yECEF2,zECEF2);
%     lats(len/2+n)=lat2;
%     lons(len/2+n)=lon2;
    lats(len-n+1)=lat2;
    lons(len-n+1)=lon2;
    
end