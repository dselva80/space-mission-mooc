function [rmin,rmax,cmin,cmax] = coverage_lineECEF( circle_lat,circle_lon,r,A,R )
% returns image A
% earth = referenceEllipsoid('earth');


%[xEast,yNorth] = geodetic2enu(lambda,lambda,0,0,0,0,earth);

%r = sqrt(xEast^2+yNorth^2);

% len = length(latg)-1;
% 
% rmin = zeros(len);
% rmax = zeros(len);
% cmin = zeros(len);
% cmax = zeros(len);
% rmin=0;
% rmax=0;
% cmin=0;
% cmax=0;
% for n =1: length(circle_lat)-1
%     %angle to latitude -> slope of groundtrack
%     theta = atand(-(long(n+1)-long(n)) / (latg(n+1) - latg(n)) );
%     
%     
%     % finding local x y z
%     y1 = -r*sind(theta)+latg(n);
%     x1 = -r*cosd(theta)+long(n);
%     y2 = r*sind(theta)+latg(n);
%     x2 = r*cosd(theta)+long(n);
    
    
    
%     [xECEF1, yECEF1, zECEF1] = geodetic2ecef(earth,y1, x1, 1e6);
% 
%     [xECEF2, yECEF2, zECEF2] = geodetic2ecef(earth,y2, x2, 1e6);
%     
%     x(n)=xECEF1;
%     y(n)=yECEF1;
%     z(n)=zECEF1;
%     
% %     x(len-n+1)=xECEF2;
% %     y(len-n+1)=yECEF2;
% %     z(len-n+1)=zECEF2;
%     x(len/2+n)=xECEF2;
%     y(len/2+n)=yECEF2;
%     z(len/2+n)=zECEF2;

    quart = round(length(circle_lon)/4);
    half = round(length(circle_lon)/2);
    tf = quart+half;
    
    y1 = circle_lon(1);
    lat1 = circle_lat(1);
    
    y2 = circle_lon(half);
    lat2 = circle_lat(half);
    
    x1 = circle_lat(quart);
    lon1 = circle_lon(quart);
    
    x2 =circle_lat(tf);
    lon2=circle_lon(tf);

    
%     % lat
%     [lat1,index1] = min(circle_lat);
%     y1 = circle_lon(index1);
%     [lat2,index2] = max(circle_lat);
%     y2 = circle_lon(index2);
%     % lon
%     [lon1,index3] = min(circle_lon);
%     x1 = circle_lat(index3);
%     [lon2,index4] = max(circle_lon);
%     x2 = circle_lat(index4);
    
    
    
    [row1,col1] = latlon2pix(R,lat1,y1);
    [row2,col2] = latlon2pix(R,lat2,y2);
    [row3,col3] = latlon2pix(R,x1,lon1);
    [row4,col4] = latlon2pix(R,x2,lon2);
    
    
    
    rmin = round(min(row1,row2));
    rmax = round(max(row1,row2));
    cmin = round(min(col3,col4));
    cmax = round(max(col3,col4));
    
   
        
        
    
     




%%%%%%%
