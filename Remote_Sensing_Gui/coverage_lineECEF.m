function [rmin,rmax,cmin,cmax] = coverage_lineECEF( circle_lat,circle_lon,r,A,R )
% returns min row 


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

     
    [row1,col1] = latlon2pix(R,lat1,y1);
    [row2,col2] = latlon2pix(R,lat2,y2);
    [row3,col3] = latlon2pix(R,x1,lon1);
    [row4,col4] = latlon2pix(R,x2,lon2);
    
    rows = [row1 row2 row3 row4];
    cols = [col3 col4 col2 col1];
    
    rmin = round(min(rows));
    rmax = round(max(rows));
    cmin = round(min(cols));
    cmax = round(max(cols));
    
   
        
        
    
     




%%%%%%%
