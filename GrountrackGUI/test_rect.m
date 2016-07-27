

lat=40;
lon=40;
h = 6.2132e+05;
lat0 = 0;
lon0=0;
h0=0;
r=10;

earth = referenceEllipsoid('earth');


[xh,yh,zh] = geodetic2ecef(earth,lat,lon,h);




   % globe on new axis
    haxes = axes('Units','Pixels'); 
    axis vis3d
    axesm ('globe','Grid', 'on','Geoid',earth);
    view(60,60)
    axis off
    % load surface of globe
    load topo.mat topolegend topo;
    hgeo = geoshow(topo,topolegend,'DisplayType','texturemap');
    demcmap(topo)
    land = shaperead('landareas','UseGeoCoords',true);
    hconts = plotm([land.Lat],[land.Lon],'Color','black');
    
    cone_handle = surf(xsurf,ysurf,zsurf,'Facecolor','m','Facealpha',.5,'LineStyle','none','Visible','off');
while(1)   
    for dist=-r:.9:3*r
        actdist=dist;
        if actdist>r
            actdist=r-(dist-r);
        end
        [latout,lonout] = reckon(lat,lon,actdist,360);
        [circlat,circlon] = scircle1(latout,lonout,5);

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

        cone_handle.Visible='on';
        cone_handle.XData = xsurf;
        cone_handle.YData = ysurf;
        cone_handle.ZData = zsurf;
        drawnow
    end
end
% [ax,ay,az] = geodetic2ecef(earth,lat,lon,h);
% 
% t1 = hgtransform('Parent',haxes);
% set(cone_handle,'Parent',t1)
% 
% for rads = .1:.1:20*pi
%     M = makehgtform('axisrotate',[ax,ay,az],rads);
%     set(t1,'Matrix',M)
%     drawnow
%     
% end