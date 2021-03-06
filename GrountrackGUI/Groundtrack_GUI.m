function Groundtrack_GUI
% Allegra Moran
% 6/6/16
% This GUI is for showing the orbit and swath of a group of satellites 

   %  Create and then hide the UI as it is being constructed.
   close all;
   figwidth = 800;
   figheight = 650;
   f = figure('Visible','off','Position',[360,500,figwidth,figheight],...
       'Name','Orbit Design','NumberTitle','off',...
       'SizeChangedFcn',@resizeui);
   
   % global variables passed between uicontrols and functions
   global a e i RAAN w v num_P lat lon alt latg long P hsat...
           hanim_gg hanim_g hanim_orb dt circ_swath_globe circ_swath_flat...
            cone_anim diam diamtext r circ_swath animating;
       
   
   halfh = round(figheight/2);
   thirdw = round(figwidth/3);
   tenthw = round(figwidth/15);
   tenthh = round(figheight/15);
   
   %radius of earth
   Re = 6.371e6;
   
   %  Construct the components.
   err_string='';
   
   animating=0;
   
   % semimajor axis edit text field with string
   atext = uicontrol('Style','text','String','Semimajor axis [m]');
   ah = uicontrol('Style','edit','Callback',@a_input);
   ah.String = 7e6;
   a=7e6;
   function [] = a_input(H,~)
        aa = get(H,'string');
        a = str2double(aa);
        if a > 50e6 || a < 6.371e6
            err_string ='Semimajor axis is generally <= 42e6 m (GEO) and > 6.37e6 m (Earth radius)';
        else
            err_string='';
        end
        herror.String = err_string;
        if animating == 1
            hanim.Value=0;
            animbutton(hanim)
            hanim.Value=1;
            animbutton(hanim)
        end
    end

   % eccentricity edit text field with string
   etext = uicontrol('Style','text','String','Eccentricity');
   eh = uicontrol('Style','edit','Callback',@e_input);
   eh.String = 0;
   e=0;
   function [] = e_input(H,~)
        ee = get(H,'string');
        e = str2double(ee);
        if e<0 || e> 1
            err_string = ('Eccentricity should be a value >=0 and <=1');
        else
            err_string='';
        end
        herror.String = err_string;
        if animating == 1
            hanim.Value=0;
            animbutton(hanim)
            hanim.Value=1;
            animbutton(hanim)
        end
    end

   % inclination edit text field with string
   itext = uicontrol('Style','text','String','Inclination [deg]');
   inch = uicontrol('Style','edit','Callback',@inc_input);
   inch.String = 10;
   i = 10;
   function [] = inc_input(H,~)
        ii = get(H,'string');
        i = str2double(ii);
        if abs(i) > 180
            err_string = ('Inclination should be between -180 and 180');
        else
            err_string='';
        end
        herror.String = err_string;
        
        if animating == 1
            hanim.Value=0;
            animbutton(hanim)
            hanim.Value=1;
            animbutton(hanim)
        end
    end
   
   % Right ascension of ascending node
   RAANtext = uicontrol('Style','text','String',...
       'Right Ascension of Ascending Node [deg]');
   RAANh = uicontrol('Style','edit','Callback',@RAAN_input);
   RAANh.String = 0;
   RAAN = 0;
   function [] = RAAN_input(H,~)
        RAAN = get(H,'string');
        RAAN = str2double(RAAN);
        if  abs(RAAN) > 180
            err_string = ('RAAN should be between -180 and 180');
        else
            err_string='';
        end
        herror.String = err_string;
        if animating == 1
            hanim.Value=0;
            animbutton(hanim)
            hanim.Value=1;
            animbutton(hanim)
        end
   end

   % True anomoly edit text field
   vtext = uicontrol('Style','text','String','True Anomoly [deg]');
   vh = uicontrol('Style','edit','Callback',@v_input);
   vh.String = 0;
   v = 0;
   function [] = v_input(H,~)
        v = get(H,'string');
        v = str2double(v);
        if  abs(v) > 180
            err_string = ('v should be between -180 and 180');
        else
            err_string='';
        end
        herror.String = err_string;
        if animating == 1
            hanim.Value=0;
            animbutton(hanim)
            hanim.Value=1;
            animbutton(hanim)
        end
   end

   % argument of perigee edit text field
   wtext = uicontrol('Style','text','String','Argument of Perigee [deg]');
   wh = uicontrol('Style','edit','Callback',@w_input);
   wh.String = 0;
   w = 0;
   function [] = w_input(H,~)
        w = get(H,'string');
        w = str2double(w);
        if  abs(w) > 180
            err_string = ('w should be between -180 and 180');
        else
            err_string='';
        end
        herror.String = err_string;
        if animating == 1
            hanim.Value=0;
            animbutton(hanim)
            hanim.Value=1;
            animbutton(hanim)
        end
   end

   % number of periods displayed
   numPtext = uicontrol('Style','text','String','Number of Periods');
   numPh = uicontrol('Style','edit','Callback',@numP_input);
   numPh.String=1;
   num_P=1;
   function [] = numP_input(H,~)
        num_P = get(H,'string');
        num_P = str2double(num_P);
        if  rem(num_P,1) ~= 0 || num_P < 1
            err_string = ('Number of periods should be an integer >0');
        else
            err_string='';
        end
        herror.String = err_string;
        if animating == 1
            hanim.Value=0;
            animbutton(hanim)
            hanim.Value=1;
            animbutton(hanim)
        end
   end

    
   % slider to control time step
   slidertext = uicontrol('Style','text','String','Time Step [s]');
   slideredit = uicontrol('Style','edit','Callback',@slidertext_input);
   slideredit.String=60;
   dt=60;
   hslider = uicontrol('Style','slider','Callback',@slider_input);
   function [] = slidertext_input(H,~)
        dt = get(H,'string');
        dt = str2double(dt);
        if animating == 1
            hanim.Value=0;
            animbutton(hanim)
            hanim.Value=1;
            animbutton(hanim)
        end
   end
   function [] = slider_input(H,~)
        dt = 60+H.Value*600;
        slideredit.String = dt;
        if animating == 1
            hanim.Value=0;
            animbutton(hanim)
            hanim.Value=1;
            animbutton(hanim)
        end
   end

    % button to control the field of view angle. 
   FOVtext = uicontrol('Style','text','String','Field of View [deg]');
   FOVedit = uicontrol('Style','edit','Callback',@FOV_input);
   FOVedit.String=30;
   FOV=30;
   function [] = FOV_input(H,~)
        FOV = get(H,'string');
        FOV = str2double(FOV);
        
         % field of view is limited by geometry of earth
        max_angle = 2*tan(Re/a)*180/pi; 
        if FOV > max_angle
            % minus 1 causes it to not crash
            FOV = max_angle-1;
            FOVedit.String = num2str(FOV);
        end
        if animating == 1
            hanim.Value=0;
            animbutton(hanim)
            hanim.Value=1;
            animbutton(hanim)
        end
   end

    
   % text field empty until error
   herror = uicontrol('Style','text','String',err_string);
    
   % push button for calculations- shows static orbit and groundtrack
   hpush = uicontrol('Style','pushbutton','String',...
       'Calculate with parameters','Callback',@pushbutton);
   
   % push button for animation of orbit and groundtrack
   hanim = uicontrol('Style','togglebutton','String',...
       'Animate','Callback',@animbutton);
   
   

   
   done_that = 0;
   function [] = animbutton(hObject, ~, ~)
        button_state = hObject.Value;
    if button_state == hObject.Max
        hanim.String = 'Stop';
        % h represents satellite
        if been_here == 1
            delete(data_g);
            delete(data_s);
            delete(data_g2);
            delete(data_swath1);
            delete(data_swath2);
        end
       if done_that==1
           delete(hsat);
           delete(hanim_gg);
           delete(hanim_g);
           delete(hanim_orb);
           delete(circ_swath_globe);
           delete(circ_swath_flat);
           delete(circ_swath);
           delete(cone_anim);
           delete(diam);
           delete(diamtext);
       end

        
        % calulates altitude, longitude and latitude of satellite orbit
       [lat, lon, alt, P] = Sat_orbit( a,e,i,RAAN,v,w,10*num_P,0,dt);
       % calculate longitude, latitude of groundtrack. alt and P not used
       [latg, long] = Sat_orbit( a,e,i,RAAN,v,w,10*num_P,1,dt);
       
       [ circle_lat,circle_lon,r ] = swath( alt(1), FOV, latg(1), long(1));
        
       heading = zeros(length(latg)-1);
       for x=1:length(latg)-1
            slat = [latg(x) latg(x+1)];
            slon = [long(x) long(x+1)];
            [heading(x), ~] = legs(slat,slon,'gc');
       end
                
            
        
        % here i set up handles so i can delete them in the loop.
        
        axes(haxes)
        % set up satellite
        hsat = plot3m(lat(1), lon(1), alt(1),'Marker','o',...
          'MarkerFaceColor','k','MarkerSize',8);
        % set up groundtrack on flat map
        hanim_g = scatter(ha,0,0,[],'ro','filled','MarkerFaceColor','r');
        % set up ground track on globe       
        hanim_gg = plot3m(latg(1),long(1),0,'LineWidth',3,'Color','r');
        % plots animated orbit
        hanim_orb = plot3m(lat(1), lon(1),alt(1),'Color','m','LineWidth',3);
        % plots circular swath on  globe
        circ_swath_globe = plot3m(circle_lat,circle_lon,0,'m.');
        %plot circular swath on flat map
        circ_swath_flat = scatter(ha,0,0,[],'m.');
        
        
        mstruct = gcm(ha);
        globestruct = gcm(haxes);
        
       
        % set up for animation loop
        
        % t counts minutes of a day in seconds
        t=0:dt:round(10*num_P*P);
        del_N =  num_P*P/3600*15;
        x = del_N/360;
        % nodal precession doesnt perfectly close up groundtack on itself at
        % larger altitudes (higher P)
        % Introducing fudgefactor which should be improved. this ff is much smaller than
        % the calculate one. not sure why that is.
        ff = round(a/6e6);
        fudgefactor = -4+ff;
        tail_length = round((num_P*P+x*num_P*P)/dt)+fudgefactor;
        if num_P ==1
            tail_length = tail_length+5;
        end
        
        % rotation of earth in deg/sec
        w_earth = 360/(24*3600);
        
        % setting up rotation transform
        t1 = hgtransform('Parent',haxes);
        set(hgeo,'Parent',t1);
        set(hconts,'Parent',t1);
        set(hanim_gg,'Parent',t1);
        set(circ_swath_globe,'Parent',t1);

        % plot satellite orbit on globe invisibly so axis doesnt shift
        % during animation
        endy = round(P/dt);
        plot3m(lat(1:endy),lon(1:endy),alt(1:endy),...
           'LineStyle','none','Marker','none');

       % setting up swath
        axes(hswath)
        ang=0:0.1:2*pi+.1;
        circ_swath = plot(0,0,'m','LineStyle','none','LineWidth',4);
        diam = plot(0,0,'m','LineStyle','none','LineWidth',3);
        diamtext = text(0,0,'','Color','m','Visible','off','FontSize',15,...
        'FontWeight','bold');
        swath_text = uicontrol('Style','text',...
            'String','Satellite Swath Screen','Visible','off',...
            'FontSize',15);
        swath_text.Position = [thirdw+1.9*tenthw,tenthh*7,200,30];
        
        downward = 0;
        k=1;
        % loop of animation
        while k <=length(t) && hObject.Value==1
            k=k+1;
            animating=1;

            
            % n is total radians rotated by earth at time
            n = w_earth/180*pi*t(k);
            Txy = makehgtform('zrotate',n);
            set(t1,'Matrix',Txy)

            % if beginning
            if k<tail_length 
                start=1;
            % else if middle
            else
                start=start+1;
            end
            stop=k;
            
            % plots animated groundtrack on globe
            [gx, gy,gz] = mfwdtran(globestruct,latg(start:stop), long(start:stop),0);
            hanim_gg.XData = gx;
            hanim_gg.YData = gy;
            hanim_gg.ZData = gz;
            
            % plots animated groundtrack on flat world map
            [flatgx, flatgy] = mfwdtran(mstruct,latg(start:stop), long(start:stop));
            hanim_g.XData = flatgx;
            hanim_g.YData = flatgy;
            
            % plot circular swath on flat map
            [ circle_lat,circle_lon,r,locx ] = ...
                swath( alt(stop),FOV, latg(stop), long(stop));
            
            %circ_swath_flat = plot3m(circle_lat,circle_lon,'m.');
            [fx,fy] = mfwdtran(mstruct,circle_lat,circle_lon);
            circ_swath_flat.XData = fx;
            circ_swath_flat.YData = fy;
            
            % plots animated orbit
            [ox, oy, oz] = mfwdtran(globestruct,lat(start:stop), ...
                lon(start:stop),alt(start:stop));
            hanim_orb.XData = ox;
            hanim_orb.YData = oy;
            hanim_orb.ZData = oz;
            
            % plot satellite
            hsat.XData = ox(end);
            hsat.YData = oy(end);
            hsat.ZData = oz(end);
           
            % plot circular swath on globe
            [cgx,cgy, cgz] = mfwdtran(globestruct,circle_lat, circle_lon,0);
            circ_swath_globe.XData = cgx;
            circ_swath_globe.YData = cgy;
            circ_swath_globe.ZData = cgz;
            drawnow
            
            
            %%%%%%% flat map of swath %%%%%%%
            
            
            % rotate swath to match heading
            
            origin = [latg(stop) long(stop) -heading(stop)]; 
            setm(hswath,'Origin',origin,'FLatLimit',[-Inf r]);
            
            
            
            
            % getting the limits of the swath
            xlimits = xlim(hswath);
            
            
            lenf = (xlimits(2)-xlimits(1));
            x1 = xlimits(1);
            x2 = x1+lenf/2;
            x3 = xlimits(2);
            ylimits = ylim(hswath);
            y1 = ylimits(1);
            y2 = y1+lenf/2;
            y3 = ylimits(2);
            
            % setting the lat ticks and tick labels for swath
            hswath.XTick=[x1 x2 x3];
            hswath.YTick=[y1 y2 y3];
            lattext = strcat(num2str(round(latg(stop))),'{\circ}');
            lat1 = round(latg(stop)-r);
            if lat1<-90
                lat1 = lat1+180;
            end
            lat3 = round(latg(stop)+r);
            if lat3>90
                lat3 =180-lat3;
            end
            lat1text = strcat(num2str(lat1),'{\circ}');
            lat3text = strcat(num2str(lat3),'{\circ}');
            
            % setting the lon ticks and tick labels for swath
            lontext = strcat(num2str(round(long(stop))),'{\circ}');
            lon1 = round(long(stop)-r);
            if lon1<-180
                lon1 = lon1+360;
            end
            lon3 = round(long(stop)+r);
            if lon3>180
                lon3 =360-lon3;
            end
            lon1text = strcat(num2str(lon1),'{\circ}');
            lon3text = strcat(num2str(lon3),'{\circ}');
            hswath.YTickLabel = {lat1text,lattext,lat3text};        
            hswath.XTickLabel = {lon1text,lontext,lon3text}; 
            
            % plot circular swath on flat swath
            rad = x2-x1;
            % this fixes the fact that the frame of the ortho projection 
            % is not filled with pixels.
            % diy frame, slightly smaller than the matlab frame would be
            xp=(rad-1/r*rad)*cos(ang);
            yp=(rad-1/r*rad)*sin(ang);
            circ_swath.XData = x2+xp;
            circ_swath.YData = y2+yp;
            circ_swath.LineStyle = '-';
            
            
            % plot title of swath
            swath_text.Visible = 'on';
            % plot diameter and diameter text
            xdiam = x1+1/r*rad:.001:x3-1/r*rad;
            ydiam = ones(size(xdiam))*y2;
            diam.XData =xdiam;
            diam.YData = ydiam;
            diam.LineStyle='-';
            
            % updating diameter km text
            str = strcat(num2str(round(abs(locx))),' km');
            diamtext.Visible = 'on';
            diamtext.String = str;
            diamtext.Position = [x2-x1/8 y2-y1/8];
            
            if downward==1
                diamtext.Position = [x2-x1/8 y2+y1/8];
            end
            
            
        end
        
    end
   if button_state == hObject.Min
     hanim.String = 'Animate';
     animating=0;
   end
        done_that = 1;
   end

   been_here = 0;
   data_g = 0;
   data_s = 0;
   data_g2 =0;
   data_swath1=0;
   data_swath2=0;
   
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CALCULATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   
   % Push button calculates groundtrack and orbits and plots them
   function [] = pushbutton(~, ~, ~)
    if strcmp(err_string,'')
       % deletes the handles of previous orbit data
       if been_here == 1
           %stop animation
            delete(data_g);
            delete(data_s);
            delete(data_g2);
            delete(data_swath1);
            delete(data_swath2);
       end
       if done_that==1
           delete(hsat);
           delete(hanim_gg);
           delete(hanim_g);
           delete(hanim_orb);
           delete(circ_swath_globe);
           delete(circ_swath_flat);
           delete(circ_swath);
           delete(cone_anim);
           delete(diam);
           delete(diamtext);
       end
       
       
       % calulates altitude, longitude and latitude of satellite orbit
       [lat, lon, alt, P] = Sat_orbit( a,e,i,RAAN,v,w,10*num_P,0,dt);
       % calculate longitude, latitude of groundtrack. alt and P not used
       [latg, long] = Sat_orbit( a,e,i,RAAN,v,w,10*num_P,1,dt);
       
       %calculate ascending nodal precession due to earths rotation
       del_N =  num_P*P/3600*15;
       x = del_N/360;
       % nodal precession doesnt perfectly close up groundtack on itself at
       % larger altitudes (higher P)
       % Introducing fudgefactor which should be changed when I
       % understand that^ calculation better.
       % fudge factor is good to +- 2 in my experience
       fudgefactor = round(a/5e6);
       fudgefactor = fudgefactor^4-fudgefactor^3-3*fudgefactor^2+fudgefactor+3;
       endy = round((num_P*P+x*num_P*P)/dt)+fudgefactor;
       
       
      
       
       % plot groundtrack on flat map
       data_g = scatterm(ha,latg(1:endy),long(1:endy),[],'r','o','filled');
       
       
       % plot groundtrack on 1st globe
       axes(haxes)
       data_g2 = plot3m(latg(1:endy),long(1:endy),zeros(1,endy),'LineWidth',3,...
            'Color','r');
        

       % plot satellite orbit on globe
       data_s = plot3m(lat(1:endy),lon(1:endy),alt(1:endy),...
           'LineWidth',3,'Color','m');
       
       
       % plots yellow circle on flat world map to represent ground coverage
       axes(ha)
       [ ~,~,r ] = swath( alt(1), FOV, latg(1), long(1));
       
       
       %outputs two lines, one line on either side of coverage area
       [ lats1, lons1] = coverage_line( latg(1:endy),long(1:endy),r );
       
       data_swath1 = plotm(lats1,lons1,0,'y.','MarkerSize',8);
       
       axes(haxes)
       data_swath2 = plot3m(lats1,lons1,0,'y.','MarkerSize',8);
       
       been_here=1;
       
       
    end  
   end

%function for when window is resized
function resizeui(~,~)
           
       % Get figure width and height
       figwidth = f.Position(3);
       figheight = f.Position(4);
       
       atext.Position = [figwidth-150 figheight-25 100 20];
       ah.Position = [figwidth-150 figheight-50 100 20];
       etext.Position = [figwidth-150 figheight-75 100 20];
       eh.Position = [figwidth-150 figheight-100 100 20];
       itext.Position = [figwidth-150 figheight-125 100 20];
       inch.Position = [figwidth-150 figheight-150 100 20];
       RAANtext.Position = [figwidth-200 figheight-175 200 20];
       RAANh.Position = [figwidth-150 figheight-200 100 20];
       vtext.Position = [figwidth-150 figheight-225 100 20];
       vh.Position = [figwidth-150 figheight-250 100 20];
       wtext.Position = [figwidth-200 figheight-275 200 20];
       wh.Position = [figwidth-150 figheight-300 100 20];
       numPtext.Position = [figwidth-150 figheight - 325 100 20];
       numPh.Position = [figwidth-150 figheight - 350 100 20];
       slidertext.Position = [figwidth-150 figheight - 375 100 20];
       slideredit.Position = [figwidth-150 figheight - 400 100 20];
       hslider.Position = [figwidth-150 figheight - 425 100 20];
       FOVtext.Position = [figwidth-150 figheight - 450 100 20];
       FOVedit.Position = [figwidth-150 figheight - 475 100 20];
       
       hpush.Position = [figwidth-175 figheight-525 150 20];
       hanim.Position = [figwidth-175 figheight-575 150 20];
       herror.Position = [25 5 figwidth 20];
       
       halfh = round(figheight/2);
       thirdw = round(figwidth/3);
       tenthw = round(figwidth/15);
       tenthh = round(figheight/15);
       ha.Position = [tenthw,halfh+5,2*thirdw,halfh-5];
       
       haxes.Position = [tenthw,tenthh/2+5,thirdw,halfh-tenthh];
       
       hswath.Position = [thirdw+2*tenthh,tenthh+5,thirdw-tenthh,halfh-2*tenthh];
end
    

    load topo.mat topolegend topo;
    
    % flat world map
    ha = axes('Units','Pixels'); 
    worldmap world;
    geoshow(topo,topolegend,'DisplayType','texturemap');
   
   % globe on new axis
    haxes = axes('Units','Pixels'); 
    axis vis3d
    axesm ('globe');
    view(60,60)
    axis off
    
    % load surface of globe
    hgeo = geoshow(topo,topolegend,'DisplayType','texturemap');
    demcmap(topo)
    land = shaperead('landareas','UseGeoCoords',true);
    hconts = plotm([land.Lat],[land.Lon],'Color','black');
   
    % create axis 
    hswath = axes('Units','pixels');
    hswath = axesm('eqdazim','Origin',[0 0],'FLatLimit',[-Inf 0], ...
            'frame','off','grid','off','FedgeColor','m','Flinewidth',3);
    
    box(hswath,'off')
    set(hswath,'color','none');
    axis on
    axis vis3d
    geoshow(topo,topolegend,'DisplayType','texturemap');
    
    
   % Make the UI visible.
   f.Visible = 'on';
   

end