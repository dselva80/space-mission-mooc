function Groundtrack_GUI
% Allegra Moran
% 6/6/16
% This GUI is for showing the orbit and swath of a group of satellites 

   %  Create and then hide the UI as it is being constructed.
   close all;
   figwidth = 800;
   figheight = 650;
   f = figure('Visible','off','Position',[360,500,figwidth,figheight],...
       'Name','Orbit and Constellation Design','NumberTitle','off',...
       'SizeChangedFcn',@resizeui);
   
   % global variables passed between uicontrols and functions
   global a e i RAAN w v num_P lat lon alt latg long P hsat...
           hanim_gg hanim_g hanim_orb dt circ_swath_globe circ_swath_flat...
           circ_swath cone_anim cone_viz diam diamtext;
       % A R
       
   
   halfh = round(figheight/2);
   thirdw = round(figwidth/3);
   tenthw = round(figwidth/15);
   tenthh = round(figheight/15);
   
   %  Construct the components.
   err_string='';
   
   stop_sign=0;
   
   swath_text = uicontrol('Style','text',...
       'String','Satellite Swath Screen','Visible','off');
   
   % semimajor axis edit text field with string
   atext = uicontrol('Style','text','String','Semimajor axis [m]');
   ah = uicontrol('Style','edit','Callback',@a_input);
   ah.String = 7e6;
   a=7e6;
   function [] = a_input(H,E)
        aa = get(H,'string');
        a = str2double(aa);
        if a > 50e6 || a < 6.371e6
            err_string ='Semimajor axis is generally <= 42e6 m (GEO) and > 6.37e6 m (Earth radius)';
        else
            err_string='';
        end
        herror.String = err_string;
    end

   % eccentricity edit text field with string
   etext = uicontrol('Style','text','String','Eccentricity');
   eh = uicontrol('Style','edit','Callback',@e_input);
   eh.String = 0;
   e=0;
   function [] = e_input(H,E)
        ee = get(H,'string');
        e = str2double(ee);
        if e<0 || e> 1
            err_string = ('Eccentricity should be a value >=0 and <=1');
        else
            err_string='';
        end
        herror.String = err_string;
    end

   % inclination edit text field with string
   itext = uicontrol('Style','text','String','Inclination [deg]');
   inch = uicontrol('Style','edit','Callback',@inc_input);
   inch.String = 10;
   i = 10;
   function [] = inc_input(H,E)
        ii = get(H,'string');
        i = str2double(ii);
        if abs(i) > 180
            err_string = ('Inclination should be between -180 and 180');
        else
            err_string='';
        end
        herror.String = err_string;
    end
   
   % Right ascension of ascending node
   RAANtext = uicontrol('Style','text','String',...
       'Right Ascension of Ascending Node [deg]');
   RAANh = uicontrol('Style','edit','Callback',@RAAN_input);
   RAANh.String = 0;
   RAAN = 0;
   function [] = RAAN_input(H,E)
        RAAN = get(H,'string');
        RAAN = str2double(RAAN);
        if  abs(RAAN) > 180
            err_string = ('RAAN should be between -180 and 180');
        else
            err_string='';
        end
        herror.String = err_string;
   end

   % True anomoly edit text field
   vtext = uicontrol('Style','text','String','True Anomoly [deg]');
   vh = uicontrol('Style','edit','Callback',@v_input);
   vh.String = 0;
   v = 0;
   function [] = v_input(H,E)
        v = get(H,'string');
        v = str2double(v);
        if  abs(v) > 180
            err_string = ('v should be between -180 and 180');
        else
            err_string='';
        end
        herror.String = err_string;
   end

   % argument of perigee edit text field
   wtext = uicontrol('Style','text','String','Argument of Perigee [deg]');
   wh = uicontrol('Style','edit','Callback',@w_input);
   wh.String = 0;
   w = 0;
   function [] = w_input(H,E)
        w = get(H,'string');
        w = str2double(w);
        if  abs(w) > 180
            err_string = ('w should be between -180 and 180');
        else
            err_string='';
        end
        herror.String = err_string;
   end

   % number of periods displayed
   numPtext = uicontrol('Style','text','String','Number of Periods');
   numPh = uicontrol('Style','edit','Callback',@numP_input);
   numPh.String=1;
   num_P=1;
   function [] = numP_input(H,E)
        num_P = get(H,'string');
        num_P = str2double(num_P);
        if  rem(num_P,1) ~= 0 || num_P < 1
            err_string = ('Number of periods should be an integer >0');
        else
            err_string='';
        end
        herror.String = err_string;
   end

    
   % slider to control time step
   slidertext = uicontrol('Style','text','String','Time Step [s]');
   slideredit = uicontrol('Style','edit','Callback',@slidertext_input);
   slideredit.String=60;
   dt=60;
   hslider = uicontrol('Style','slider','Callback',@slider_input);
   function [] = slidertext_input(H,E)
        dt = get(H,'string');
        dt = str2double(dt);
   end
   function [] = slider_input(H,E)
        dt = 60+H.Value*600;
        slideredit.String = dt;
   end

    % button to control the field of view angle. 
   FOVtext = uicontrol('Style','text','String','Field of View [deg]');
   FOVedit = uicontrol('Style','edit','Callback',@FOV_input);
   FOVedit.String=30;
   FOV=30;
   function [] = FOV_input(H,E)
        FOV = get(H,'string');
        FOV = str2double(FOV);
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
   function [] = animbutton(hObject, eventdata, handles)
        button_state = hObject.Value;
    if button_state == hObject.Max
        hanim.String = 'Stop';
        swath_running=0;
       
        % h represents satellite
        if been_here == 1
            delete(data_g);
            delete(data_s);
            delete(data_g2);
            delete(data_swath1);
            delete(data_swath2);
            delete(cone_viz);
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
        
        
        % here i set up handles so i can delete them in the loop.
        
        axes(haxes)
        % set up satellite
        hsat = plot3m(lat(1), lon(1), alt(1),'Marker','o',...
          'MarkerFaceColor','k');
        % set up groundtrack on flat map
        hanim_g = scatterm(ha,latg(1),long(1),[],'r','o','filled');
        % set up ground track on globe       
        hanim_gg = plot3m(latg(1),long(1),0,'LineWidth',3,'Color','r');
        % plots animated orbit
        hanim_orb = plot3m(lat(1), lon(1),alt(1),'Color','m');
        % plots circular swath on  globe
        circ_swath_globe = plot3m(circle_lat,circle_lon,0,'y--');
        %plot circular swath on flat map
        circ_swath_flat = scatterm(ha,circle_lat,circle_lon,[],'y','.');
        
       
        % set up for animation loop
        
        % t counts minutes of a day in seconds
        t=0:dt:round(10*num_P*P);
        tail_length = round(num_P*P/dt);
        
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
        inv = plot3m(lat(1:endy),lon(1:endy),alt(1:endy),...
           'LineStyle','none','Marker','none');

       % reset axis limits
%        lat2 = max(inv.XData);
%        lon2 = max(inv.YData);
%        lat1 = min(inv.XData);
%        lon1 = min(inv.YData);
%        axis([lat1-1 lat2+1 lon1-1 lon2+1]);
        
        % loop of animation
        for k = 1:length(t)
            axes(haxes)
            % checks if stop was pushed
            if stop_sign ==1
                stop_sign=0;
                break
            end
            
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
            delete(hanim_gg)
            hanim_gg = plot3m(latg(start:stop),long(start:stop),...
                    0,'LineWidth',3,'Color','r');
                % animation/stop animation should be one button
                % default values in gui
            set(hanim_gg,'Parent',t1);
            % plots animated groundtrack on flat world map
            delete(hanim_g)
            axes(ha)
            hanim_g = plot3m(latg(start:stop),long(start:stop),'ro',...
                'MarkerFaceColor','r');
            
            % plot circular swath on flat map
            delete(circ_swath_flat);
            
            [ circle_lat,circle_lon,r,locx ] = ...
                swath( alt(stop),FOV, latg(stop), long(stop));
            
            circ_swath_flat = plot3m(circle_lat,circle_lon,'y.');
            
            axes(haxes)
            % plots animated orbit
            delete(hanim_orb)
            hanim_orb = plot3m(lat(start:stop), lon(start:stop),...
                alt(start:stop),'Color','m');
            % plot satellite
            delete(hsat);
            hsat = plot3m(lat(stop), lon(stop), alt(stop),...
                'Marker','o','MarkerFaceColor','k');
            
            % plot circular swath on globe
            delete(circ_swath_globe);
            circ_swath_globe = plotm(circle_lat,circle_lon,'y.');
            set(circ_swath_globe,'Parent',t1);
            drawnow
            
            
             % flat map of swath
            axes(hswath)
            hswath.Visible = 'on' ;
            lat1 = min(circle_lat);
            lat2 = max(circle_lat);
            lon1 = min(circle_lon);
            lon2 = max(circle_lon);
            worldmap([lat1 lat2],[lon1 lon2]);
            load topo.mat topolegend topo;
%             geoshow(A,R)
            geoshow(topo,topolegend,'DisplayType','texturemap');
            if swath_running==1
                delete(circ_swath);
            end
            circ_swath = plot3m(circle_lat,circle_lon,'y','LineWidth',3);
            swath_text.Position = [thirdw+1.9*tenthw,tenthh*6.5,200,30];
            swath_text.Visible = 'on';
            
            quart = round(length(circle_lon)/4);
            half = round(length(circle_lon)/2);
            tf = quart+half;
            latscale = [circle_lat(quart+1) circle_lat(tf)];
            lonscale = [circle_lon(quart+1) circle_lon(tf)];
            diam = plot3m(latscale, lonscale,'y');
            str = strcat(num2str(round(abs(locx))),' km');
            diamtext = textm(latg(stop)+.5,long(stop),str,'Color','y','FontSize',15);
            
            swath_running =1;
            
            
        end
        
       elseif button_state == hObject.Min
         hanim.String = 'Animate';
         stop_sign=1;
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
   function [] = pushbutton(hObject, eventdata, handles)
    if strcmp(err_string,'')
       % deletes the handles of previous orbit data
       if been_here == 1
           %stop animation
            delete(data_g);
            delete(data_s);
            delete(data_g2);
            delete(data_swath1)
            delete(data_swath2)
            delete(cone_viz);
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
       
       endy = round(num_P*P/dt);
       
       % plot groundtrack on flat map
       data_g = scatterm(ha,latg(1:endy),long(1:endy),[],'r','o','filled');
       
       
       % plot groundtrack on 1st globe
       axes(haxes)
       data_g2 = plot3m(latg(1:endy),long(1:endy),zeros(1,endy),'LineWidth',3,...
            'Color','r');
        

       % plot satellite orbit on globe
       data_s = plot3m(lat(1:endy),lon(1:endy),alt(1:endy),...
           'LineWidth',2,'Color','m');
       
       
       % pass axis limits with orbit already drawn to animbutton
       %lim = axis;
        
       % calculate streets of coverage
       % num_sats per plane
       %num_sats=1;
       
       % plots yellow circle on flat world map to represent ground coverage
       axes(ha)
       [ circle_lat,circle_lon,r ] = swath( alt(1), FOV, latg(1), long(1));
       
       
       %outputs two lines, one line on either side of coverage area
       [ lats1, lons1] = coverage_line( latg(1:endy),long(1:endy),r );
       
       data_swath1 = plotm(lats1,lons1,0,'y.','MarkerSize',8);
       
       axes(haxes)
       data_swath2 = plot3m(lats1,lons1,0,'y.','MarkerSize',8);
       
       
       
      been_here=1;
       
        
       
    end  
   end

%function for when window is resized
function resizeui(hObject,callbackdata)
           
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
       %haxes.Position = [tenthw,tenthh/2+5,2*thirdw,halfh-tenthh];
       
       haxes.Position = [tenthw,tenthh/2+5,thirdw,halfh-tenthh];
       
       hswath.Position = [thirdw+2*tenthh,tenthh+5,thirdw-tenthh,halfh-2*tenthh];
end
    % load nasa map
%    nasa = wmsfind('nasa','SearchField','serverurl');
%    layer = nasa.refine('bluemarbleng','SearchField','layername','MatchType','exact');
%    [A,R] = wmsread(layer);

    load topo.mat topolegend topo;
    

  
   
   
     
   % flat world map
   ha = axes('Units','Pixels'); 
   worldmap world;
%    geoshow(A,R);
 geoshow(topo,topolegend,'DisplayType','texturemap');
   
   % globe on new axis
    haxes = axes('Units','Pixels'); 
    axis vis3d
    axesm ('globe');
    view(60,60)
    axis off
    
    % load surface of globe
    % loading high resolution blue marble image from nasa
 
%    hgeo=geoshow(A,R);
    hgeo = geoshow(topo,topolegend,'DisplayType','texturemap');
    demcmap(topo)
    land = shaperead('landareas','UseGeoCoords',true);
    hconts = plotm([land.Lat],[land.Lon],'Color','black');
   
   hswath = axes('Units','Pixels'); 
   hswath.Visible='off';
   
   
   %not sure
   %zoom on
   
   % Make the UI visible.
   f.Visible = 'on';
   
   % for next time smaller more detailed terrain map
   % lambda d h > compute spatial resolution num_pixel off nadir pushbroom
   % whiskbroom conical crossbreak
   % displays outputs pixel size vg  z swath area covered(acr) km^2/s

end