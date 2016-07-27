function Swath_GUI
% Allegra Moran
% 6/6/16
% This GUI is for showing the orbit and swath of a group of satellites 

   %  Create and then hide the UI as it is being constructed.
   close all;
   figwidth = 1000;
   figheight = 600;
   f = figure('Visible','off','Position',[360,500,figwidth,figheight],...
       'Name','Remote Sensing Design','NumberTitle','off',...
       'SizeChangedFcn',@resizeui);
   
   % global variables passed between uicontrols and functions
   global a e i RAAN w v num_P lat lon alt latg long P hsat...
           hanim_orb dt sat_z cone_z...
           circ_swath cone_anim f_orbit diam diamtext...
           D h lambda A R C A_orig vel ground_v int_t Nx Ny val del_x locx...
           hzoom ;
       
       % global handles of button
   global ah atext etext eh itext inch RAANtext RAANh vtext vh wtext...
       wh numPtext numPh slidertext slideredit hslider Nytext Nyedit...
       Nxedit Nxtext ground_voutput ground_vtext int_toutput int_ttext;
   
   halfh = round(figheight/2);
   thirdw = round(figwidth/3);
   tenthw = round(figwidth/15);
   tenthh = round(figheight/15);
   
   
   %radius of earth
   Re = 6.371e6;
   
   %  Construct the components.
   err_string='';
   
   stop_sign=0;
   
   
   
   Nxtext = uicontrol('Style','text','String','# Crosstrack Pixels, Nx');
   Nxedit = uicontrol('Style','edit','Callback',@nx_input);
   Nxedit.String=100;
   Nx=100;
  function [] = nx_input(H,~)
        nx = get(H,'string');
        Nx = str2double(nx);
  end

   Nytext = uicontrol('Style','text','String','# Alongtrack Pixels, Ny');
   Nyedit = uicontrol('Style','edit','Callback',@ny_input);
   Nyedit.String=100;
   Ny=100;
  function [] = ny_input(H,~)
        ny = get(H,'string');
        Ny = str2double(ny);
   end
   
   lambdatext = uicontrol('Style','text','String','Wavelength [nm]');
   lambdaedit = uicontrol('Style','edit','Callback',@lam_input);
   lambdaedit.String=500;
   lambda=500e-9;
   function [] = lam_input(H,~)
        lam = get(H,'string');
        lambda = str2double(lam)*1e-9;
   end
   
   Dtext = uicontrol('Style','text','String','Aperture Diameter [m]');
   Dedit = uicontrol('Style','edit','Callback',@D_input);
   Dedit.String=.1;
   D=.1;
  function [] = D_input(H,~)
        d = get(H,'string');
        D = str2double(d);
   end
   
   hvartext = uicontrol('Style','text','String','Altitude [m]');
   hvaredit = uicontrol('Style','edit','Callback',@hvar_input);
   h=7e6-Re;
   a=7e6;
   hvaredit.String=num2str(h,'%10.3e');
   function [] = hvar_input(H,~)
        h = get(H,'string');
        h = str2double(h);
        a = h+Re;
%         ah.String =num2str(a,'%10.3e'); 
        
   end
        
    % Create sensor pop up menu
    sensortext = uicontrol('Style','text','String','Select Sensor Option');
    sensor = uicontrol('Style','popupmenu','String',...
        {'Conical','Whisk broom','Push broom'},'Callback',@sensor_input);
    val=1;
    function sensor_input(source,~)
        val = source.Value;
        
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
            FOV = max_angle;
            FOVedit.String = num2str(FOV);
        end
   end

   Nbandtext = uicontrol('Style','text','String','Number of bands, Nband');
   Nbandedit = uicontrol('Style','edit','Callback',@Nband_input);
   Nbandedit.String=30;
   Nband=30;
   function [] = Nband_input(H,~)
        Nbandz = get(H,'string');
        Nband = str2double(Nbandz);
   end

   nadirtext = uicontrol('Style','text','String','Off Nadir Angle [deg]');
   nadiredit = uicontrol('Style','edit','Callback',@nadir_input);
   nadiredit.String=10;
   nadir=10;
   function [] = nadir_input(H,~)
        nad = get(H,'string');
        nadir = str2double(nad);
   end

   
   

   bitstext = uicontrol('Style','text','String','Bits / pixel');
   bitsedit = uicontrol('Style','edit','Callback',@bits_input);
   bitsedit.String=12;
   bits=12;
   function [] = bits_input(H,~)
        bit = get(H,'string');
        bits = str2double(bit);
   end

   tracktext = uicontrol('Style','text','String','Track Angle [deg]');
   trackedit = uicontrol('Style','edit','Callback',@track_input);
   trackedit.String=100;
   track=100;
   function [] = track_input(H,~)
        tra = get(H,'string');
        track = str2double(tra);
   end

   alongtext = uicontrol('Style','text','String','Along Track Angle [deg]');
   alongedit = uicontrol('Style','edit','Callback',@along_input);
   alongedit.String=10;
   along=10;
   function [] = along_input(H,~)
        alo = get(H,'string');
        along = str2double(alo);
   end
   alongtext.Visible='off';
   alongedit.Visible='off';

   crosstext = uicontrol('Style','text','String','Cross Track Angle [deg]');
   crossedit = uicontrol('Style','edit','Callback',@cross_input);
   crossedit.String=10;
   cross=10;
   function [] = cross_input(H,~)
        cro = get(H,'string');
        cross = str2double(cro);
   end
   crosstext.Visible='off';
   crossedit.Visible='off';
   

    %%%%% OUTPUTS %%%%%%
   outputtext = uicontrol('Style','text','String','OUTPUTS');
   del_xtext = uicontrol('Style','text','String','Spatial Resolution');
   del_xoutput = uicontrol('Style','text','String','');
   ground_vtext = uicontrol('Style','text','String','Ground Velocity');
   ground_voutput = uicontrol('Style','text','String','');
   int_ttext = uicontrol('Style','text','String','Integration Time');
   int_toutput = uicontrol('Style','text','String','');
   swath_area_text = uicontrol('Style','text','string','Swath Area');
   swath_area_output = uicontrol('Style','text','string','');
   data_rate_text = uicontrol('Style','text','string','Data Rate');
   data_rate_output = uicontrol('Style','text','string','');
  
    
   % text field empty until error
   herror = uicontrol('Style','text','String',err_string);
    
   
   
   % push button for animation of orbit and groundtrack
   hanim = uicontrol('Style','togglebutton','String',...
       'Animate','Callback',@animbutton);
  
 
   % push button for animation of orbit and groundtrack
   hedit = uicontrol('Style','pushbutton','String',...
       'Edit Orbit Parameters','Callback',@editbutton);
   

   a=7e6;
   i=10;
   e=0;
   RAAN=0;
   v=0;
   w=0;
   num_P=1;
   
   
   %%%%%%%%%%%%%%%%% NEW FIGURE edit orbit parameters %%%%%%%%%%%%%%%%%%%%%
   function [] = editbutton(~,~,~)
       % makes new figure for orbit parameters
        f_orbit = figure('Position',[200,200,300,450],...
       'Name','Orbit Parameters','NumberTitle','off');
   
       % semimajor axis edit text field with string
       atext = uicontrol('Style','text','String','Semimajor axis [m]');
       ah = uicontrol('Style','edit','Callback',@a_input);
       ah.String = 7e6;
       a=7e6;
       function [] = a_input(H,~)
            aa = get(H,'string');
            a = str2double(aa);
            h=a-Re;
%             hvaredit.String = num2str(h,'%10.3e');
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
       function [] = e_input(H,~)
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
       function [] = inc_input(H,~)
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
       function [] = RAAN_input(H,~)
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
       function [] = v_input(H,~)
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
       function [] = w_input(H,~)
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
       function [] = numP_input(H,~)
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
       function [] = slidertext_input(H,~)
            dt = get(H,'string');
            dt = str2double(dt);
       end
       function [] = slider_input(H,~)
            dt = 60+H.Value*600;
            slideredit.String = dt;
       end
       
       %set positions for all uicontrols
       spacing = 25;
       y = 20;
       x = 100;
       
       atext.Position =      [x y+16*spacing 100 20];
       ah.Position =         [x y+15*spacing 100 20];
       etext.Position =      [x y+14*spacing 100 20];
       eh.Position =         [x y+13*spacing 100 20];
       itext.Position =      [x y+12*spacing 100 20];
       inch.Position =       [x y+11*spacing 100 20];
       RAANtext.Position =   [x-50 y+10*spacing 200 20];
       RAANh.Position =      [x y+9*spacing 100 20];
       vtext.Position =      [x y+8*spacing 100 20];
       vh.Position =         [x y+7*spacing 100 20];
       wtext.Position =      [x-50 y+6*spacing 200 20];
       wh.Position =         [x y+5*spacing 100 20];
       numPtext.Position =   [x y+4*spacing 100 20];
       numPh.Position =      [x y+3*spacing 100 20];
       slidertext.Position = [x y+2*spacing 100 20];
       slideredit.Position = [x y+1*spacing 100 20];
       hslider.Position =    [x y 100 20];
        
       figure(f);
   end


   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ANIMATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   done_that = 0;
   function [] = animbutton(hObject, ~, ~)
        button_state = hObject.Value;
        figure(f);
    if button_state == hObject.Max
        hanim.String = 'Stop';
       
        
       if done_that==1
           delete(hsat);
           delete(hanim_orb);
           delete(circ_swath);
           delete(cone_anim);
           delete(sat_z);
           delete(cone_z);
           delete(diam);
           delete(diamtext);
            A = A_orig;
       end
       
       [ del_x ] = spatial_resolution(h,lambda,D );
       del_xoutput.String=strcat(num2str(del_x),' m');
       
       if val == 1
            %conical swath
            ground_voutput.String = '';
            int_toutput.String = '';
            area = round(pi*(locx/2)^2);
            swath_area_output.String = strcat(num2str(area,'%10.2e'),' km^2');
            nadirtext.Visible = 'on';
            nadiredit.Visible = 'on';
            tracktext.Visible='on';
            trackedit.Visible='on';
            alongtext.Visible = 'off';
            alongedit.Visible = 'off';
            crosstext.Visible='off';
            crossedit.Visible='off';
           
       elseif val ==2
            %whiskbroom
            [ int_t, ground_v,data_rate ] = whiskbroom( h,vel(1),del_x,Nx,Ny,locx,bits,Nband);
            ground_voutput.String = strcat(num2str(round(ground_v)),' m/s');
            int_toutput.String = strcat(num2str(int_t*1000),' ms');
            data_rate_output.String = strcat(num2str(data_rate/1e6,'%10.2e'),' Mbps');
            nadirtext.Visible = 'off';
            nadiredit.Visible = 'off';
            tracktext.Visible='off';
            trackedit.Visible='off';
            alongtext.Visible = 'on';
            alongedit.Visible = 'on';
            crosstext.Visible='on';
            crossedit.Visible='on';
        else
            %pushbroom
            [ int_t, ground_v,data_rate ] = pushbroom( h,vel(1),del_x,Nx,Nband,bits);
            ground_voutput.String = strcat(num2str(round(ground_v)),' m/s');
            int_toutput.String = strcat(num2str(int_t*1000),' ms');
            data_rate_output.String = strcat(num2str(data_rate/1e6,'%10.2e'),' Mbps');
            nadirtext.Visible = 'off';
            nadiredit.Visible = 'off';
            tracktext.Visible='off';
            trackedit.Visible='off';
            alongtext.Visible = 'on';
            alongedit.Visible = 'on';
            crosstext.Visible='on';
            crossedit.Visible='on';
       end

        
        % calulates altitude, longitude and latitude of satellite orbit
       [lat, lon, alt, P,vel] = Sat_orbit( a,e,i,RAAN,v,w,10*num_P,0,dt);
       alt=alt*Re;
       % calculate longitude, latitude of groundtrack. alt and P not used
       [latg, long] = Sat_orbit( a,e,i,RAAN,v,w,10*num_P,1,dt);
       
       [ ~,~,r ] = swath( alt(1)/Re, FOV, latg(1), long(1));
        
        
        % here i set up handles so i can delete them in the loop.
        
        axes(haxes)
        % set up satellite
        hsat = plot3m(lat(1), lon(1), alt(1),'Marker','o',...
          'MarkerFaceColor','k','MarkerSize',8);
        % plots animated orbit
        hanim_orb = plot3m(lat(1), lon(1),alt(1),'Color','m','LineWidth',3);
        
        globestruct = gcm(haxes);


        % calculates and runs conical swath animation on globe
        [xECEF, yECEF, zECEF] = conical_animation(r,alt(1),latg(1),long(1));
            
        cone_anim = surf(xECEF,yECEF,zECEF,'Facecolor','m','Facealpha',.4,'LineStyle','none');
 

       
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
        set(cone_anim,'Parent',t1);

        % plot satellite orbit on globe invisibly so axis doesnt shift
        % during animation
        endy = round(P/dt);
        plot3m(lat(1:endy),lon(1:endy),alt(1:endy),...
           'LineStyle','none','Marker','none');
       
       % determine heading
       heading = zeros(length(latg)-1);
       for x=1:length(latg)-1
            slat = [latg(x) latg(x+1)];
            slon = [long(x) long(x+1)];
            [heading(x), ~] = legs(slat,slon,'gc');
       end
       
       % setting up swath
        axes(hswath)
        axis tight
        ang=0:0.1:2*pi+.1;
        circ_swath = plot(0,0,'m','Visible','off','LineWidth',4);
        diam = plot(0,0,'m','Visible','off','LineWidth',3);
        diamtext = text(0,0,'','Color','m','Visible','off','FontSize',15,...
        'FontWeight','bold');
    
        % setting up zoom axis
        axes(hzoom)
        axis tight
        sat_z = plot3m(0,0,0,'ks','Markerfacecolor','k','MarkerSize',15,'Visible','off');
        rix = [0 0;0 0];
        cone_z = surf(rix,rix,rix,'Visible','off','Facecolor','m',...
            'Facealpha',.4,'LineStyle','none');
    
        zlim(hzoom, [0 .1]);
    
        % loop of animation
        for k = 1:length(t)
%             axes(haxes)
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
            
            %%%%%%%% globe
            
            % convert lat and lon to x,y,z units to update globe objects
            [ox, oy, oz] = mfwdtran(globestruct,lat(start:stop), ...
                lon(start:stop),alt(start:stop));
            
            % plot satellite on globe
            hsat.XData = ox(end);
            hsat.YData = oy(end);
            hsat.ZData = oz(end);
            
            %plot animated orbit on globe
            hanim_orb.XData = ox;
            hanim_orb.YData = oy;
            hanim_orb.ZData = oz;

            % get circle lat and lon
            [ circle_lat,circle_lon,r,locx] = swath( ...
                alt(stop)/Re,FOV, latg(stop), long(stop));

            % plot cone under sat on globe
           [xECEF, yECEF, zECEF] = conical_animation(r,alt(stop),latg(stop),long(stop));
           cone_anim.XData = xECEF;
           cone_anim.YData = yECEF;
           cone_anim.ZData = zECEF;
           

            
            [rmin,rmax,cmin,cmax] = coverage_lineECEF( circle_lat,circle_lon,r,A,R);
             if abs(cmin-cmax) > 512/2;
                % if at  edge of image
                if latg(stop)+r >= 90
                    %if at north pole
                    disp('north')
                    hgeo.CData(1:rmax,cmin:cmax,:) = C(1:rmax,cmin:cmax,:);
                    hgeo.CData(1:rmin,cmin:cmax,:) = C(1:rmin,cmin:cmax,:);
%                     A(1:rmin,:,:) = C(1:rmin,:,:);
                elseif latg(stop)-r <= -90
                        disp('south')
                else
                    % else crossing prime meridian
                    hgeo.CData(rmin:rmax,cmax:512,:) = C(rmin:rmax,cmax:512,:);
                    hgeo.CData(rmin:rmax,1:cmin,:) = C(rmin:rmax,1:cmin,:);
                    disp('true')
                    
                end
            
             else
                hgeo.CData(rmin:rmax,cmin:cmax,:) = C(rmin:rmax,cmin:cmax,:);
             end
            
            
            
           
            
            %%%%%%% flat map of swath %%%%%%%
            
            origin = [round(latg(stop)) round(long(stop) -heading(stop))];
            setm(hswath,'Origin',origin,'FLatLimit',[-Inf r])
            
            % if latitude starts decreasing
            % compare lats 2 pts apart. pts next door will give
            % errors around poles occasionally because its not accurate to
            % tell if increasing/decreasing. Needs more testing to be sure
            if stop>2 && latg(stop-1)>latg(stop+1)
                view(hswath, [180,90])
            else
                view(hswath,[0,90])
            end
            
            
            % getting the limits of the swath
            xlimits = xlim(hswath);
            x1 = xlimits(1);
            x2 = x1+(xlimits(2)-xlimits(1))/2;
            x3 = xlimits(2);
            ylimits = ylim(hswath);
            y1 = ylimits(1);
            y2 = y1+(ylimits(2)-ylimits(1))/2;
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
            xp=(rad-rad/r)*cos(ang);
            yp=(rad-rad/r)*sin(ang);
            circ_swath.XData = x2+xp;
            circ_swath.YData = y2+yp;
            circ_swath.Visible = 'on';
            
            
            % plot diameter and diameter text
            xdiam = x1+1/r*rad:.001:x3-1/r*rad;
            ydiam = ones(size(xdiam))*y2;
            diam.XData = xdiam;
            diam.YData = ydiam;
            diam.Visible='on';
            
            % updating diameter km text
            str = strcat(num2str(round(abs(locx))),' km');
            diamtext.Visible = 'on';
            diamtext.String = str;
            diamtext.Position = [x2-x1/8 y2-y1/4];
            
            swathgeo.CData(rmin:rmax,cmin:cmax,:) = C(rmin:rmax,cmin:cmax,:);

            
            
            %%%% hzoom axis
            setm(hzoom,'Origin',origin,'FLatLimit',[-Inf r*3])
            
            hzoom.Visible = 'on' ;
            
           
            
           
            xlimits2 = xlim(hzoom);
            x4 = xlimits2(1);
            x6 = xlimits2(2);
            x5 = x4+(x6-x4)/2;


            ylimits2 = ylim(hzoom);
            y4 = ylimits2(1);
            y6 = ylimits2(2);
            y5 = y4+(y6-y4)/2;
            
            zoomgeo.CData = hgeo.CData;
            
            % approximate altitude
            zeta = alt(stop)/1e7;
            
            zlim(hzoom, [0 zeta]);
            
            % updating sat on zoomie axis
            sat_z.XData = x5;
            sat_z.YData = y5;
            sat_z.ZData = zeta;
            sat_z.Marker = 's';
            sat_z.Visible='on';
            
            
            
            % plot cone on zoom
            axis(hzoom,'off')
            xcone=(rad+rad/r)*cos(ang);
            ycone=(rad+rad/r)*sin(ang);
            xcirc2 = xcone+x5;
            ycirc2 = ycone+y5;
            xsurf = zeros(2,length(xcirc2));
            xsurf(1,:) = xcirc2; 
            xsurf(2,:) = x5;

            ysurf = zeros(2,length(ycirc2));
            ysurf(1,:) = ycirc2;
            ysurf(2,:) = y5;

            zsurf = zeros(2,length(xcirc2));
            zsurf(1,:) = 0;
            zsurf(2,:) = zeta;
            
            % updating cone surface
            cone_z.XData=xsurf;
            cone_z.YData=ysurf;
            cone_z.ZData=zsurf;

            cone_z.Visible = 'on';
            drawnow
        end
        
       elseif button_state == hObject.Min
         hanim.String = 'Animate';
         stop_sign=1;
    end
        done_that = 1;
   end



%function for when window is resized
function resizeui(~,~)
    
        figure(f)
           
       % Get figure width and height
       figwidth = f.Position(3);
       figheight = f.Position(4);
       
       %1st col on main figure
       hedit.Position = [figwidth-212 figheight-25 150 20];
       Nxtext.Position = [figwidth-275 figheight - 50 150 20];
       Nxedit.Position = [figwidth-250 figheight - 75 100 20];
       Nytext.Position = [figwidth-275 figheight - 100 150 20];
       Nyedit.Position = [figwidth-250 figheight - 125 100 20];
       lambdatext.Position = [figwidth-250 figheight - 150 100 20];
       lambdaedit.Position = [figwidth-250 figheight - 175 100 20];
       Dtext.Position = [figwidth-250 figheight - 200 100 20];
       Dedit.Position = [figwidth-250 figheight - 225 100 20];
       hvartext.Position = [figwidth-250 figheight - 250 100 20];
       hvaredit.Position = [figwidth-250 figheight - 275 100 20];
       Nbandtext.Position = [figwidth-275 figheight - 300 150 20];
       Nbandedit.Position = [figwidth-250 figheight - 325 100 20];
       bitstext.Position = [figwidth-250 figheight - 350 100 20];
       bitsedit.Position = [figwidth-250 figheight - 375 100 20];
       FOVtext.Position = [figwidth-250 figheight - 400 100 20];
       FOVedit.Position = [figwidth-250 figheight - 425 100 20];
       
       
       %2nd column
       sensortext.Position = [figwidth-125 figheight - 50 100 20];
       sensor.Position = [figwidth-130 figheight - 75 110 20];
       
       nadirtext.Position = [figwidth-125 figheight - 100 100 20];
       nadiredit.Position = [figwidth-125 figheight - 125 100 20];
       tracktext.Position = [figwidth-125 figheight - 150 100 20];
       trackedit.Position = [figwidth-125 figheight - 175 100 20];
       
       alongtext.Position = [figwidth-150 figheight - 200 150 20];
       alongedit.Position = [figwidth-125 figheight - 225 100 20];
       crosstext.Position = [figwidth-150 figheight - 250 150 20];
       crossedit.Position = [figwidth-125 figheight - 275 100 20];
       
       %animate button
       hanim.Position = [figwidth-212 figheight-575 150 20];
       herror.Position = [25 5 figwidth 20];
       
       
       halfh = round(figheight/2);
       thirdw = round(figwidth/3);
       tenthw = round(figwidth/15);
       tenthh = round(figheight/15);
       
       % outputs
       ih = figheight-375;
       step=25;
       outputtext.Position =        [tenthw ih-step 150 20];
       del_xtext.Position =         [tenthw ih-step*2 150 20];
       del_xoutput.Position =       [tenthw+125 ih-step*2 100 20];
       ground_vtext.Position =      [tenthw ih-step*3 150 20];
       ground_voutput.Position =    [tenthw+125 ih-step*3 100 20];
       int_ttext.Position =         [tenthw ih-step*4 150 20];
       int_toutput.Position =       [tenthw+125 ih-step*4 100 20];
       swath_area_text.Position =   [tenthw ih-step*5 150 20];
       swath_area_output.Position = [tenthw+125 ih-step*5 100 20];
       data_rate_text.Position =    [tenthw ih-step*6 150 20];
       data_rate_output.Position =  [tenthw+125 ih-step*6 100 20];
       
    
       % maps
       haxes.Position = [thirdw+tenthw,halfh,thirdw,halfh-tenthh];
       hswath.Position = [tenthw,halfh,thirdw,halfh-tenthh];
       hzoom.Position = [thirdw+tenthw,tenthh,thirdw,halfh-2*tenthh];
end




   nasa = wmsfind('nasa','SearchField','serverurl');
   layer = nasa.refine('bluemarbleng','SearchField','layername','MatchType','exact');
   [A,R] = wmsread(layer);
   A_orig = A;
   C = (A+30)*2;
  
   % globe on new axis
    earth = referenceEllipsoid('earth','m');
    haxes = axes('Units','Pixels'); 
    axis vis3d
    axesm ('globe','Grid', 'on','Geoid',earth);
    view(60,60)
    axis off
    hgeo = geoshow(A,R);
    

    
   
       % create axis 
    hswath = axes('Units','pixels');
    hswath = axesm('ortho','Origin',[0 0],'FLatLimit',[-Inf 0], ...
            'frame','off','grid','off');
    
    box(hswath,'off')
    set(hswath,'color','none');
    axis on
    swathgeo = geoshow(A,R);
   
    % zoom axis
   hzoom = axes('Units','Pixels'); 
   
   hzoom = axesm('ortho','Origin',[0 0],'FLatLimit',[-Inf 0], ...
            'frame','off','grid','off');
   view(60,40)
   zoomgeo = geoshow(A,R);
%    axis tight
   zoom(hzoom,2)
   box(hzoom,'off')
   axis off
   
   
   % Make the UI visible.
   f.Visible = 'on';
   
end