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
           hanim_gg hanim_g hanim_orb dt circ_swath_globe circ_swath_flat...
           circ_swath cone_handle cone_anim cover_anim f_orbit diam diamtext...
           D h lambda A R C A_orig vel ground_v int_t Nx Ny val del_x locx...
           hzoom sat_z cone_z;
       
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
   
   swath_running=0;
   
   Nxtext = uicontrol('Style','text','String','# Crosstrack Pixels, Nx');
   Nxedit = uicontrol('Style','edit','Callback',@nx_input);
   Nxedit.String=100;
   Nx=100;
  function [] = nx_input(H,E)
        nx = get(H,'string');
        Nx = str2double(nx);
  end

   Nytext = uicontrol('Style','text','String','# Alongtrack Pixels, Ny');
   Nyedit = uicontrol('Style','edit','Callback',@ny_input);
   Nyedit.String=100;
   Ny=100;
  function [] = ny_input(H,E)
        ny = get(H,'string');
        Ny = str2double(ny);
   end
   
   lambdatext = uicontrol('Style','text','String','Wavelength [nm]');
   lambdaedit = uicontrol('Style','edit','Callback',@lam_input);
   lambdaedit.String=500;
   lambda=500e-9;
   function [] = lam_input(H,E)
        lam = get(H,'string');
        lambda = str2double(lam)*1e-9;
   end
   
   Dtext = uicontrol('Style','text','String','Aperture Diameter [m]');
   Dedit = uicontrol('Style','edit','Callback',@D_input);
   Dedit.String=.1;
   D=.1;
  function [] = D_input(H,E)
        d = get(H,'string');
        D = str2double(d);
   end
   
   hvartext = uicontrol('Style','text','String','Altitude [m]');
   hvaredit = uicontrol('Style','edit','Callback',@hvar_input);
   h=7e6-Re;
   a=7e6;
   hvaredit.String=num2str(h,'%10.3e');
   function [] = hvar_input(H,E)
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
    function sensor_input(source,callbackdata)
        val = source.Value;
        maps = source.String; 
        


    end
   

    % button to control the field of view angle. 
   FOVtext = uicontrol('Style','text','String','Field of View [deg]');
   FOVedit = uicontrol('Style','edit','Callback',@FOV_input);
   FOVedit.String=30;
   FOV=30;
   function [] = FOV_input(H,E)
        FOV = get(H,'string');
        FOV = str2double(FOV);
        % field of view is limited by geometry of earth
        max_angle = 2*tan(Re/a)*180/pi; 
        if FOV > max_angle
            FOV = max_angle;
            FOVedit.String = num2str(FOV);
        end
   end

    %%%%% OUTPUTS %%%%%%
   outputtext = uicontrol('Style','text','String','OUTPUTS');
   del_xtext = uicontrol('Style','text','String','Spatial Resolution');
   del_xoutput = uicontrol('Style','text','String','');
   ground_vtext = uicontrol('Style','text','String','Ground Velocity');
   ground_voutput = uicontrol('Style','text','String','');
   int_ttext = uicontrol('Style','text','String','Integration Time');
   int_toutput = uicontrol('Style','text','String','');

    
   % text field empty until error
   herror = uicontrol('Style','text','String',err_string);
    
   
   
   % push button for animation of orbit and groundtrack
   hanim = uicontrol('Style','togglebutton','String',...
       'Animate','Callback',@animbutton);
  
   
   figure(f);
   % push button for animation of orbit and groundtrack
   hedit = uicontrol('Style','pushbutton','String',...
       'Edit Orbit Parameters','Callback',@editbutton);
   
   %%%% set positions of buttons
           
   % Get figure width and height
   hedit.Position = [figwidth-175 figheight-50 150 20];
   Nxtext.Position = [figwidth-175 figheight - 100 150 20];
   Nxedit.Position = [figwidth-150 figheight - 125 100 20];
   Nytext.Position = [figwidth-175 figheight - 150 150 20];
   Nyedit.Position = [figwidth-150 figheight - 175 100 20];
   lambdatext.Position = [figwidth-150 figheight - 200 100 20];
   lambdaedit.Position = [figwidth-150 figheight - 225 100 20];
   Dtext.Position = [figwidth-150 figheight - 250 100 20];
   Dedit.Position = [figwidth-150 figheight - 275 100 20];
   hvartext.Position = [figwidth-150 figheight - 300 100 20];
   hvaredit.Position = [figwidth-150 figheight - 325 100 20];
   FOVtext.Position = [figwidth-150 figheight - 350 100 20];
   FOVedit.Position = [figwidth-150 figheight - 375 100 20];
   sensortext.Position = [figwidth-150 figheight - 450 100 20];
   sensor.Position = [figwidth-175 figheight - 475 150 20];
   hanim.Position = [figwidth-175 figheight-550 150 20];
   herror.Position = [25 5 figwidth 20];


   halfh = round(figheight/2);
   thirdw = round(figwidth/3);
   tenthw = round(figwidth/15);
   tenthh = round(figheight/15);

   % outputs
   outputtext.Position = [tenthw figheight-500 150 20];
   del_xtext.Position = [tenthw figheight-525 150 20];
   del_xoutput.Position = [tenthw+125 figheight-525 100 20];
   ground_vtext.Position = [tenthw figheight-550 150 20];
   ground_voutput.Position = [tenthw+125 figheight-550 100 20];
   int_ttext.Position = [tenthw figheight-575 150 20];
   int_toutput.Position = [tenthw+125 figheight-575 100 20];
       
    
       
   
   a=7e6;
   i=10;
   e=0;
   RAAN=0;
   v=0;
   w=0;
   num_P=1;
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NEW FIGURE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   function [] = editbutton(hObject, eventdata, handles)
       % makes new figure for orbit parameters
        f_orbit = figure('Position',[200,200,300,450],...
       'Name','Orbit Parameters','NumberTitle','off');
   
       % semimajor axis edit text field with string
       atext = uicontrol('Style','text','String','Semimajor axis [m]');
       ah = uicontrol('Style','edit','Callback',@a_input);
       ah.String = 7e6;
       a=7e6;
       function [] = a_input(H,E)
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
   function [] = animbutton(hObject, eventdata, handles)
        button_state = hObject.Value;
        figure(f);
    if button_state == hObject.Max
        hanim.String = 'Stop';
       
        
       if done_that==1
           delete(hsat);
           delete(hanim_orb);
%            delete(circ_swath_flat);
           delete(circ_swath);
           delete(cone_anim);
%            delete(cover_anim);
            A = A_orig;
       end
       
       [ del_x ] = spatial_resolution(h,lambda,D );
       del_xoutput.String=strcat(num2str(del_x),' m');
       
       if val == 1
            %conical swath
            ground_voutput.String = '';
            int_toutput.String = '';
       elseif val ==2
            %whiskbroom
            [ int_t, ground_v ] = whiskbroom( h,vel(1), del_x,Nx,Ny,locx);
            ground_voutput.String = strcat(num2str(round(ground_v)),' m/s');
            int_toutput.String = strcat(num2str(int_t*1000),' ms');
        else
            %pushbroom
            [ int_t, ground_v ] = pushbroom( h,vel(1),del_x);
            ground_voutput.String = strcat(num2str(round(ground_v)),' m/s');
            int_toutput.String = strcat(num2str(int_t*1000),' ms');
       end

        
        % calulates altitude, longitude and latitude of satellite orbit
       [lat, lon, alt, P,vel] = Sat_orbit( a,e,i,RAAN,v,w,10*num_P,0,dt);
       alt=alt*Re;
       % calculate longitude, latitude of groundtrack. alt and P not used
       [latg, long] = Sat_orbit( a,e,i,RAAN,v,w,10*num_P,1,dt);
       
       [ circle_lat,circle_lon,r ] = swath( alt(1)/Re, FOV, latg(1), long(1));
        
        
        % here i set up handles so i can delete them in the loop.
        
        axes(haxes)
        % set up satellite
        hsat = plot3m(lat(1), lon(1), alt(1),'Marker','s',...
          'MarkerFaceColor','k');
        % plots animated orbit
        hanim_orb = plot3m(lat(1), lon(1),alt(1),'Color','m');
        %plot circular swath on flat map
%         circ_swath_flat = scatterm(ha,circle_lat,circle_lon,[],'y','.');


        % calculates and runs conical swath animation
        [xECEF, yECEF, zECEF] = conical_animation(r,alt(1),latg(1),long(1));
            
        cone_anim = surf(xECEF,yECEF,zECEF,'Facecolor','m','Facealpha',.4,'LineStyle','none');
 

       
        % set up for animation loop
        
        % t counts minutes of a day in seconds
        t=0:dt:round(10*num_P*P);
        tail_length = round(num_P*P/dt);
        
        % rotation of earth in deg/sec
        w_earth = 360/(24*3600);
        
        % setting up rotation transform
        t1 = hgtransform('Parent',haxes);
        set(hgeo,'Parent',t1);
%         set(hconts,'Parent',t1);
%         set(hanim_gg,'Parent',t1);
        set(cone_anim,'Parent',t1);
%         set(cover_anim,'Parent',t1);

        % plot satellite orbit on globe invisibly so axis doesnt shift
        % during animation
        endy = round(P/dt);
        plot3m(lat(1:endy),lon(1:endy),alt(1:endy),...
           'LineStyle','none','Marker','none');
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

            [ circle_lat,circle_lon,r,locx] = swath( ...
                alt(stop)/Re,FOV, latg(stop), long(stop));

            
            

            
            % plot cone under sat 
           [xECEF, yECEF, zECEF] = conical_animation(r,alt(stop),latg(stop),long(stop));
            
           cone_anim.XData = xECEF;
           cone_anim.YData = yECEF;
           cone_anim.ZData = zECEF;
           

            
            [rmin,rmax,cmin,cmax] = coverage_lineECEF( circle_lat,circle_lon,r,A,R );
            
             if abs(cmin-cmax) > 512/2;
                A(rmin:rmax,cmax:512,:) = C(rmin:rmax,cmax:512,:);
                minny = cmin;
                cmin=1;
                cmax = minny;
            end
            
            A(rmin:rmax,cmin:cmax,:) = C(rmin:rmax,cmin:cmax,:);
            
            
            delete(hgeo)
            hgeo = geoshow(A,R);
%             demcmap(A)
            set(hgeo,'Parent',t1);
            
            % plots animated orbit
            delete(hanim_orb)
            hanim_orb = plot3m(lat(start:stop), lon(start:stop),...
                alt(start:stop),'Color','m','LineWidth',3);
            % plot satellite
            delete(hsat);
            hsat = plot3m(lat(stop), lon(stop), alt(stop),...
                'Marker','o','MarkerFaceColor','k');
            
           
            
            %%%%%%% flat map of swath
            axes(hswath)
            hswath.Visible = 'on' ;
            
            % converts lat and lon of swath to rows and cols of pixels
            circ_rows = zeros(size(circle_lat));
            circ_cols = zeros(size(circle_lat));
            for p=1:length(circle_lat)
                [r1,c1] = setpostn(A,R,circle_lat(p),circle_lon(p));
                circ_rows(p)=r1;
                circ_cols(p)=c1;
            end
            
           
            % creates image to be displayed by swath
            cmap = hgeo.CData;
            
            [B] = flat_swath( circ_rows,circ_cols,latg(stop),r,256,512,A );
             
            
            %displays image
            swat = geoshow(B,R);
            colormap(winter)
            
            % changing tick and axis values on swath
            
%             image(A(minrow:maxrow,mincol:maxcol),'CDataMapping','scaled')
%             colormap(A.CData)
            hold on
            xlimits = xlim;
            x1 = xlimits(1);
            x2 = x1+(xlimits(2)-xlimits(1))/2;
            x3 = xlimits(2);
            ylimits = ylim;
            y1 = ylimits(1);
            y2 = y1+(ylimits(2)-ylimits(1))/2;
            y3 = ylimits(2);
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
            hswath.YTickLabel = {lat3text,lattext,lat1text};        
            hswath.XTickLabel = {lon1text,lontext,lon3text}; 
            
            
            if swath_running==1
                delete(circ_swath);
                delete(diam);
                delete(diamtext);
                delete(sat_z);
                delete(cone_z);
%                 delete(test);
%                 delete(test2);
%                 delete(test3);
            end
            
            %plot circular swath on flat swath
            ang=0:0.1:2*pi+.1; 
            xp=(x2-x1)*cos(ang);
            yp=(y2-y1)*sin(ang);
            circ_swath = plot(x2+xp,y2+yp,'y','LineWidth',4);
            % plot title
            swath_text.Position = [thirdw+1.9*tenthw,tenthh*7,200,30];
            swath_text.Visible = 'on';
            % plot diameter and diameter text
            xdiam = x1:1:x3;
            ydiam = ones(size(xdiam))*y2;
            diam = plot(xdiam,ydiam,'y','LineWidth',3);
            str = strcat(num2str(round(abs(locx))),' km');
            diamtext = text(x2+x2/4,y2+y2/4,str,'Color','y','FontSize',15);
            hold off
            
            
            
            %%%% set up hzoom axis
            axes(hzoom)
            hzoom.Visible = 'on' ;
            hold on
            %creates image to be displayed by swath
            [ circle_lat2,circle_lon2,r2,locx2] = swath( ...
                alt(stop)/Re,FOV*2, latg(stop), long(stop));
            circ_rows2 = zeros(size(circle_lat2));
            circ_cols2 = zeros(size(circle_lat2));
            for p=1:length(circle_lat2)
                [r1,c1] = setpostn(A,R,circle_lat2(p),circle_lon2(p));
                circ_rows2(p)=r1;
                circ_cols2(p)=c1;
            end
            [Q] = flat_swath( circ_rows2,circ_cols2,latg(stop),r2,256,512,A);
            
            
            
            xlimits2 = xlim;
            x4 = xlimits2(1);
            x6 = xlimits2(2);
            x5 = x4+(x6-x4)/2;


            ylimits2 = ylim;
            y4 = ylimits2(1);
            y6 = ylimits2(2);
            y5 = y4+(y6-y4)/2;
            
           
            
            geoshow(Q,R);
            view(60,40)
            axis off
            
            % approximate altitude in pixels
            % used 512 pixels (longitude) = 4e7m cirfumference at equator of earth;
            zeta = alt(stop)/4e7*512;
            sat_z = plot3(x5,y5,zeta,'ks','MarkerFaceColor','k','MarkerSize',15);
%             test = plot3(x4,y4,0,'rs','MarkerFaceColor','r','MarkerSize',10);
%             test2 = plot3(x5,y5,0,'cs','MarkerFaceColor','g','MarkerSize',10);
%             test3 =plot3(x6,y6,0,'gs','MarkerFaceColor','g','MarkerSize',10);
           
            
            
            
            % plot cone on swath
            xcirc2 = xp+x5-3;
            ycirc2 = yp+y5+2;
            xsurf = zeros(2,length(xcirc2));
            xsurf(1,:) = xcirc2; 
            xsurf(2,:) = x5;

            ysurf = zeros(2,length(ycirc2));
            ysurf(1,:) = ycirc2;
            ysurf(2,:) = y5;

            zsurf = zeros(2,length(xcirc2));
            zsurf(1,:) = 0;
            zsurf(2,:) = zeta;
            
            cone_z = surf(xsurf,ysurf,zsurf,'Facecolor','m','Facealpha',.4,'LineStyle','none');
            
            swath_running=1;
            drawnow
        end
        
       elseif button_state == hObject.Min
         hanim.String = 'Animate';
         stop_sign=1;
    end
        done_that = 1;
   end



%function for when window is resized
function resizeui(hObject,callbackdata)
    
        figure(f)
           
       % Get figure width and height
       figwidth = f.Position(3);
       figheight = f.Position(4);
       
       hedit.Position = [figwidth-175 figheight-50 150 20];
       Nxtext.Position = [figwidth-175 figheight - 100 150 20];
       Nxedit.Position = [figwidth-150 figheight - 125 100 20];
       Nytext.Position = [figwidth-175 figheight - 150 150 20];
       Nyedit.Position = [figwidth-150 figheight - 175 100 20];
       lambdatext.Position = [figwidth-150 figheight - 200 100 20];
       lambdaedit.Position = [figwidth-150 figheight - 225 100 20];
       Dtext.Position = [figwidth-150 figheight - 250 100 20];
       Dedit.Position = [figwidth-150 figheight - 275 100 20];
       hvartext.Position = [figwidth-150 figheight - 300 100 20];
       hvaredit.Position = [figwidth-150 figheight - 325 100 20];
       
       
       FOVtext.Position = [figwidth-150 figheight - 350 100 20];
       FOVedit.Position = [figwidth-150 figheight - 375 100 20];
       
       
       
       sensortext.Position = [figwidth-150 figheight - 450 100 20];
       sensor.Position = [figwidth-175 figheight - 475 150 20];
       
       
       
       hanim.Position = [figwidth-175 figheight-550 150 20];
       herror.Position = [25 5 figwidth 20];
       
       
       halfh = round(figheight/2);
       thirdw = round(figwidth/3);
       tenthw = round(figwidth/15);
       tenthh = round(figheight/15);
       
       % outputs
       outputtext.Position = [tenthw figheight-500 150 20];
       del_xtext.Position = [tenthw figheight-525 150 20];
       del_xoutput.Position = [tenthw+125 figheight-525 100 20];
       ground_vtext.Position = [tenthw figheight-550 150 20];
       ground_voutput.Position = [tenthw+125 figheight-550 100 20];
       int_ttext.Position = [tenthw figheight-575 150 20];
       int_toutput.Position = [tenthw+125 figheight-575 100 20];
       
    
       % maps
       haxes.Position = [thirdw+2*tenthw,halfh,thirdw-tenthw,halfh-tenthh];
       hswath.Position = [tenthw,halfh,thirdw-tenthw,halfh-tenthh];
end




   nasa = wmsfind('nasa','SearchField','serverurl');
   layer = nasa.refine('bluemarbleng','SearchField','layername','MatchType','exact');
   [A,R] = wmsread(layer);
   A_orig = A;
   C = (A+30)*2;
  
     
   % flat world map
%    ha = axes('Units','Pixels'); 
%    worldmap world;
%    load topo.mat topolegend topo;
%    geoshow(topo,topolegend,'DisplayType','texturemap')
   
   % globe on new axis
    earth = referenceEllipsoid('earth','m');
    haxes = axes('Units','Pixels'); 
    axis vis3d
    axesm ('globe','Grid', 'on','Geoid',earth);
    view(60,60)
    axis off
    
    hgeo = geoshow(A,R);
    
    % load surface of globe
%     load topo.mat topolegend topo;
%     hgeo = geoshow(topo,topolegend,'DisplayType','texturemap');
%     demcmap(topo)
%     land = shaperead('landareas','UseGeoCoords',true);
%     hconts = plotm([land.Lat],[land.Lon],'Color','black');
    
   
   hswath = axes('Units','Pixels'); 
   hswath.Visible='off';
   
   hzoom =  axes('Units','Pixels'); 
   axis off
   
   haxes.Position = [thirdw+2*tenthw,halfh,thirdw-tenthw,halfh-tenthh];
   hswath.Position = [tenthw,halfh,thirdw-tenthw,halfh-tenthh];
   hzoom.Position = [thirdw+2*tenthw,tenthh,thirdw-tenthw,halfh-2*tenthh];
   
   % Make the UI visible.
   f.Visible = 'on';
   
end