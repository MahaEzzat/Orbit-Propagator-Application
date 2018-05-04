function Orbit_propagator
%%figure design
clear all
scrsz = get(0,'ScreenSize');
f = figure('Visible','off','Position',[scrsz(1),scrsz(2),scrsz(3)-200,scrsz(4)],'Color','k');
% Construct the components.
hR_V = uicontrol('Style','pushbutton','String','R_V',...
'Position',[140,scrsz(4)-70,70,25],...
'Callback',{@R_Vbutton_Callback});
h6OEs = uicontrol('Style','pushbutton','String','6OEs',...
'Position',[15,scrsz(4)-70,70,25],...
'Callback',{@OEbutton_Callback});
hTLE = uicontrol('Style','pushbutton',...
'String','TLE',...
'Position',[300,scrsz(4)-70,70,25],...
'Callback',{@TLEbutton_Callback});

orbit_plot = uicontrol('Style','pushbutton',...
'String','plot',...
'Position',[430,scrsz(4)-70,70,25],...
'Callback',{@orbit_plotbutton_Callback});

num = uicontrol('Style','pushbutton',...
'String','Numerical','Visible','off',...
'Position',[140,scrsz(4)-300,70,25],...
'Callback',{@num_button_Callback});

keplar = uicontrol('Style','pushbutton',...
'String','Keplar','Visible','off',...
'Position',[50,scrsz(4)-300,70,25],...
'Callback',{@keplar_button_Callback});
hinterplan = uicontrol('Style','pushbutton',...
'String','Interplantary',...
'Position',[600,scrsz(4)-70,70,25],...
'Callback',{@Intbutton_Callback});
    
slider = uicontrol('Style', 'slider',...
        'Min',15,'Max',300,'Value',30,...
        'Position', [400 scrsz(4)-130 120 15]); %animation speed slider

slider_title = uicontrol('Units','pixels','Position',[400 scrsz(4)-160 120 20],...
         'String','Animation Speed','style','text','backgroundcolor','k',...
         'foregroundcolor','w','fontsize',9); 
     
strings = {'a','e','i','RAAN','AOP','TA(v)'};
tags = {'a','e','i','RAAN','AOP','v'};
values = {'8000','0.1','30','345','115','80'};
strings1 = {'R','V'};
tags1 = {'R','V'};
values1 = {'[7016 5740 638]','[0.24 -0.79 -7.11]'};
strings2 = {'TLE','Planet name mars or venus','Earth parking orbit raduis in Km','Planet parking orbit raduis in Km'};
tags2 = {'TLE','planet','R_p_e','R_p_p'};
values2 = {'2 23455  99.0090 272.6745 0008546 223.1686 136.8816 14.11711747148495','mars','6600','6400'};
camzoom(0.5);
set([f,hR_V,h6OEs,hTLE,orbit_plot,hinterplan,slider_title,slider,keplar,num],...
'Units','normalized');
set(f,'Name','Orbit Propagator')

global ppp pppp slv;

for kk = 1:length(strings),
      ppp(kk) = uicontrol('Units','pixels','Position',[2 (20*(kk-1)+scrsz(4)-210) 50 20],...
         'String',strings{kk},'style','text','backgroundcolor','k',...
         'foregroundcolor','w','fontsize',9);
     pppp(kk)= uicontrol('Units','pixels','Position',[55 (20*(kk-1)+scrsz(4)-210) 50 20],...
         'tag',tags{kk},'style','edit','string',values{kk},'backgroundcolor',[1 1 1]);
    

end

ppp(7)= uicontrol('Units','pixels','Position',[125 (20*(5)+scrsz(4)-210) 50 20],...
         'String',strings1{1},'style','text','backgroundcolor','k',...
         'foregroundcolor','w','fontsize',9,'Visible','on');
      pppp(7)=uicontrol('Units','pixels','Position',[180 (20*(5)+scrsz(4)-210) 100 20],...
         'tag',tags1{1},'style','edit','string',values1{1},'backgroundcolor',[1 1 1]);

    ppp(8)= uicontrol('Units','pixels','Position',[125 (20*(4)+scrsz(4)-210) 50 20],...
         'String',strings1{2},'style','text','backgroundcolor','k',...
         'foregroundcolor','w','fontsize',9);
   pppp(8)=   uicontrol('Units','pixels','Position',[180 (20*(4)+scrsz(4)-210) 100 20],...
         'tag',tags1{2},'style','edit','string',values1{2},'backgroundcolor',[1 1 1]);
     ppp(9)=   uicontrol('Units','pixels','Position',[280 (20*(5)+scrsz(4)-210) 50 20],...
         'String',strings2{1},'style','text','backgroundcolor','k',...
         'foregroundcolor','w','fontsize',9);
     pppp(9)= uicontrol('Units','pixels','Position',[335 (20*(5)+scrsz(4)-210) 400 20],...
         'tag',tags2{1},'style','edit','string',values2{1},'backgroundcolor',[1 1 1]);
     ppp(10) = uicontrol('Units','pixels','Position',[40,scrsz(4)-260,80,20],...
         'String','TOF in Hours','style','text','backgroundcolor','k',...
         'foregroundcolor','w','fontsize',9);
     pppp(10)= uicontrol('Units','pixels','Position',[125,scrsz(4)-260,20,20],...
         'tag','TOF1','style','edit','string','2','backgroundcolor',[1 1 1]);
 
set([ppp,pppp],...
 'Units','normalized','Visible','off') 
set([ppp(10),pppp(10)],...
'Visible','on')   

%%Graphics 
image_file = 'BlueMarble1.bmp'; %Texture 2d image
%earth dimensions
erad    = 6371008.7714*10^-3; % equatorial radius (meters)
prad    = 6371008.7714*10^-3; % polar radius (meters)
erot    = 7.2921158553e-5; % earth rotation rate (radians/sec)
npanels = 180;   % Number of globe panels around the equator deg/panel = 180/npanels
alpha   = 1; %opaque
GMST0 = 4.89496121282306;
wa=0; %counter
wd=0; %counter
count = 0; %counter
step=100; %time step required for solving
cdata = imread(image_file); 



set(gca, 'NextPlot','add', 'Visible','off','dataaspectratio',[1 1 1]);
view(30,10);%camera view

%%plane construction
pointA = [-10000,-10000,0];
pointB = [-10000,10000,0];
pointC = [10000,10000,0];
pointD = [10000,-10000,0];
points=[pointA' pointB' pointC' pointD']; % using the data given in the question
fill3(points(1,:),points(2,:),points(3,:),[0.9 0.8 0.6],'FaceAlpha', 0.7);


%Earth construction
[x, y, z] = ellipsoid(0, 0, 0, erad, erad, prad, npanels); %generate the mesh
hold on
globe = surf(x, y, -z, 'FaceColor', 'none', 'EdgeColor', 0.5*[1 1 1]); %sphere drawing
axis equal
xl = xlim(); %x_axis
yl = ylim(); %y_axis
zl = zlim(); %z_axis
hold on;
%Earth-centered inertial (ECI)axes
 x_arrow=line(2*xl, [0,0], [0,0], 'LineWidth', 1, 'Color', [0.56, 0.90, 0.93]);
 y_arrow=line([0,0], 2*yl, [0,0], 'LineWidth', 1, 'Color', [0.7, 0, 0]);
 z_arrow=line([0,0], [0,0], 2*zl, 'LineWidth', 1, 'Color', [1, 1, 0.76]);

if ~isempty(GMST0)
    hgx = hgtransform;
    set(hgx,'Matrix', makehgtform('zrotate',GMST0));
    set(globe,'Parent',hgx);
end
hold on

%earth texture
set(globe, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', alpha, 'EdgeColor', 'none');
set(f,'Visible','on');



%%functions Callback
    function num_button_Callback(source,eventdata)
        slv = 2; %slv=2 for numerical solution
    end

    function keplar_button_Callback(source,eventdata)
        slv = 1; %slv=2 for keplar equations solution
    end

function OEbutton_Callback(source,eventdata)
    figure(f);
    %hide unnecessary elements
  mqq=find(ishandle(ppp)~=0);
 set([ppp(mqq),pppp(mqq)],...
 'Visible','off')
 set([ppp(1:6),pppp(1:6),num,keplar,ppp(10),pppp(10)],...
 'Visible','on')
end

function R_Vbutton_Callback(source,eventdata)
    figure(f);
    %hide unnecessary elements
   mqq=find(ishandle(ppp)~=0);
    set([ppp(mqq),pppp(mqq)],...
 'Visible','off')
 set([ppp(7:8),pppp(7:8),num,keplar,ppp(10),pppp(10)],...
'Visible','on')
end

function TLEbutton_Callback(source,eventdata)
    figure(f);
    %hide unnecessary elements
    mqq=find(ishandle(ppp)~=0);
    if length(mqq)>1
set([ppp(mqq),pppp(mqq)],...
 'Visible','off')
 set([ppp(9),pppp(9),num,keplar,ppp(10),pppp(10)],...
 'Visible','on')
    end
end


function Intbutton_Callback(source,eventdata) %Interplanetary function
    
    f2 = figure('Visible','off','Position',scrsz,'Color','k'); %create new figure
    set(f2,'Name','Interplanetary'); %name it as Interplanetary
   %%interplanetary figure construction
    Start = uicontrol('Style','pushbutton',...
'String','Start',...
'Position',[200,scrsz(4)-70,70,25],...
'Callback',{@startbutton_Callback});
 dd(1)=   uicontrol('Units','pixels','Position',[1 (20*(5)+scrsz(4)-210) 300 20],...
         'String',strings2{2},'style','text','backgroundcolor','k',...
         'foregroundcolor','w','fontsize',8);
     ddd(1)= uicontrol('Units','pixels','Position',[250 (20*(5)+scrsz(4)-210) 50 20],...
         'tag',tags2{2},'style','edit','string',values2{2},'backgroundcolor',[1 1 1]);
     dd(2)=   uicontrol('Units','pixels','Position',[1 (20*(3)+scrsz(4)-210) 300 20],...
         'String',strings2{3},'style','text','backgroundcolor','k',...
         'foregroundcolor','w','fontsize',8);
     ddd(2)= uicontrol('Units','pixels','Position',[250 (20*(3)+scrsz(4)-210) 50 20],...
         'tag',tags2{3},'style','edit','string',values2{3},'backgroundcolor',[1 1 1]);
     dd(3)=   uicontrol('Units','pixels','Position',[1 (20*(1)+scrsz(4)-210) 300 20],...
         'String',strings2{4},'style','text','backgroundcolor','k',...
         'foregroundcolor','w','fontsize',8);
     ddd(3)= uicontrol('Units','pixels','Position',[250 (20*(1)+scrsz(4)-210) 50 20],...
         'tag',tags2{4},'style','edit','string',values2{4},'backgroundcolor',[1 1 1]);
     haa = axes('Units','Pixels','Position',[1,1,scrsz(3),scrsz(4)],'Color','k');

set([dd,ddd,Start,haa],...
'Units','normalized');
set(f2,'Visible','on');

    global f2  
end


function startbutton_Callback(source,eventdata)
    global f2 s1 sun ell venus earth venus1 park_earth park_planet hyper_earth hyper_planet globe1

    count=count+1;
   if count > 1
  set(ell, 'Visible', 'off');
set(earth, 'Visible', 'off');
set(venus1, 'Visible', 'off');
set(park_earth, 'Visible', 'off');
set(park_planet, 'Visible', 'off');
set(hyper_earth, 'Visible', 'off');
set(hyper_planet, 'Visible', 'off');
set(globe1, 'Visible', 'off');
set(s1, 'Visible', 'off');
set(venus, 'Visible', 'off');
set(sun, 'Visible', 'off');
    end
   
    
    %getting data from the user
planet = get(findobj('tag','planet'),'string');
R_p_e = str2num(get(findobj('tag','R_p_e'),'string'));
R_p_p = str2num(get(findobj('tag','R_p_p'),'string'));
R_e = 1.496*10^8;
m_s = 1.327*10^11;
m_e = 398600;
 
%the choosen planet
switch planet
    case 'mars'
image_file_venus = 'mars.bmp';
m_p = 43050;
R_p = 2.278*10^8;
    case 'venus'
        image_file_venus = 'venus.bmp';
       m_p = 3.257*10^5;
       R_p = 1.081*10^8;
    end



%calculations
%1st region
a_trans = (R_e+R_p)/2;
E_trans = -m_s/(2*a_trans);


V_e = (m_s/R_e)^(1/2);
V_trans_e = (2*(m_s/R_e + E_trans))^(1/2);
V_inf_e = norm(V_trans_e - V_e);

V_m = (m_s/R_p)^(1/2);
V_trans_m = (2*(m_s/R_p + E_trans))^(1/2);
V_inf_m = norm(V_trans_m - V_m);

%second region
V_p_e = (m_e/R_p_e)^(1/2);
V_trans_e = (2*(m_e/R_p_e + (V_inf_e^2)/2))^(1/2);
AV1 = norm(V_trans_e - V_p_e);
%third region
V_p_m = (m_p/R_p_p)^(1/2);
V_trans_m = (2*(m_p/R_p_p + (V_inf_m^2)/2))^(1/2);
AV2 = norm(V_trans_m - V_p_m);

AV = AV1+AV2;

%Time and phase 
TOF = (pi)*((a_trans^3)/m_s)^(1/2);
T_O_F = TOF/(30*24*60^2);
n_p = (m_s/R_p^3)^(1/2);
n_e = (m_s/R_e^3)^(1/2);
alpha_lead = n_p*TOF;
phi = 180 - 180*alpha_lead/pi;
Time1 = 2*pi/((norm(n_p-n_e))*(30*24*60^2));
Time2 = 2*Time1;
AV1T = num2str(AV1);
AV2T = num2str(AV2);
phiT = num2str(phi);
Time1T = num2str(Time1);
Time2T = num2str(Time2);
TOFT = num2str(T_O_F);
% strings3 = {'AV1','AV2','phi','TOFT','Time1','Time2'};
% values3 =  {AV1T,AV2T,phiT,TOFT,Time1T,Time2T};

    


%simulation calulations

if R_e>R_p
    R_per = R_p;
    R_ap = R_e;
else
    R_per = R_e;
    R_ap = R_p;
end

%transfer orbit calculations
e1=double(solve(R_per == a_trans*(1-sym('e')^2)/(1+sym('e')),sym('e')));
c = a_trans*e1;
b =(-c^2+ a_trans^2)^(1/2);
n_park_p = (m_p/R_p_p^3)^(1/2);
n_park_e = (m_e/R_p_e^3)^(1/2);
n_trans = (m_s/a_trans^3)^(1/2);

%graphics
space_color = 'k';

image_file_sun='sun.bmp'; %sun Image
image_file = 'BlueMarble1.bmp'; %Texture 2d image


%% Create figure
hold on;
set(gca, 'NextPlot','add', 'Visible','off','dataaspectratio',[1 1 1]);% Turn off the normal axes
view(0,90); % Set initial view
axis vis3d; %set 3d axes
axis([-1.5*a_trans 1.5*a_trans -1.7*a_trans 1.7*a_trans ]); %set axes limits
%% Create wireframe globe
coords_ell = calculateEllipse(0, 0, a_trans, b, 0, 300);
coords_ell = coords_ell(1:150,:);

%hyper_earth_orbit
a = R_e;
b1 = coords_ell(1,2)/(coords_ell(1,1)^2/a^2 - 1)^(1/2);
x5 = linspace(a,coords_ell(1,1),5);
y5 = b1*(x5.^2/a^2 -1).^(1/2);

%hyper_planet_orbit
a1 = -R_p;
b2 = coords_ell(150,2)/(coords_ell(150,1)^2/a1^2 - 1)^(1/2);
x6 = linspace(coords_ell(150,1),a1,5);
y6 = b2*(x6.^2/a1^2 -1).^(1/2);

% Create a 3D meshgrid of the spheres points using the ellipsoid function
%orbits mesh
coords_earth = calculateEllipse(0, 0, R_e, R_e, 0, 300);
coords_venus = calculateEllipse(0, 0, R_p, R_p, 0, 300);
coords_park_earth = calculateEllipse(R_e,0,R_p_e*10^3,R_p_e*10^3,0,50);
coords_park_planet = calculateEllipse(R_p,0,R_p_p*10^3,R_p_p*10^3,0,ceil(50*n_park_p/n_park_e));
 %Earth
 %earth dimensions
erad    = 6371008.7714*10^-3; % equatorial radius (meters)
prad    = 6371008.7714*10^-3; % polar radius (meters)
erot    = 7.2921158553e-5; % earth rotation rate (radians/sec)
npanels = 180;   % Number of globe panels around the equator deg/panel = 180/npanels
alpha   = 1; %opaque
GMST0 = 4.89496121282306;
[x, y, z] = ellipsoid(0, R_e, 0, erad*10^3, erad*10^3, prad*10^3, npanels);
hold on
globe1 = surf(x, y, -z, 'FaceColor', 'none', 'EdgeColor', 0.5*[1 1 1]);
hold all

%venus
[x2, y2, z2] = ellipsoid(0, R_p, 0, 6052*10^3, 6052*10^3, 6052*10^3, npanels);
venus = surf(x2, y2, -z2, 'FaceColor', 'none', 'EdgeColor', 0.5*[1 1 1]);

%sun
[x3, y3, z3] = ellipsoid(0, 0, 0, 69570000, 69570000, 69570000, npanels);
sun = surf(x3, y3, -z3, 'FaceColor', 'none', 'EdgeColor', 0.5*[1 1 1]);

%transfer orbit
 ell=plot3(coords_ell(:,1),coords_ell(:,2),zeros(length(coords_ell),1));
 earth=plot3(coords_earth(:,1),coords_earth(:,2),zeros(length(coords_earth),1));
 venus1=plot3(coords_venus(:,1),coords_venus(:,2),zeros(length(coords_venus),1));
 %parking orbits
park_earth = plot3(coords_park_earth(:,1),coords_park_earth(:,2),zeros(length(coords_park_earth),1));
park_planet = plot3(coords_park_planet(:,1),coords_park_planet(:,2),zeros(length(coords_park_planet),1));
%hyperbolic  orbits
hyper_earth = plot3(x5,y5,zeros(length(x5)));
 hyper_planet = plot3(x6,y6,zeros(length(x6)));

 
 %Colors
set(ell, 'Color', [0.8 ,0 ,0]);
set(earth, 'Color', [0 ,0 ,0.8]);
set(venus1, 'Color', [0 ,0.8 ,0]);
set(park_earth, 'Color', [0.8 ,0 ,0]);
set(park_planet, 'Color', [0 ,0 ,0.8]);
set(hyper_earth, 'Color', [0.8 ,0.8 ,0]);
set(hyper_planet, 'Color', [0.8 ,0.5 ,0.2]);

if ~isempty(GMST0)
    hgx = hgtransform;
    set(hgx,'Matrix', makehgtform('zrotate',GMST0));
    set(globe1,'Parent',hgx);
    set(venus,'Parent',hgx);
    set(sun,'Parent',hgx);
end

%% Texturemap the globe,Sun and the planet

% Load Earth image for texture map
cdata = imread(image_file);
cdata_venus = imread(image_file_venus);
cdata_sun = imread(image_file_sun);
% Set image as color data (cdata) property, and set face color to indicate
% a texturemap, which Matlab expects to be in cdata. Turn off the mesh edges.
hold on
set(globe1, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', 1, 'EdgeColor', 'none');
set(venus, 'FaceColor', 'texturemap', 'CData', cdata_venus, 'FaceAlpha', alpha, 'EdgeColor', 'none');
set(sun, 'FaceColor', 'texturemap', 'CData', cdata_sun, 'FaceAlpha', alpha, 'EdgeColor', 'none');

switch (planet)
    case 'mars'
        alph_e = n_e/n_trans; %1 degree every 1 loop
alph_p = n_p/n_trans; %n_p/n_e degree every 1 loop
alph_p_planet = (phi+21.2)*alph_p ;
alph_p_park = (phi+35)*alph_p;
i=length(coords_park_planet)/2;
    case 'venus'
         alph_e = n_e/n_trans; %1 degree every 1 loop
alph_p = n_p/n_trans; %n_p/n_e degree every 1 loop
alph_p_planet = (phi+27)*alph_p ;
alph_p_park = (phi+35)*alph_p;
i=1;
end
l=1;
r_c = 2000000; %sat_raduis
rotate(globe1,[0,0,1],-10.5);
rotate(venus,[0,0,1],alph_p_planet);
rotate(park_planet,[0,0,1],alph_p_park);
rotate(sun,[0,1,0],-90);
set(globe1, 'Visible','on');
set(f2,'Visible','on');
%% animation loop
while l>0
   
    %first phase (hyper_earth)
    if l<6
x0 = x5(l); y0 = y5(l); z0 = 0;
[x1, y1, z1] = sphere(50);
x1 = x1*r_c + x0;
y1 = y1*r_c + y0;
z1 = z1*r_c + z0;
s1 = surf(x1,y1,z1);
    end
     %second phase (transfer orbit)
    if l>5 && l < length(coords_ell)+5
x0 = coords_ell(l-5,1); y0 = coords_ell(l-5,2); z0 = 0;
[x1, y1, z1] = sphere(50);
x1 = x1*r_c + x0;
y1 = y1*r_c + y0;
z1 = z1*r_c + z0;
s1 = surf(x1,y1,z1);
    end
    
 %third phase (hyper_planet)
    if l>length(coords_ell)+4 && l<length(coords_ell)+9
x0 = x6(l-length(coords_ell)-4); y0 =y6(l-length(coords_ell)-4); z0 = 0;
[x1, y1, z1] = sphere(50);
x1 = x1*r_c + x0;
y1 = y1*r_c + y0;
z1 = z1*r_c + z0;
s1 = surf(x1,y1,z1);
end
     %fourth phase (planet_park)
     if l>length(coords_ell)+8    
x0 = coords_park_planet(i,1); y0 =coords_park_planet(i,2); z0 = 0;
[x1, y1, z1] = sphere(50);
x1 = x1*r_c + x0;
y1 = y1*r_c + y0;
z1 = z1*r_c + z0;
s1 = surf(x1,y1,z1);
rotate(s1,[0,0,1],l*alph_p+alph_p_park);
i=i+1;
if i>length(coords_park_planet)
    i=1;
end
     end

set(s1, 'EdgeColor', [0 ,0.75 ,0.75]);
rotate(globe1,[0,0,1],alph_e);
rotate(park_earth,[0,0,1],alph_e);
rotate(venus,[0,0,1],alph_p);
rotate(park_planet,[0,0,1],alph_p);
rotate(sun,[0,0,1],5);
 

  pause(0.1);
      
    if ishandle(globe1)
        
    else
       count=0;
      break
    end
    
l = l+1;
     set(s1,'Visible','off');
 
end
    
end

    function orbit_plotbutton_Callback(source,eventdata) %plot button Callback
        wa=wa+1;
        global a e i v RAAN AOP rot f2 s Tp n M
        if wa>1
            set(s,'Visible','off');
        end
  Mu = 398600;
   % get orbital elements from GUI
   if length(findobj(gcf,'style','edit','Visible','on'))-1==6
   a = str2num(get(findobj('tag','a'),'string'));
   e = str2num(get(findobj('tag','e'),'string'));
   i = str2num(get(findobj('tag','i'),'string'));
   RAAN = str2num(get(findobj('tag','RAAN'),'string'));
   AOP = str2num(get(findobj('tag','AOP'),'string'));
   v = str2num(get(findobj('tag','v'),'string'));
   Tp = 2*pi*((a^(3))/Mu)^(1/2);
   n = 2*pi/Tp;
   E=acos((e+cosd(v))/(1+e*cosd(v)));
   M=E-e*sin(E);
    end

   if length(findobj(gcf,'style','edit','Visible','on'))-1==2
   R=get(findobj('tag','R'),'string');
   R = str2num(R);
   V = get(findobj('tag','V'),'string');
   V = str2num(V);
   r = norm(R);
v11 = norm(V);
Ene = (v11^2)/2 - Mu/r;
a = -Mu/(2*Ene);
H = cross(R,V);
h = norm(H);
Tp = 2*pi*((a^(3))/Mu)^(1/2);
i = (acos(H(3)/h))*180/pi;
e_v = (1/Mu)*((v11^2 - (Mu/r)).*R - dot(R,V).*V);
e = norm(e_v);
v = (acosd((dot(e_v,R))/(e*r))); 
N = cross([0 0 1], H);
N_scaler = norm(N);
RAAN = 360-(acos(N(1)/N_scaler))*180/pi;
AOP = 360-(acos((dot(e_v,N))/(e*N_scaler)))*180/pi;
n=2*pi/Tp;
E=acos((e+cosd(v))/(1+e*cosd(v)));
   M=E-e*sin(E);
   end

   if length(findobj(gcf,'style','edit','Visible','on'))-1==1
   TLE = get(findobj('tag','TLE'),'string');
n = str2double(TLE(53:63)); %rev./day
n = n*(2*pi/(24*60^2)); % rad/sec
Tp = 2*pi/n; %sec
a = (Mu/n^2)^(1/3); %from center of the earth
i = str2double(TLE(9:16)); %deg
RAAN = str2double(TLE(18:25)); %deg
e = TLE(27:33);
power = length(e);
e = str2double(e);
e = e*10^(-power);  %decimal
AOP = str2double(TLE(35:42)); %deg
M = str2double(TLE(44:51)); 
M = M*pi/180; %rad
E = double(solve(M == sym('tt')-e*sin(sym('tt')),sym('tt')));
v = acosd((cos(E)-e)/(1-e*cos(E))); %intialy in deg
   end 

   TOF = 0:step:Tp+200;
    p=a*(1-e^2);
    Q_Xx = [cosd(RAAN)*cosd(AOP)-sind(RAAN)*sind(AOP)*cosd(i) sind(RAAN)*cosd(AOP)+cosd(RAAN)*cosd(i)*sind(AOP)  sind(i)*sind(AOP);-cosd(RAAN)*sind(AOP)-sind(RAAN)*cosd(AOP)*cosd(i) -sind(RAAN)*sind(AOP)+cosd(RAAN)*cosd(i)*cosd(AOP)  sind(i)*cosd(AOP); sind(RAAN)*sind(i) -cosd(RAAN)*sind(i) cosd(i)];

   if slv==1;  %Solving with keplar 
       M2 = n*(TOF)+ M;
for u = 1:length(M2)
while (M2(u)>2*pi)
    M2(u) = M2(u) - 2*pi;
end
E2(u) = double(solve(M2(u) == sym('x')-e*sin(sym('x')),sym('x')));
v1(u) = acos((cos(E2(u))-e)/(1-e*cos(E2(u)))); %rad
if E2(u)>pi
    v1(u) = 2*pi-v1(u);
end
r_new(u) = double(p/(1+e.*cos(v1(u)))); %mag.
V_PQnew(u,:) = transpose(((Mu/p)^(1/2))*[ -sin(v1(u)); e+cos(v1(u)) ;0]); %Vector in PQ
R_PQnew(u,:) = transpose(r_new(u)*[ cos(v1(u)) ; sin(v1(u)) ; 0]); %Vector in PQ
R_IJKnew(u,:) = transpose(double(transpose(Q_Xx)*transpose(R_PQnew(u,:))));
V_IJKnew(u,:) = transpose(double(transpose(Q_Xx)*transpose(V_PQnew(u,:))));

%ground tracks in keplar soultion
ne=7.2921*10^-5;
d2(u)=dot([0 0 1],R_IJKnew(u,1:3))/norm(R_IJKnew(u,1:3));
lat(u) = 90 - acosd(d2(u));
d1(u)=dot([1 0],R_IJKnew(u,1:2))/norm(R_IJKnew(u,1:2));
lon(u) = acosd(d1(u));
if R_IJKnew(u,2) < 0
    lon(u) = -lon(u);
end
lon(u)=lon(u)-ne*180*step*u/pi;
if lon(u)<-180
    lon(u)=(lon(u)+360);
end
end
TOF1 = 60^2*str2num(get(findobj('tag','TOF1'),'String'));
Mfinal= TOF1*n+M;
Efinal = double(solve(Mfinal == sym('x')-e*sin(sym('x')),sym('x')));
vfinal = acos((cos(Efinal)-e)/(1-e*cos(Efinal))); %rad

if Efinal>pi
    vfinal = 2*pi-vfinal;
end
r_final= double(p/(1+e.*cos(vfinal))); %mag.
V_PQnew= transpose(((Mu/p)^(1/2))*[ -sin(vfinal); e+cos(vfinal) ;0]); %Vector in PQ
R_PQnew = transpose(r_final*[ cos(vfinal) ; sin(vfinal) ; 0]); %Vector in PQ
R_IJKnew1 = transpose(double(transpose(Q_Xx)*transpose(R_PQnew)));
V_IJKnew1= transpose(double(transpose(Q_Xx)*transpose(V_PQnew)));
   end
   
   
   
   if slv==2; %numerical Intergrator
 r = p/(1+e*cosd(v));
V_PQ = ((Mu/p)^(1/2))*[ -sind(v); e+cosd(v) ;0];
R_PQ = r*[ cosd(v) ; sind(v) ; 0];
R_IJK0 = transpose(Q_Xx)*R_PQ;
V_IJK0 = transpose(Q_Xx)*V_PQ;

options = odeset('RelTol',1e-4,'AbsTol',[1e-8 1e-8 1e-8 1e-5 1e-5 1e-5]);
[T,R_V_IJKnew] = ode45(@TwoBody,TOF,[R_IJK0 V_IJK0],options); 
       ne=7.2921*10^-5;
       
       TOF1 = 60^2*str2num(get(findobj('tag','TOF1'),'String'));
       TOF1 = 0:15:TOF1;
       [T,R_V_IJKnew1] = ode45(@TwoBody,TOF1,[R_IJK0 V_IJK0],options); 
  R_IJKnew1=R_V_IJKnew1(length(R_V_IJKnew1),1:3);
  V_IJKnew1=R_V_IJKnew1(length(R_V_IJKnew1),4:6);

for u = 1:length(R_V_IJKnew)
%ground tracks in numerical soultion
d2(u)=dot([0 0 1],R_V_IJKnew(u,1:3))/norm(R_V_IJKnew(u,1:3));
lat(u) = 90 - acosd(d2(u));
d1(u)=dot([1 0],R_V_IJKnew(u,1:2))/norm(R_V_IJKnew(u,1:2));
lon(u) = acosd(d1(u));

if R_V_IJKnew(u,2) < 0
    lon(u) = -lon(u);
end
lon(u)=lon(u)-ne*180*step*u/pi;
if lon(u)<-180
    lon(u)=(lon(u)+360);
end
R_IJKnew = R_V_IJKnew(:,1:3);
V_IJKnew = R_V_IJKnew(:,4:6);
end
   end

R_IJKnew1 = num2str(R_IJKnew1(1,:)); %output R
V_IJKnew1 = num2str(V_IJKnew1(1,:)); %output V
   hold on
   axis([-a a -3*a 3*a ]);
   hold on
 orbit = plot3(R_IJKnew(:,1),R_IJKnew(:,2),R_IJKnew(:,3),'Color',[rand(1,1) rand(1,1) rand(1,1)]); %orbit plot 
 if ishandle(orbit)
    wd=wd+1;
end
 if wd==1;
     camzoom(2);
 end
 
 %output text boxes for R & V
 output(1) = uicontrol('Units','pixels','Position',[60,scrsz(4)-360,250,15],...
         'String',R_IJKnew1,'style','text','backgroundcolor',[0.1 0.1 0.1],...
         'foregroundcolor','w','fontsize',9);
     output(2) = uicontrol('Units','pixels','Position',[60,scrsz(4)-390,250,15],...
         'String',V_IJKnew1,'style','text','backgroundcolor',[0.1 0.1 0.1],...
         'foregroundcolor','w','fontsize',9);
     output(3) = uicontrol('Units','pixels','Position',[20,scrsz(4)-395,20,20],...
         'String','V','style','text','backgroundcolor','k',...
         'foregroundcolor','w','fontsize',9);
      output(4) = uicontrol('Units','pixels','Position',[20,scrsz(4)-365,20,20],...
         'String','R','style','text','backgroundcolor','k',...
         'foregroundcolor','w','fontsize',9);
     set([output],...
'Units','normalized','Visible','on');
 
%% Ground tracks Figure
cdata = imread(image_file);
f3=figure(3);
set(f3,'Name','Ground Tracks');
movegui(f3,'east');
clf('reset') %clear all old objects

 h2=surface([-180 180; -180 180], [-90 -90; 90 90], [-1 -1; -1 -1], ...
    'FaceColor', 'texturemap','FaceAlpha', 0.9, 'CData', cdata ); %plane for 2d map
rotate(h2,[1 0 0],180); %flipping the image as a correction
axis([-180 180 -90 90]);
w=0;%counter
bb=0;%counter
%seperation of data as it the track get out from one side and enter from
%the other
for u=1:length(lon)-1
    if abs(lon(u+1) - lon(u)) > 200
        bb=bb+1;
           w(bb) = u+1;
      end
end
if bb~=0
    w(bb+1)=length(lon+1);
    zz=1;
    for ii=1:bb+1
   A{ii}=lon(zz:w(ii)-1);
   B{ii}=lat(zz:w(ii)-1);  
   zz = w(ii);
   hold on
  plot3(A{ii},B{ii},zeros(1,length(A{ii})),'Color',[1 1 0]); %ground track plot as it get out
    end
else
    hold on
    plot3(lon,lat,zeros(1,length(lon)),'Color',[1 1 0]); %ground track plot if it dosen't get out from the map
end 

%figure &axes are on
 set(gca, 'Visible','on');
 grid minor
 axh = gca;
set(axh,'XGrid','on')
set(axh,'YGrid','on')
j=1;
l=0;
r_c = 300; %sat raduis

%animation of satallite and Earth
 figure(f);
while j>0
    rot = get(findobj('style','slider'),'Value');
alph = ne*180*rot/(pi);  
alphsat=ceil(rot/step);
      l=l+1;
      if l*alphsat > length(R_IJKnew)-alphsat
          k=length(R_IJKnew);
        l=2;
x0 = R_IJKnew(k,1); y0 = R_IJKnew(k,2); z0 = R_IJKnew(k,3);
      else
          x0 = R_IJKnew(l*alphsat,1); y0 = R_IJKnew(l*alphsat,2); z0 = R_IJKnew(l*alphsat,3);
      end
[x1, y1, z1] = sphere(50);
x1 = x1*r_c + x0;
y1 = y1*r_c + y0;
z1 = z1*r_c + z0;

   if gcf == f 
  rotate(globe,[0,0,1],alph);
    s = surf(x1,y1,z1); %satallite shpere
 set(s, 'EdgeColor', [0 ,0.75 ,0.75]);
   end
pause(0.1);
      if ishandle(s) 
     set(s,'Visible','off');
     else
      break
    end
    j = j+1;
end
end

end

