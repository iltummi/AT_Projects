%% Progetto ESE

clc
clear all
close all
tic

addpath(genpath("M-Files for Orbital Mechanics for Engineering Students, 3e"));
addpath(genpath("ESE_functions"));
addpath(genpath("images"));
%% Initialization

%Sun gravitational parameter
global mu
mu = 1.327565122000000e+11; %[km^3/s^2]

%   planet_id - planet identifier:
%                1 = Mercury
%                2 = Venus
%                3 = Earth
%                4 = Mars
%                5 = Jupiter
%                6 = Saturn
%                7 = Uranus
%                8 = Neptune
%                9 = Pluto
%               10 = Enceladus
%               11 = Sun

colors = ["g"          %green
      "m"          %magenta
      "b"          %blue
      "r"          %red
      "#A2142F"    %darker red
      "#7E2F8E"    %purple
      "#4DBEEE"    %darker cyan
      "c"          %(bright) cyan
      "#D95319"    %orange
      "#77AC30"    %darker green
      "#D95319"];  %orange, not visible due to Sun orbit dimensions

[xx,yy,zz] = sphere(10);

masses = 10^24 * [0.330104
                  4.86732
                  5.97219
                  0.641693
                  1898.13
                  568.319
                  86.8103
                  102.410
                  0.01309
                  8.6*1e-5
                  1989100]; %[kg]

radii = [2439.7
         6051.8 
         6371
         3389.5
         69911
         58232
         25362
         24622
         1151
         249.9
         695508]; %[km] 
     
distances = [57909227
             108209475
             149598262
             227943824
             778340821
             1426666422
             2870658186
             4498396441
             1426904442
             413690250
                     0];%[km]

G = 6.6742e-20; %[km^3/kg/s^2]
    
%% Data
%Distances from the Sun
Earth_to_Sun = distances(3);%149598262; %[km]
Venus_to_Sun = distances(2);%108209475; %[km]
Saturn_to_Sun = distances(6);%1426666422 %[km]
Enc_to_Sun = distances(10); %1426904442[km]
Enc_to_Sat = 238020; %[km]

%Masses
Earth_mass = masses(3);%2.589*10^20; %[kg]
Venus_mass = masses(2);%8.958*10^20; %[kg]
Saturn_mass= masses(6); %[Kg]
Enc_mass= masses(10); %[Kg]
Sun_mass = masses(11);%1.989*10^30; %[kg]

%SOI: (m_planet/m_Sun)^(2/5) * distance_from_Sun
Earth_SOI = (Earth_mass/Sun_mass)^(2/5)*Earth_to_Sun;%9.24*10^5; %[km]
Venus_SOI = (Venus_mass/Sun_mass)^(2/5)*Venus_to_Sun;%6.1620*10^5; %[km]
Saturn_SOI = (Saturn_mass/Sun_mass)^(2/5)*Saturn_to_Sun; %5.4538*10^7 %[km]
Enceladus_SOI = (Enc_mass/Saturn_mass)^(2/5)*Enc_to_Sat; %4.452213284110794e+02 %[km]


%Parking orbits (not considering body radius)
Epark_radius = 200; %[km]
Epark_inclination = 0; %[rad]
Spark_radius=400; %[km]
Encpark_radius=100; %[km]

%%%Date%%%
%Data Parcheggio 02/07/2022
data_p.year=2022;
data_p.month=7;
data_p.day=2;

%Data 0 02/08/2022 Departure from Terra.
data0.year=2022;
data0.month=8;
data0.day=2;

%Data 1 02/12/2022 Fly-by on Venere.
data1.year=2022;
data1.month=12;
data1.day=2;

%Data 2 02/05/2023 Fly-by on Earth and travel to Saturn.
data2.year=2023;
data2.month=5;
data2.day=2;

%Data 3 02/02/2028  Entering Saturn parking orbit.
data3.year=2028;
data3.month=2;
data3.day=2;

%Data 4 28/05/2028 Exiting Saturn parking orbit and travel to Enceladus.
data4.year=2028;
data4.month=5;
data4.day=28;

%Data 5 29/05/2028 Entering Enceladus parking orbit.
data5.year=2028;
data5.month=5;
data5.day=29;
%% Delta v vectors

delta_v = zeros(4, 1);
% delta_v(1) escape from Earth
% delta_v(2) capture on Saturn
% delta_v(3) escape from Saturn
% delta_v(4) capture on Enceladus


delta_change = zeros(4, 1);
% delta_change(1) change of plane on Earth parking orbit.
% delta_change(2) change of plane on 
% delta_change(3) change of plane on 
% delta_change(4) change of plane on Enceladus parking orbit.


%% Planets configurations
%{
    The mission is subdivided in three parts:
    - Earth to Venus;
    - Venus to Earth;
    - Earth to Saturn;
    - Saturn to Enceladus.
    Variables for each part have the same name but different suffix.
    
    Planet configurations have the following numerical indications:
    - 0: at Earth departure (to Venus) -> 02/08/2022
    - 1: at Venus arrival/departure (for the first fly-by) -> 02/12/2022 
    - 2: at Earth arrival/departure (for the second fly-by) -> 02/05/2023
    - 3: at Saturn arrival (from Earth) -> 02/02/2028
    - 4: at Saturn departure (to Enceladus) -> 28/05/2028 
    - 5: at Enceladus arrival (from Saturn) -> 29/05/2028
%}

%Earth
%02/08/2022
[Earth_coe0, Earth_r0, Earth_v0, ~] =...
                        planet_elements_and_svMOD(3,data0.year,data0.month,data0.day,0,0,0);
%02/12/2022                    
[Earth_coe1, Earth_r1, Earth_v1, ~] =...
                         planet_elements_and_svMOD(3,data1.year,data1.month,data1.day,0,0,0);
%02/05/2023                     
[Earth_coe2, Earth_r2, Earth_v2, ~] =...
                        planet_elements_and_svMOD(3,data2.year,data2.month,data2.day,0,0,0);
%02/02/2028                    
[Earth_coe3, Earth_r3, Earth_v3, ~] =...
                        planet_elements_and_svMOD(3,data3.year,data3.month,data3.day, 0, 0, 0);

%Venus
%02/08/2022
[Venus_coe0, Venus_r0, Venus_v0, ~] =...
                        planet_elements_and_svMOD(2,data0.year,data0.month,data0.day,0,0,0);
%02/12/2022
[Venus_coe1, Venus_r1, Venus_v1, ~] =...
                        planet_elements_and_svMOD(2,data1.year,data1.month,data1.day,0,0,0);
%02/05/2023
[Venus_coe2, Venus_r2, Venus_v2, ~] =...
                        planet_elements_and_svMOD(2,data2.year,data2.month,data2.day,0,0,0);

%Saturn
%02/02/2028
[Saturn_coe3, Saturn_r3, Saturn_v3, ~] =...
                        planet_elements_and_svMOD(6, data3.year,data3.month,data3.day, 0, 0, 0);
%28/05/2028
[Saturn_coe4, Saturn_r4, Saturn_v4, ~] =...
                        planet_elements_and_svMOD(6, data4.year,data4.month,data4.day,0,0,0);
%29/05/2028
[Saturn_coe5, Saturn_r5, Saturn_v5, ~] =...
                        planet_elements_and_svMOD(6, data5.year,data5.month,data5.day,0,0,0);
                    
%Enceladus
%28/05/2028
[Enceladus_coe4, Enceladus_r4, Enceladus_v4, ~] =...
                        planet_elements_and_svMOD(10,data4.year,data4.month,data4.day,0,0,0);
%29/05/2028
[Enceladus_coe5, Enceladus_r5, Enceladus_v5, ~] =...
                        planet_elements_and_svMOD(10,data5.year,data5.month,data5.day,0,0,0);
%% Earth - Venus travel 02/08/2022 -> 02/12/2022
if exist('figure2') == 0  %#ok<*EXIST>
    figure()
else
    figure2()
end

title("Interplanetary orbits")
hold on

%Interplanetary orbit
fprintf('\n\n EARTH TO VENUS \n\n')
%02/08/2022 -> 02/12/2022
[body_pos1, sp_v1, body_posf1, sp_vf1,tof1, orb_elem1] = ...
                gen_orbitMOD(3,2,[data0.year data0.month data0.day 0 0 0],[data1.year data1.month data1.day 0 0 0],0);
            
Ev_orbit = intpl_orbit(tof1,Earth_r0,sp_v1);

%Planet orbits
plot_orbitMOD(3,data0.year)
plot_orbitMOD(2,data1.year)

%Planet positions
plot3(Earth_r0(1),Earth_r0(2),Earth_r0(3),'o','Color',colors(3))
plot3(Earth_r1(1),Earth_r1(2),Earth_r1(3),'x','Color',colors(3))
plot3(Venus_r0(1),Venus_r0(2),Venus_r0(3),'o','Color',colors(4))
plot3(Venus_r1(1),Venus_r1(2),Venus_r1(3),'x','Color',colors(4))

xlabel('x')
ylabel('y')
zlabel('z')
view(-10,45)
grid

%% Venus - Earth travel 02/12/2022 -> 02/05/2023

fprintf('\n\n EARTH TO VENUS \n\n')

%02/12/2022 -> 02/05/2023
[body_pos2, sp_v2, body_posf2, sp_vf2,tof2, orb_elem2] = ...
                gen_orbitMOD(2,3,[data1.year data1.month data1.day,0,0,0],[data2.year data2.month data2.day,0,0,0],0);
Ve_orbit = intpl_orbit(tof2,Venus_r1,sp_v2);

%Planet orbits
plot_orbitMOD(2,data1.year)
plot_orbitMOD(3,data2.year)

%Planet positions
plot3(Earth_r2(1),Earth_r2(2),Earth_r2(3),'d','Color',colors(3))
plot3(Venus_r2(1),Venus_r2(2),Venus_r2(3),'d','Color',colors(4))


%% Earth - Saturn travel 02/05/2023 -> 02/02/2028

fprintf('\n\n EARTH TO SATURN \n\n')

%02/05/2023 -> 02/02/2028
[body_pos3, sp_v3, body_posf3, sp_vf3,tof3, orb_elem3] = ...
                gen_orbitMOD(3,6,[data2.year data2.month data2.day,0,0,0],[data3.year data3.month data3.day,0,0,0],0);
Es_orbit = intpl_orbit(tof3,Earth_r2,sp_v3);

%Planet orbits

plot_orbitMOD(3,data2.year)
plot_orbitMOD(6,data3.year)
%plot_orbitMOD(10, data3.year)

%Planet positions
plot3(Earth_r3(1),Earth_r3(2),Earth_r3(3),'+','Color',colors(3))
plot3(Saturn_r3(1),Saturn_r3(2),Saturn_r3(3),'+','Color',colors(6))
hold on
plot3(Enceladus_r5(1), Enceladus_r5(2), Enceladus_r5(3), '+', 'Color', colors(10));
% view(108, -5)
%% Earth close-up 02/07/2022 -> 02/08/2022
if exist('figure2') == 0
    figure()
else
    figure2()
end

xlabel('x')
ylabel('y')
zlabel('z')
xlim([1.1826*10^8, 1.1829*10^8])
ylim([8.985*10^7, 8.988*10^7])
zlim([-1.5*10^4, 10^4])
view(-10,45)
grid
title("Earth close-up")
hold on

%02/07/2022 -> 02/08/2022
[Earth_esc, delta_v(1)] = escape_hyp_MY(3,Ev_orbit(1:2,1:3),[data0.year data0.month data0.day 0 0 0],...
                                         Epark_radius, orb_elem1, sp_v1);

%02/07/2022 -> 02/08/2022
[park_E0,t_E0] = park_out_MY(3, Earth_r0, Epark_radius, orb_elem1,...
                  Earth_esc(1,1:3), [data_p.year,data_p.month, data_p.day,0,0,0], [data0.year,data0.month,data0.day,0,0,0]);
% %plot change orbit          
% [park_E1,t_E1] = park_out_MY(3, Earth_r0, Epark_radius, zeros(1, 6),...
%                   zeros(1, 3), [data_p.year,data_p.month, data_p.day,0,0,0], [data0.year,data0.month,data0.day,0,0,0]);

[~, ~, delta_change(1)] = ...
    change_of_plane_MY(Epark_inclination, orb_elem1(4), orb_elem1(3), orb_elem1(3), sqrt((Earth_mass*G)/(Epark_radius + radii(3))));

%% Earth SOI close-up
if exist('figure2') == 0
    figure()
else
    figure2()
end

xlabel('x')
ylabel('y')
zlabel('z')
view(-10,45)
grid
title("Earth SOI close-up")
hold on

surface(Earth_r0(1)+Earth_SOI*xx, Earth_r0(2)+Earth_SOI*yy,...
        Earth_r0(3)+Earth_SOI*zz,'FaceColor','none','EdgeColor',colors(3))

%02/07/2022 -> 02/08/2022
hyperbola = escape_hyp_MY(3,Ev_orbit(1:2,1:3),[data0.year data0.month data0.day 0 0 0],...
                        Epark_radius, orb_elem1, norm(sp_vf1));
                    
%02/07/2022 -> 02/08/2022
park_out_MY(3, Earth_r0, Epark_radius, orb_elem1,...
                  Earth_esc(1,1:3), [data_p.year,data_p.month, data_p.day,0,0,0], [data0.year,data0.month,data0.day,0,0,0]);

%% Venus close-up (first fly-by) 02/12/2022

v1=Ev_orbit(end, 4:end);
v2=Ve_orbit(1, 4:end);

[teta, alpha_x,beta_y,gamma_z, R]=angolocompreso(v1,v2);
offset=-100.15;

%02/12/2022
% flyby(planet_id,theta_inf,altitude,flag,year,month,day,hour,minute,second)
[par,fly] = flyby_MY(2, teta,300,1,data1.year,data1.month,data1.day,0,0,0, 0, 0, 0, offset);

% % Per evidenziare il flyby
% figure(1)
% hold on
% xlim([ 5.9*10^7, 6.1*10^7])
% ylim([8.9*10^7, 9.1*10^7])
% zlim([-2.5*10^6, -2.2*10^6])

if exist('figure2') == 0
    figure()
else
    figure2()
end

xlabel('x')
ylabel('y')
zlabel('z')
xlim([ 6.0245*10^7, 6.027*10^7])
ylim([8.9615*10^7, 8.964*10^7])
zlim([-2.254*10^6, -2.24*10^6])
view(-10,45)
grid
title("Venus close-up")
hold on

body_sphere_MY(2,Venus_r1);

plot3(fly(:,1), fly(:,2),fly(:,3),'g-')

[~, ~, delta_change(2)] = change_of_plane_MY(orb_elem1(4), orb_elem2(4), orb_elem1(3), orb_elem2(3), norm(sp_v2)); 

%% Venus SOI close-up
if exist('figure2') == 0
    figure()
else
    figure2()
end

xlabel('x')
ylabel('y')
zlabel('z')
view(-10,45)
grid
title("Venus SOI close-up")
hold on

body_sphere_MY(2,Venus_r1)

surface(Venus_r1(1)+Venus_SOI*xx, Venus_r1(2)+Venus_SOI*yy,...
        Venus_r1(3)+Venus_SOI*zz,'FaceColor','none','EdgeColor',colors(2))
    
plot3(fly(:,1), fly(:,2),fly(:,3),'k')

%% Earth close-up (second fly-by) 02/05/2023
v1=Ve_orbit(end,4:6);
v2=Es_orbit(1,4:6);

[teta2, ~, ~, ~, ~]=angolocompreso(v1,v2);
offset=-100.07;

% flyby(planet_id,theta_inf,altitude,flag,year,month,day,hour,minute,second)
[par2,fly2] = flyby_MY(3,teta2,1000,-1,data2.year,data2.month,data2.day,0,0,0,0,0,0, offset);

%Evidenzia l'iperbole 
figure(1)
hold on
xlim([ 9.05*10^7, 9.3*10^7])
ylim([-1.22*10^8, -1.195*10^8])
zlim([-5*10^4, 5*10^4])

if exist('figure2') == 0
    figure()
else
    figure2()
end

xlabel('x')
ylabel('y')
zlabel('z')
xlim([9.1898e+07, 9.1924e+07])
ylim([-1.2091e+08, -1.2089e+08])
zlim([-6.2687e+03, 1.9215e+04])
view(-10,45)
%view(3)

grid
title("Earth close-up 2")
hold on

body_sphere_MY(3,Earth_r2);

plot3(fly2(:,1), fly2(:,2),fly2(:,3),'g-')

[~, ~, delta_change(3)] = change_of_plane_MY(orb_elem2(4), orb_elem3(4), orb_elem2(3), orb_elem3(3), norm(sp_v3)); 

%% Earth SOI close-up 2
if exist('figure2') == 0
    figure()
else
    figure2()
end

xlabel('x')
ylabel('y')
zlabel('z')
view(-10,45)
grid
title("Earth SOI close-up")
hold on

body_sphere_MY(3,Earth_r2)

surface(Earth_r2(1)+Earth_SOI*xx,Earth_r2(2)+Earth_SOI*yy,...
        Earth_r2(3)+Earth_SOI*zz,'FaceColor','none','EdgeColor',colors(3))
    
plot3(fly2(:,1), fly2(:,2),fly2(:,3),'r')

%% Saturn close-up (arrival) 02/02/2028
if exist('figure2') == 0
    figure()
else
    figure2()
end

xlabel('x')
ylabel('y')
zlabel('z')
xlim([1.187314511476635*10^9, 1.187547439476635*10^9])
ylim([0.713998758042353*10^9, 0.714231686042353*10^9])
zlim([-0.059781141998104*10^9, -0.059548213998104*10^9])
view(-10,45)
grid
title("Saturn close-up (arrival)")
hold on

[Saturn_cap, delta_v(2)] = capture_hyp_MY(6,Es_orbit(end-1:end,1:3),[data3.year data3.month data3.day 0 0 0],...
                        Spark_radius,orb_elem3,sp_vf3);

[park_S,t_S] = park_in_MY(6, Saturn_r3, Spark_radius, orb_elem3,...
                Saturn_cap(end,1:3), [data3.year,data3.month,data3.day,0,0,0], [data4.year,data4.month,data4.day,0,0,0]);

plot3(park_S(end, 1),park_S(end, 2),park_S(end, 3),'o','Color',colors(6))

[~, ~, delta_change(4)] = change_of_plane_MY(orb_elem3(4), 0, orb_elem3(3), orb_elem3(3), sqrt((masses(6)*G)/(Spark_radius+radii(6))));
hold off
%% Saturn close-up (departure) 28/05/2028 -> 29/05/2028
if exist('figure2') == 0
    figure()
else
    figure2()
end

tf_SE=(data5.day-data4.day)*24*3600; %[s]

[park_S4,t_S4] = park_in_MY(6, Saturn_r4, Spark_radius, orb_elem3,...
                Saturn_cap(end,1:3), [data4.year,data4.month,data4.day,0,0,0], [data5.year,data5.month,data5.day,0,0,0]);
hold on
%[orb, t, deltav] = transfer_orbitSE(obj_id, body_pos, r1, r2, tf, grade)
[orb_SE, t_SE, delta_v(3), v_init, v_f]=transfer_orbitSE(6, Saturn_r4, park_S4(1, :) , Enceladus_r5, tf_SE, 'pro');

%Plot planet
body_sphere_MY(10, Enceladus_r5);
body_sphere_MY(6, Saturn_r4);
view(20, 45)
grid
hold off
%% Enceladus close-up 29/05/2028
if exist('figure2') == 0
    figure()
else
    figure2()
end

mu_sat = 3.793074669799999e+07;
%orb_vel=st2vel(orb_SE, t_SE, v_init);

sat2enc_coe=coe_from_sv(orb_SE(1, :), v_f, mu_sat);
[Enceladus_cap, delta_v(4)] = capture_hyp_MY(10,orb_SE(end-1:end, :),[data5.year data5.month data5.day 0 0 0],...
                        Encpark_radius,sat2enc_coe, v_f);
                    
%Data 6
data6.year=2028;
data6.month=6;
data6.day=1;
                    
[park_S5,t_S5] = park_in_MY(10, Enceladus_r5, Encpark_radius, sat2enc_coe,...
                Enceladus_cap(end, :), [data5.year,data5.month,data5.day,0,0,0], [data5.year,data5.month,data5.day,0,28,0]);

[park_S6,t_S6] = park_in_MY_eq(10, Enceladus_r5, Encpark_radius, zeros(1, 6),...
                park_S5(end, :), [data5.year,data5.month,data5.day,1,42,0], [data6.year,data6.month,data6.day,0,0,0]);

plot3(park_S5(end, 1),park_S5(end, 2),park_S5(end, 3),'o','Color',colors(10))
plot3(park_S6(end, 1),park_S6(end, 2),park_S6(end, 3),'o','Color',colors(7))

grid
xlim([1.128676004518143e09 1.128677004118143e09])
ylim([0.795635954495918e09 0.795636954095918e09])
zlim([-0.058740801023396e09 -0.058739801423396e09])
view(-130, 5)
hold off
