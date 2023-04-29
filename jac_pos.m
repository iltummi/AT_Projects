    syms lat lon h vN vE vD ba1  ba2 ba3 ax ay az phi theta psi

    r_e = 6378137; % m %equatorial radius
    r_p = 6356752.3142; %m %polar radius
    e = sqrt(1-(r_p/r_e)^2); % eccentricity

    N = r_e/sqrt(1-(e^2)*(sin(lat))^2); %normal radius
    M = N*(1-e^2)/(1-(e^2)*(sin(lat))^2)^2; %meridian radius

    C_ned_b = [cos(theta)*cos(psi) sin(theta)*sin(phi)*cos(psi)-cos(phi)*sin(psi) sin(theta)*cos(phi)*cos(psi)+sin(phi)*sin(psi);
               cos(theta)*sin(psi) sin(theta)*sin(phi)*sin(psi)+cos(phi)*cos(psi) sin(theta)*cos(phi)*sin(psi)-sin(phi)*cos(psi);
               -sin(theta) sin(phi)*cos(theta) cos(phi)*cos(theta)]; %body to ned
    
    C_ned_e = [-sin(lat)*cos(lon) -sin(lat) cos(lat);
               -sin(lon)  cos(lon) 0;
               -cos(lat)*cos(lon) -cos(lat)*sin(lon) -sin(lat)];

    w_e_ie = [0 0 7.292115*10^-5]; %rad/s velocità ang della terra in ecef
    w_n_ie = C_ned_e*w_e_ie'; %ruoto la vel angolare della terra da ecef a NED

    w_n_en = [vE/(N+h); -vN/(M+h); -vE*tan(lat)/(N+h)];%transport rate

 % State transition function, xdot = f(x, u)
    inv_D = diag([1/(M+h) 1/((N+h)*cos(lat)) -1]);
        f1 = inv_D*[vN;vE;vD];

    pluto = (skew(2*w_n_ie+w_n_en))*[vN;vE;vD];

    f2 = C_ned_b*[ax;ay;az] -pluto + [0;0;9.8] - [ba1;ba2;ba3];  %navigation equation
  
  f = [f1;f2;zeros(3,1)];

  A_pos = [diff(f(1,1),lat) diff(f(1,1),lon) diff(f(1,1),h) diff(f(1,1),vN) diff(f(1,1),vE) diff(f(1,1),vD) diff(f(1,1),ba1) diff(f(1,1),ba2) diff(f(1,1),ba3);
           diff(f(2,1),lat) diff(f(2,1),lon) diff(f(2,1),h) diff(f(2,1),vN) diff(f(2,1),vE) diff(f(2,1),vD) diff(f(2,1),ba1) diff(f(2,1),ba2) diff(f(2,1),ba3);
           diff(f(3,1),lat) diff(f(3,1),lon) diff(f(3,1),h) diff(f(3,1),vN) diff(f(3,1),vE) diff(f(3,1),vD) diff(f(3,1),ba1) diff(f(3,1),ba2) diff(f(3,1),ba3);
           diff(f(4,1),lat) diff(f(4,1),lon) diff(f(4,1),h) diff(f(4,1),vN) diff(f(4,1),vE) diff(f(4,1),vD) diff(f(4,1),ba1) diff(f(4,1),ba2) diff(f(4,1),ba3);
           diff(f(5,1),lat) diff(f(5,1),lon) diff(f(5,1),h) diff(f(5,1),vN) diff(f(5,1),vE) diff(f(5,1),vD) diff(f(5,1),ba1) diff(f(5,1),ba2) diff(f(5,1),ba3);
           diff(f(6,1),lat) diff(f(6,1),lon) diff(f(6,1),h) diff(f(6,1),vN) diff(f(6,1),vE) diff(f(6,1),vD) diff(f(6,1),ba1) diff(f(6,1),ba2) diff(f(6,1),ba3);
           diff(f(7,1),lat) diff(f(7,1),lon) diff(f(7,1),h) diff(f(7,1),vN) diff(f(7,1),vE) diff(f(7,1),vD) diff(f(7,1),ba1) diff(f(7,1),ba2) diff(f(7,1),ba3);
           diff(f(8,1),lat) diff(f(8,1),lon) diff(f(8,1),h) diff(f(8,1),vN) diff(f(8,1),vE) diff(f(8,1),vD) diff(f(8,1),ba1) diff(f(8,1),ba2) diff(f(8,1),ba3);
           diff(f(9,1),lat) diff(f(9,1),lon) diff(f(9,1),h) diff(f(9,1),vN) diff(f(9,1),vE) diff(f(9,1),vD) diff(f(9,1),ba1) diff(f(9,1),ba2) diff(f(9,1),ba3)]




