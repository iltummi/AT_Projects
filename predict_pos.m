function [x_p,P_pos] = predict_pos(x_p, dt, accel, q0, q1, q2, q3, P_pos, Q_pos, N, M)

    ax = accel(1);
    ay = accel(2);
    az = accel(3);

    lat = x_p(1); %sono posizioni di gps 
    lon = x_p(2); %sono posizioni di gps
    h = x_p(3); %sono posizioni di gps
    vN = x_p(4);
    vE = x_p(5);
    vD = x_p(6);
    ba1 = x_p(7);
    ba2 = x_p(8);
    ba3 = x_p(9);

%     C_ned_b = [cos(theta)*cos(psi) sin(theta)*sin(phi)*cos(psi)-cos(phi)*sin(psi) sin(theta)*cos(phi)*cos(psi)+sin(phi)*sin(psi);
%                cos(theta)*sin(psi) sin(theta)*sin(phi)*sin(psi)+cos(phi)*cos(psi) sin(theta)*cos(phi)*sin(psi)-sin(phi)*cos(psi);
%                -sin(theta) sin(phi)*cos(theta) cos(phi)*cos(theta)]; %body to ned
    
    C_ned_b = [q0^2+q1^2-q2^2-q3^2 2*(q1*q2-q0*q3) 2*(q1*q3+q0*q2);
               2*(q1*q2+q0*q3) q0^2-q1^2+q2^2-q3^2 2*(q2*q3-q0*q1);
               2*(q1*q3-q0*q2) 2*(q2*q3+q0*q1) q0^2-q1^2-q2^2+q3^2];
    
    C_ned_e = [-sin(lat)*cos(lon) sin(lat)*sin(lon) cos(lat);
               sin(lon)  cos(lon) 0;
               -cos(lat)*cos(lon) cos(lat)*sin(lon) -sin(lat)];

    w_e_ie = [0 0 7.292115*10^-5]; %rad/s velocità ang della terra in ecef
    w_n_ie = C_ned_e*w_e_ie'; %ruoto la vel angolare della terra da ecef a NED

    w_n_en = [vE/(N+h); -vN/(M+h); -vE*tan(lat)/(N+h)];%transport rate

    pluto = (skew(2*w_n_ie+w_n_en))*[vN;vE;vD];

 % State transition function, xdot = f(x, u)
    inv_D = diag([1/(M+h) 1/((N+h)*cos(lat)) -1]);  
    f1 = inv_D*[vN;vE;vD];

    f2 = C_ned_b*(accel') - pluto + [0;0;9.8] - [ba1;ba2;ba3];  %navigation equation

    f = [f1;f2;zeros(3,1)];

    x_p = x_p + f*dt;

    lat = x_p(1); 
    lon = x_p(2); 
    h = x_p(3);   
    vN = x_p(4);
    vE = x_p(5);
    vD = x_p(6);
    ba1 = x_p(7);
    ba2 = x_p(8);
    ba3 = x_p(9);

  A_pos = substitution(lat,vE,vN,vD,h);

  P_pos = P_pos + dt * (A_pos * P_pos + P_pos * A_pos' + Q_pos);

end
