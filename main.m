clc
clear all
close all

syms t real

format long

imuFs=100; %imu frequency
gpsFs=10; %gps frequency

localOrigin = [38.090826 13.302356 195.035000]; % Palermo

rng('default');

%% Create a circular trajectory
        
 % Trajectory parameters
    r = 20; % radius(m)
    t_max = 10; %seconds %simulation time
    campioni = imuFs*t_max; %samples
    
    angxsec = 2*pi/campioni;
    
    alfa = 0:angxsec:2*pi;

    %%Real local position
    x_pos = r .* cos(alfa); 
    y_pos = r .* sin(alfa);
    z_pos = zeros(size(x_pos));

    position = [x_pos' y_pos' z_pos'];

%% Moto circolare uniforme (MCU) (solo traslatorio)

        angvel_mcu_matr=zeros(imuFs*t_max,3); %angular velocity matrix init (MCU)                              
        acc_mcu_matr=zeros(imuFs*t_max,3); %%acceleration matrix init (MCU)
        trajectory = zeros(imuFs*t_max,3); %%trajectory of waypoints init

        t1 = 0:1/imuFs:t_max; %time vector
        
        idx=0; %index init
        
        %%NO ATTITUDE 
        att_eul_matr_mcu = zeros(imuFs*t_max+1,3); 
        orientation_mcu = quaternion(att_eul_matr_mcu, 'euler', 'ZYX', 'frame');

        toa = t1(1,:); %time of arrival

        % Generate trajectory
        groundTruth_mcu = waypointTrajectory('SampleRate', imuFs, ...
                                             'SamplesPerFrame',1, ...
                                             'Waypoints', position, ...
                                             'TimeOfArrival', toa, ...
                                             'Orientation', orientation_mcu);

        while ~isDone(groundTruth_mcu)
            %% Ricavo POSIZIONE, ORIENTAMENTO, VELOCITA, VEL_ANGOLARE, ACCELERAZIONE VERI
            [cerchio_waypoints,~,~,acc_mcu,angvel_mcu] = groundTruth_mcu();
            
            idx=idx+1;
       
            trajectory(idx,:)=cerchio_waypoints;    
            angvel_mcu_matr(idx,:)=angvel_mcu;                                                                                               
            acc_mcu_matr(idx,:)=acc_mcu;   

        end

%% Moto roto-traslatorio

        % Define orientation IN BODY
        yaw = deg2rad(15)*sin(2*pi*t+120*pi/180);
        yaw(t)=yaw;

        pitch = deg2rad(10)*sin(2*pi*t+240*pi/180);
        pitch(t)=pitch;

        roll = deg2rad(8)*sin(2*pi*t+360*pi/180); 
        roll(t)=roll;

        att_eul_matr = zeros(imuFs*t_max,3); %attitude euler angles matrix init

        t = 0:1/imuFs:t_max; 
        
        %%Euler angles
        for j=1:1:imuFs*t_max+1
        
                    psi = yaw(t(j)); %radians
                    theta = pitch(t(j));
                    phi  = roll(t(j));
            
            att_eul_matr(j,:)=[psi theta phi]; 
        
        end
        
        %%attitude in quaternion form
        orientation = quaternion([att_eul_matr(:,1), att_eul_matr(:,2), att_eul_matr(:,3)], 'euler', 'ZYX', 'frame');
    
        toa = t(1,:);

        % Generate trajectory
        groundTruth = waypointTrajectory('SampleRate', imuFs, ...
                                         'SamplesPerFrame',1, ...
                                         'Waypoints', position, ...
                                         'TimeOfArrival', toa, ...
                                         'Orientation', orientation);

        
        %REAL POSITION INIT
        posvera_matr = zeros(imuFs*t_max,3);                                  % NED                                 
        %REAL ATTITUDE INIT
        yprvera_matr=zeros(imuFs*t_max,3);
        %REAL VELOCITY INIT
        velvera_matr=zeros(imuFs*t_max,3);                                    % NED
        %REAL ANGULAR VELOCITY INIT
        angvelvera_matr=zeros(imuFs*t_max,3);                                 % NED                                
        %REAL ACCELERATION INIT
        accvera_matr=zeros(imuFs*t_max,3);                                    
        q_from_am_matr=zeros(imuFs*t_max,3);

        %IMU settings
        imu = imuSensor('accel-gyro-mag', 'SampleRate', imuFs); 

        imu.Gyroscope.ConstantBias = [0.02 0.04 0.08];
        imu.Gyroscope.NoiseDensity = [0.004 0.004 0.004];

        imu.Accelerometer.ConstantBias = [0.03 0.06 0.09];
        imu.Accelerometer.NoiseDensity = [0.004 0.004 0.004];

        %%GPS settings
        gps = gpsSensor('SampleRate',gpsFs,'PositionInputFormat','Local','ReferenceFrame', 'NED');  %NED  default
        gps.ReferenceLocation = localOrigin;
        gps.HorizontalPositionAccuracy = 1.0;
        gps.VerticalPositionAccuracy = 1.0;
        gps.VelocityAccuracy = 0.1;

        %%IMU & GPS measurements matrix init
        accelmis_matr = zeros(imuFs*t_max,3); % tangential acceleration
        gyromis_matr = zeros(imuFs*t_max,3);
        magmis_matr = zeros(imuFs*t_max,3);

        pos_gps = zeros(imuFs*t_max,3);
        vel_gps = zeros(imuFs*t_max,3);

        %%Total body acceleration matrix init
        accel_body_true_matr = zeros(imuFs*t_max,3);
        accel_body_meas_matr = zeros(imuFs*t_max,3);
        accel_ned_meas_matr=zeros(imuFs*t_max,3);

        %%Angular velocity matrix init
        angvel_body_true_matr = zeros(imuFs*t_max,3);
        angvel_body_meas_matr = zeros(imuFs*t_max,3);
        angvel_ned_meas_matr=zeros(imuFs*t_max,3);
        
        %%Estimated bias gyroscope matrix init
        bias_gyro_matr = zeros(imuFs*t_max,3);
        %%Estimated bias accelerometer matrix init
        bias_acc_matr = zeros(imuFs*t_max,3);

        %%Real geodetic position matrix
        real_geo_pos = zeros(imuFs*t_max,3);

        %%Estimated geodetic and local position matrix init
        est_geo_pos = zeros(imuFs*t_max,3);
        est_local_pos = zeros(imuFs*t_max,3);
        
        %%Estimated velocity matrix init
        vel_est = zeros(imuFs*t_max,3);

        %%Estimated attitude euler and quaternion matrix init
        eul_est_radians=zeros(imuFs*t_max,3);
        eul_est_degree=zeros(imuFs*t_max,3);
        est_quat_matr=zeros(imuFs*t_max,4);

        %%Q-TRIAD measurement attitude euler and quaternion matrix init
        triad_meas_eul=zeros(imuFs*t_max,3);
        triad_meas_quat=zeros(imuFs*t_max,4);

        dt=0.01; %integration step

%% ATTITUDE EKF init
% Get the initial ground truth pose from the first sample of the trajectory
% and release the ground truth trajectory to ensure the first sample is not 
% skipped during simulation.
     [initialPos, initialAtt, initialVel] = groundTruth();  
     reset(groundTruth);
    
        idx=0; %index init

        initstate = zeros(7,1);
        initstate(1:4) = compact(initialAtt)';  
        initstate(5:7) = 0;  
        x = initstate; %initial state vector for Attitude EKF               
       
        % Initialise covariance matrix
        cAtt0 = 1e-7;
        cBias0 = 1;
        P = diag([[1 1 1 1] * cAtt0, [1e-6 1e-6 1e-6] * cBias0]);
        
        % Process noise matrix
        nProcAtt  = 1e-17;
        nProcBias = 1;
        Q = diag([[1 1 1 1] * nProcAtt, [1e-9 1e-15 1e-15] * nProcBias]);
        
        % Measurement noise matrix
        R = diag([0.134276651403430 0.00200164425615 0.001021361257225 0.000708066823338].^2);
        
%% POSITION EKF init

        %%Initialise position for state vector  
        [lat_init,lon_init,h_init] = local2latlon(posvera_matr(1,1),posvera_matr(1,2),posvera_matr(1,3),localOrigin);

        initstate_pos = zeros(9,1);
        initstate_pos(1:3) = [lat_init lon_init h_init]';
        initstate_pos(4:6) = 0; 
        initstate_pos(7:9) = 0;
        x_p = initstate_pos; %initial state vector for Position EKF 

        R_pos = diag([0.000000000007321 0.000000000029513 0.277123557327953 (0.009725562107564^3)/3 0.009949335669347^3 0.009927574788694]);

        % Initialise covariance matrix
        nPos0 = 1e-9;
        nVel0 = 1e-2;
        nBiasAcc0 = 1;
        P_pos = diag([[1 1 1] * nPos0, [1 1 1] * nVel0, [1e-2 1e-2 1] * nBiasAcc0]);

        % Process noise matrix
        nProcPos = 1e-7;
        nProcVel = 1e-10;
        nProcBiasAcc = 1;
        Q_pos  = diag([[1 1 1] * nProcPos, [1 1 1] * nProcVel, [1e-7 1e-6 1e-7] * nProcBiasAcc]);

    contatore = 1; %counter init

    r_e = 6378137; % m %equatorial radius
    r_p = 6356752.3; %m %polar radius
    e = sqrt(1-(r_p/r_e)^2); % eccentricity

    %Parameters for gravity model wgs84
    gws0=9.7803267714;
    gws1=0.00193185138639;

    while ~isDone(groundTruth)

        %%Real position, attitude, velocity, acceleration, angular
        %%velocity from trajectory 
        [truePosition, trueOrientation, trueVel, trueAcc, trueAngVel] = groundTruth();

        idx=idx+1;
        
        %%IMU measurements
        [accel, gyro, mag] = imu(trueAcc-acc_mcu_matr(idx,:)+[0 0 2*9.8],trueAngVel-angvel_mcu_matr(idx,:),trueOrientation);

        accelmis_matr(idx,:) = accel;                                       %measured rotational acceleration - BODY 
        gyromis_matr(idx,:) = gyro;                                         %measured angular velocity in BODY
        magmis_matr(idx,:) = mag;                                           %measured magnetic field  

        %% Attitude estimation

            [x,P] = att_ekf(x,dt,gyro,accel,mag,P,R,Q); %Attitude EKF
                   
        %%GPS Measurements
        [lla, gpsVel] = gps(truePosition, trueVel);
        
        pos_gps(idx,:)=lla;                                                  %GPS geodetic position and altitude 
        vel_gps(idx,:)=gpsVel;                                               %GPS velocity

        %%Gravity model
        N = r_e/sqrt(1-(e^2)*(sin(x_p(1,1)))^2);                             %normal radius
        M = N*(1-e^2)/(1-(e^2)*(sin(x_p(1,1)))^2)^2;                         %meridian radius
        Radius = sqrt(M*N);

        gamma = 9.8*(Radius/(Radius+x_p(3,1)))^2;                            %gravity that change with altitude
        %gamma = gws0*((1+gws1*(sin(x_p(1,1)))^2)/sqrt(1-(e^2)*(sin(x_p(1,1)))^2))-2*gws0*x_p(3,1)/Radius;

        C_ned_b = [x(1)^2+x(2)^2-x(3)^2-x(4)^2 2*(x(2)*x(3)-x(1)*x(4)) 2*(x(2)*x(4)+x(1)*x(3));
                    2*(x(2)*x(3)+x(1)*x(4)) x(1)^2-x(2)^2+x(3)^2-x(4)^2 2*(x(3)*x(4)-x(1)*x(2));
                    2*(x(2)*x(4)-x(1)*x(3)) 2*(x(3)*x(4)+x(1)*x(2)) x(1)^2-x(2)^2-x(3)^2+x(4)^2]; %rotation matrix body to ned with estimated quaternion

        %%Total acceleration and angular velocity in BODY frame
        %trueAcc-[0 0 gamma] true acceleration NED
        accel_body = C_ned_b'*(trueAcc-[0 0 gamma])';                        %true total acceleration in BODY (rotated with estimated quaternion)
        angvel_body = C_ned_b'*trueAngVel';                                  %true angular velocity in BODY  

        [accel_body_meas,angvel_meas,~] = imu(-accel_body' + [0 0 9.8],angvel_body'); 

%         [accel_ned_meas,angvel_ned_meas] = imu(-(trueAcc-[0 0 gamma])+[0 0 9.8],trueAngVel);

        %% Prediction Position EKF - Prediction faster than update

            [x_p,P_pos] = predict_pos(x_p, dt, accel_body_meas, x(1), x(2), x(3), x(4), P_pos, Q_pos, N, M);
        
        contatore = contatore + 1;
    
       %% Update Position EKF
       if contatore == 5
    
             [x_p,P_pos] = update_pos(x_p, gpsVel, lla, P_pos, R_pos); 
    
             contatore = 0;
       end   

        posvera_matr(idx,:)=truePosition;                                    %real position NED 

        accel_body_true_matr(idx,:) = accel_body;                            %true total acceleration in BODY
        accel_body_meas_matr(idx,:) = accel_body_meas;                       %measured total acceleration in BODY
%         accel_ned_meas_matr(idx,:) = accel_ned_meas;                         %measured total acceleration in NED

        angvel_body_true_matr(idx,:) = angvel_body;                          %true angular velocity in BODY
        angvel_body_meas_matr(idx,:) = angvel_meas;                          %measured angular velocity in BODY
%         angvel_ned_meas_matr(idx,:) = angvel_ned_meas;                       %measured angular velocity in NED
        
        trueOrientation_vett = compact(trueOrientation);                   
        ypr_true = quat2eul(trueOrientation,'ZYX');                                                                                                                                                                    
        yprvera_matr(idx,:)=rad2deg(ypr_true);                               %true attitude in euler angles in degrees                                     
         
        velvera_matr(idx,:)=trueVel;                                         %true velocity NED
        angvelvera_matr(idx,:)=trueAngVel;                                   %true angular velocity NED 
        accvera_matr(idx,:)=trueAcc;                                         %true total acceleration NED

        est_quat_matr(idx,:) = [x(1) x(2) x(3) x(4)];                        %estimated attitude in quaternion form
        est_eul_angles = quat2eul([x(1) x(2) x(3) x(4)],'ZYX');

        eul_est_radians(idx,:) = est_eul_angles;                             %estimated attitude in euler angles in radians
        eul_est_degree(idx,:) = rad2deg(est_eul_angles); %degree             %estimated attitude in euler angles in degrees

        psi_est = eul_est_radians(idx,1);                                    %estimated yaw, pitch, roll in radians
        theta_est = eul_est_radians(idx,2);
        phi_est = eul_est_radians(idx,3);

        bias_gyro_hat = [x(5) x(6) x(7)];                   
        bias_gyro_matr(idx,:)=bias_gyro_hat;                                 %estimated bias gyroscope
        
        q_from_am=qtriad(accel',mag');                                       %attitude measurements in quaternion form with Q-TRIAD
        triad_meas_quat(idx,:) = q_from_am;
        triad_att=quat2eul(q_from_am,'ZYX');
        triad_meas_eul(idx,:)=triad_att;                                     %attitude measurements in euler angles with Q-TRIAD
        
        est_geo_pos(idx,:) = [x_p(1,1) x_p(2,1) x_p(3,1)];                   %estimated geodetic position
        vel_est(idx,:) = [x_p(4,1) x_p(5,1) x_p(6,1)];                       %estimated velocity
        bias_acc_matr(idx,:) = [x_p(7,1) x_p(8,1) x_p(9,1)];                 %estimated bias accelerometer
        
        [real_lat,real_lon,real_h] = local2latlon(posvera_matr(idx,1),posvera_matr(idx,2),posvera_matr(idx,3),localOrigin);
        real_geo_pos(idx,:) = [real_lat real_lon real_h];                    %true geodetic position
        
        [est_x, est_y, est_z] = latlon2local(x_p(1,1),x_p(2,1),x_p(3,1),localOrigin);
        est_local_pos(idx,:) = [est_x est_y est_z];                          %estimated local position

        err_att_eul_degrees = eul_est_degree-yprvera_matr;                   %attitude error in euler angles in degrees
        err_vel = vel_est-velvera_matr;                                      %velocity error 
        err_pos_local = est_local_pos - posvera_matr;                        %local position error
        err_bias_gyro = bias_gyro_matr - ones(imuFs*t_max,3).*[imu.Gyroscope.ConstantBias];   %error bias gyro
        err_bias_acc = bias_acc_matr - ones(imuFs*t_max,3).*[imu.Accelerometer.ConstantBias]; %error bias acc
   
    end

%% PLOT
 
% Real and estimated attitude
figure        
plot(eul_est_degree) % in degrees
hold on
plot(rad2deg(att_eul_matr))  %in degrees
legend('yaw_s','pitch_s','roll_s','yaw_v','pitch_v','roll_v')
xlabel('samples')
ylabel('degrees')
title('True and estimated attitude (euler angles in degrees)')

% Attitude measurements in euler angles with Q-TRIAD
figure
plot(triad_meas_eul*180/pi) %in degrees
legend('yaw triad','pitch triad','roll triad')
xlabel('samples')
ylabel('degrees')
title('Attitude measurements with Q-Triad')

% Real and estimated bias gyro
figure
plot(bias_gyro_matr)
hold on 
plot(ones(imuFs*t_max,3).*[imu.Gyroscope.ConstantBias])
legend('bwx_s','bwy_s','bwz_s','bwx_v','bwy_v','bwz_v')
xlabel('samples')
title('Estimated and true bias gyroscope')

% Real and estimated geodetic position
figure 
geoplot(est_geo_pos(:,1),est_geo_pos(:,2))
hold on
geoplot(real_geo_pos(:,1),real_geo_pos(:,2))
legend('est geo pos','real geo pos')
title('True and estimated geodetic position')

% Real and estimated local position
figure 
plot(est_local_pos)
hold on
plot(posvera_matr)
legend('posx_s','posy_s','posz_s','posx_v','posy_v','posz_v')
xlabel('samples')
ylabel('meters')
title('True and estimated local position')

% Estimated and real altitude
figure
plot(est_geo_pos(:,3))
hold on
plot(real_geo_pos(:,3))
legend('est altitude','real altitude')
xlabel('samples')
ylabel('meters')
title('Estimated and real altitude')

% Estimated, true and measured velocity in NED frame
figure
plot(vel_est)
hold on
plot(velvera_matr)
plot(vel_gps)
legend('velx_s','vely_s','velz_s','velx_v','vely_v','velz_v','velx_m','vely_m','velz_m')
xlabel('samples')
ylabel('m/s')
title('Estimated, true and measured velocity in NED frame')

%Estimated and real bias accelerometer
figure
subplot(3,1,1)
plot(bias_acc_matr(:,1),'b')
hold on
plot(ones(imuFs*t_max,3).*imu.Accelerometer.ConstantBias(1,1),'r')
legend('bax_s','bax_v','Location','southeast')
xlabel('samples')
title('Estimated and true bias x accelerometer')
subplot(3,1,2)
plot(bias_acc_matr(:,2),'m')
hold on
plot(ones(imuFs*t_max,3).*imu.Accelerometer.ConstantBias(1,2),'c')
legend('bay_s','bay_v','Location','southeast')
xlabel('samples')
title('Estimated and true bias y accelerometer')
subplot(3,1,3)
plot(bias_acc_matr(:,3),'g')
hold on
plot(ones(imuFs*t_max,3).*imu.Accelerometer.ConstantBias(1,3),'k')
legend('baz_s','baz_v','Location','southeast')
xlabel('samples')
title('Estimated and true bias z accelerometer')
hold off

%%Errors

%Attitude errors in euler angles in degrees
%%Velocity errors
figure
subplot(2,1,1)
plot(err_att_eul_degrees)
legend('yaw error','pitch error','roll error')
xlabel('samples')
ylabel('degrees')
title('Attitude errors in euler angles (degrees)')
subplot(2,1,2)
plot(err_vel)
legend('vel x err','vel y err','vel z err')
xlabel('samples')
ylabel('m/s')
title('Velocity errors')

%Bias gyroscope errors
figure 
plot(err_bias_gyro)
legend('bias gyro x err','bias gyro y err','bias gyro z err','Location','southeast')
xlabel('samples')
title('Bias gyroscope errors')

%Bias accelerometer errors
figure
subplot(3,1,1)
plot(err_bias_acc(:,1))
legend('bias acc x err')
xlabel('samples')
title('Bias accelerometer x errors')
subplot(3,1,2)
plot(err_bias_acc(:,2))
legend('bias acc y err')
xlabel('samples')
title('Bias accelerometer y errors')
subplot(3,1,3)
plot(err_bias_acc(:,3))
legend('bias acc z err')
xlabel('samples')
title('Bias accelerometer z errors')

%Real and measured acceleration
figure
plot(accel_body_true_matr)
hold on
plot(accel_body_meas_matr)  
legend('acc x true','acc y true','acc z true','acc x est','acc y est','acc z est')
xlabel('samples')
ylabel('m/s^2')
title('True and measured acceleration in BODY frame')

%Real and measured angular velocity
figure
plot(angvel_body_true_matr)
hold on
plot(angvel_body_meas_matr)
legend('ang vel x true','ang vel y true','ang vel z true','ang vel x meas','ang vel y meas','ang vel z meas')
xlabel('samples')
ylabel('rad/s')
title('True and measured angular velocity in BODY frame')

% figure
% plot(bias_acc_matr)
% hold on 
% plot(ones(imuFs*t_max,3).*[imu.Accelerometer.ConstantBias])
% legend('bax_s','bay_s','baz_s','bax_v','bay_v','baz_v','Location','southeast')
% title('bias acc stimato e vero')
% hold off

%%PLOT x relazione

%%True data (vel NED, acc NED, ang vel NED)
% figure
% subplot(3,1,1)
% plot(velvera_matr)
% legend('velx_v','vely_v','velz_v')
% xlabel('samples')
% ylabel('m/s')
% title('True velocity in NED frame')
% subplot(3,1,2)
% plot(accvera_matr)
% legend('acc_x_v','acc_y_v','acc_z_v')
% xlabel('samples')
% ylabel('m/s^2')
% title('True acceleration in NED frame')
% subplot(3,1,3)
% plot(angvelvera_matr)
% legend('angvel_x_v','angvel_y_v','angvel_z_v')
% xlabel('samples')
% ylabel('rad/s')
% title('True angular velocity in NED frame')
% 
%%IMU measuremets (acc NED, ang vel NED, magnetic field)
% figure
% subplot(3,1,1)
% plot(accel_ned_meas_matr)
% hold on
% plot(accvera_matr)
% legend('acc_x meas','acc_y meas','acc_z meas','acc_x true','acc_y true','acc_z true')
% xlabel('samples')
% ylabel('m/s^2')
% title('Measured and true acceleration in NED frame')
% subplot(3,1,2)
% plot(angvel_ned_meas_matr)
% hold on
% plot(angvelvera_matr)
% legend('angvel_x meas','angvel_y meas','angvel_z meas','angvel_x true','angvel_y true','angvel_z true')
% xlabel('samples')
% ylabel('rad/s')
% title('Measured and true angular velocity in NED frame')
% subplot(3,1,3)
% plot(magmis_matr)
% legend('mag_x','mag_y','mag_z')
% xlabel('samples')
% ylabel('μT')
% title('Magnetic field')
% 
% %%GPS measurements
% geoplot(pos_gps(:,1),pos_gps(:,2))
% title('Measured GPS position')
% 
% 
% figure
% plot(pos_gps(:,3))
% legend('altitude')
% xlabel('samples')
% ylabel('meters')
% title('gps altitude')
% 
% figure
% plot(vel_gps)
% legend('vel gps_x','vel gps_y','vel gps_z')
% xlabel('samples')
% ylabel('m/s')
% title('Measured GPS velocity in NED frame')
% 
% 
% 
% 
% 
%             