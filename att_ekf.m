function [x,P] = att_est(x, dt, gyro, accel, mag, P, R, Q)
   
    w1 = gyro(1);
    w2 = gyro(2);
    w3 = gyro(3);
    
    q0 = x(1);
    q1 = x(2);
    q2 = x(3);
    q3 = x(4);
    bw1 = x(5);
    bw2 = x(6);
    bw3 = x(7);

    ax = accel(1);
    ay = accel(2);
    az = accel(3);

    mx = mag(1);
    my = mag(2);
    mz = mag(3); 

 % State transition function, xdot = f(x, u)
        C_body2qdot = 0.5 * [-q1, -q2, -q3; ...
                              q0, -q3,  q2; ...
                              q3,  q0, -q1; ...
                             -q2,  q1,  q0];
                         
        f = [ C_body2qdot * [w1 - bw1; w2 - bw2; w3 - bw3]; 0; 0; 0];

        
        x = x + f*dt;

            q0 = x(1);
            q1 = x(2);
            q2 = x(3);
            q3 = x(4);
            bw1 = x(5);
            bw2 = x(6);
            bw3 = x(7);

     A = 0.5 * [0 -(w1-bw1) -(w2-bw2) -(w3-bw3) q1 q2 q3;
               (w1-bw1) 0 (w3-bw3) -(w2-bw2)  -q0 q3 -q2;
               (w2-bw2) -(w3-bw3) 0 (w1-bw1)  -q3 -q0 q1;
               (w3-bw3) (w2-bw2) -(w1-bw1) 0  q2 -q1 -q0;
                0 0 0 0 0 0 0;
                0 0 0 0 0 0 0;
                0 0 0 0 0 0 0];
                
         P = P + dt * (A * P + P * A' + Q);

         H = [eye(4) zeros(4,3)];   

         K = P * H' * inv(H * P * H' + R);
         P = (eye(7) - K*H)*P;

         q_from_am = qtriad(accel',mag'); 

         triad_eul = quat2eul(q_from_am,'ZYX');
         Hx_eul=quat2eul((H*x)','ZYX');
         
         %DIFFERENZA ANGOLI DI EULERO
         e_att = triad_eul - Hx_eul;

         err_att_quat=eul2quat(e_att,'ZYX');

         x = x + K*(err_att_quat');
        
            qNorm = sqrt(x(1) * x(1) + x(2) * x(2) + x(3) * x(3) + x(4) * x(4));
            x(1) = x(1) / qNorm;
            x(2) = x(2) / qNorm;
            x(3) = x(3) / qNorm;
            x(4) = x(4) / qNorm;


end
