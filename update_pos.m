function [x_p,P_pos] = update_pos(x_p, gpsVel, lla, P_pos, R_pos)
    
    z = [lla';gpsVel'];

    H = [eye(6) zeros(6,3)]; 

    K = P_pos * H' * inv(H * P_pos * H' + R_pos);
        
    P_pos = (eye(9) - K*H)*P_pos;

    err = z - H*x_p;

    x_p = x_p + K*err;
    
end