%%% GATTO Model
function [Xdot] = gatto(tt,X)

global delta_E delta_P sigma eta gamma_I gamma_A gamma_Q gamma_H alpha_I alpha_H xi Lvect beta_P beta_I beta_A theta 

S = X(1);
E = X(2);
P = X(3);
I = X(4);
A = X(5);
H = X(6);
Q = X(7);
R = X(8);
D = X(9);


L = Lvect(fix(tt)+1);

Xdot = zeros(9,1);
 
    Xdot(1) = -((beta_P*P+beta_A*A+beta_I*I)/(S+E+P+I+A+R))*S*(1-theta*L)^2;
    Xdot(2) = ((beta_P*P+beta_A*A+beta_I*I)/(S+E+P+I+A+R))*S*(1-theta*L)^2 - delta_E*E;
    Xdot(3) = delta_E*E - delta_P*P;
    Xdot(4) = sigma*delta_P*P - (eta+gamma_I+alpha_I)*I;
    Xdot(5) = (1-sigma)*delta_P*P - gamma_A*A;
    Xdot(6) = (1-xi)*eta*I - (gamma_H + alpha_H)*H;
    Xdot(7) = xi*eta*I - gamma_Q*Q;
    Xdot(8) = gamma_I*I + gamma_A*A + gamma_H*H+gamma_Q*Q;
    Xdot(9) = alpha_I*I + alpha_H*H;