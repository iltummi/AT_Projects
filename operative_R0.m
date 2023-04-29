%% Gatto Simulated Annealing

clear 
close all
clc

%% load parameters
R = 0;
global delta_E delta_P sigma eta gamma_I gamma_A gamma_Q gamma_H alpha_I alpha_H ... 
    xi X0 beta_P beta_I beta_A N Lvect PIL_ITA_pc VSL theta u 

 for R0 = 2:0.1:4

            parameters;

            %% Model simulation for 200 days
            N=200;
            Lvect = zeros(1,N);
            time = 0:1:N-1; 

            Xnolk = ode3('gatto',time,X0);

            %% doubling time evaluation
            step = 2; %time interval to evaluate m
            m = (log(Xnolk(17, 2))-log(Xnolk(15, 2)))/step;
            tdr = log(2)/m;
            
            if tdr>=2.5 && tdr<=3.5
                R = [R, R0];
            end
 end
 
 %semilogarithmic plot
 
 time2 = 10:1:30;
 figure
 semilogy(time2,Xnolk(time2,2:7).*100);grid on;
 legend({'E' 'P' 'I' 'A' 'H' 'Q'});