%% MAIN 

clc; clear all;

%% Global

    global delta_E delta_P sigma eta gamma_I gamma_A gamma_Q gamma_H alpha_I alpha_H...
        xi X0 beta_P beta_I beta_A N Lvect PIL_ITA_pc VSL theta u 

    R0 = 3.6; %basic reproduction number from Gatto paper

%% Simulation of Gatto Model for 200 days

    parameters; % load parameters

    N=200;
    Lvect = zeros(1,N);
    time = 0:1:N-1; 

    Xnolk = ode3('gatto',time,X0);

    figure('Name', 'FREE RESPONSE GATTO')
    hold on; grid on;
    plot([0:N-1],Xnolk(:,1).*100,'k');
    plot([0:N-1],Xnolk(:,2).*100,'b');
    plot([0:N-1],Xnolk(:,3).*100,'g'); 
    plot([0:N-1],Xnolk(:,4).*100,'r'); 
    plot([0:N-1],Xnolk(:,5).*100,'c');  
    plot([0:N-1],Xnolk(:,6).*100,'m');  
    plot([0:N-1],Xnolk(:,7).*100,'y');  
    plot([0:N-1],Xnolk(:,8).*100,'o');  
    plot([0:N-1],Xnolk(:,9).*100,'-');  
    legend({'S' 'E' 'P' 'I' 'A' 'H' 'Q' 'R' 'D'});
    title('No lockdown case'); xlabel('Time (days)'); ylabel('Population (%)');
    
%% Theoric R0 evaluation

    lambda = (beta_P*Xnolk(5,3)+beta_A*Xnolk(5,5)+beta_I*Xnolk(5,4))/(sum(Xnolk(5,1:5))+Xnolk(5,8));

    R0_t = lambda*Xnolk(5,1)/(delta_E*Xnolk(5,2));

%% Operative R0 evaluation

    operative_R0; % computation of R0 which satisfies the doubling time condition 
                  % between 2.5 days-3.5 days
    R0 = R(end); 
    
    parameters;   %parameters evaluation with new R0

%% Optimal Lockdown with Simulated Annealing

    addpath('Optimal Lockdown')

% lockdown strategy with simulated annealing using a constant U0
    lockdown_sa 

% comparison between lockdown strategy evaluated for different u (cost function weight)
    open('cost_fun_weight_variation.fig')
    
%% GLOBAL SEARCH using Multistart Algorithm

% Searching an optimal U0
    
    %     options = optimoptions('fmincon','Display','iter-detailed','Algorithm','active-set','FunValCheck','on','MaxFunctionEvaluations',4000,'MaxIterations',300);
    %      
    %     %local minima MultiStart 
    %     rng default % For reproducibility
    %     problem = createOptimProblem('fmincon','objective',@cost_fun_sa,'x0',U0,'lb',lb,'ub',ub,'options',options);
    %     ms = MultiStart('MaxTime',1500,'StartPointsToRun','bounds','Display','iter','FunctionTolerance',1e-4);
    %     [U0_ms,fval] = run(ms,problem,10);
   
    load('U0_ms.mat', 'Uvec'); %loading optimal U0
    U0_ms=Uvec;

    optimal_lockdown

%% CYCLIC LOCKDOWN STRATEGY

    addpath('Cyclic Lockdown')

%%%%CYCLICAL STRATEGY
% stepped lockdown strategy where the S.A. finds the optimal step widths every 14 days

    step_function_strategy; 

%% MPC 

% Evaluation of the most important parameters using Sensitivity Analisis
    
    addpath('Simulink')
    
    %par_sim; %load simulink parameters
    %open gattosim;
    
    open('scatter_plot.fig')
    open('corr_lin_nonlin.fig')
    open('regr_lin_nonlin.fig')

 % MPC Strategy

    addpath('MPC')
     
    %MPC_main
    load('MPC_data.mat'); %loading optimal U0

    plot_mpc

%% Add realism with ODE series

    addpath('ODE series')

% Comparison Gatto Model with ODE series and Gatto Model 

    free_response_ODE

%Lockdown Strategy with simulated annealing using constant U0 and a decreasing over time U0


    %%% CONSTANT U0
    
    const_lockdown_ODE 

    %%% DECREASING U0
    
    addpath('Optimal Lockdown')
    
    load('U0_ms.mat', 'Uvec'); %loading optimal U0
    U0_ms=Uvec;
    
    addpath('ODE series')
    
    decreasing_lockdown_ODE

% Cyclic Lockdown Strategy using stepped strategy

    cyclic_lockdown_ODE
