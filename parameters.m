%% GATTO parameters
 
syms beta_P beta_I beta_A

delta_E = 1/3.32; %rate at which the EXPOSED people become PRESYMPTOMATIC
delta_P = 1/0.75; % rate at which PRESYMPTOMATICS become SERIOUSLY INFECTED
sigma = 0.25; % probability with which PRESYMPTOMATICS become SERIOUS INFECTED
eta = 1/4.05; % rate at which SERIOUS INFECTED are isolated
xi = 0.4; % fraction of SERIOUS INFECTIONS going to quarantine

%% In order to reduce the number of parameters to be estimated it is necessary: 
 
    % gamma_H=gamma_Q=gamma_I 
    % gamma_A = 2*gamma_I
    % alfa_I=alfa_H
    
gamma_I  = 1/14.32; % rate at which SERIOUS INFECTIONS recover from infection
gamma_H = 1/14.32; % rate at which HOSPITALIZEDs recover from infection
gamma_Q =  1/14.32; % rate at which those in QUARANTINE recover from infection
gamma_A = 2*gamma_I; % rate at which ASYMPTOMATICS recover from infection
 
alpha_I = 1/24.23 ; % rate at which SERIOUS INFECTED dies
alpha_H = 1/24.23 ; % rate at which HOSPITALIZEDs die
 
%% INITIAL CONDITIONS

X0 = [1-0.0000000833 0.0000000833 0 0 0 0 0 0 0];  %S E P I A H Q R D 
 
%% BETA PARAMETERS

R0_P = beta_P/delta_P;
R0_I = sigma*(beta_I/(eta+alpha_I+gamma_I));
R0_A = (1-sigma)*(beta_A/gamma_A);

 
eqns = [beta_A - 0.036 * beta_P == 0, beta_I - 1.38 * beta_A == 0, R0_P + R0_I + R0_A - R0 == 0];
S = solve(eqns, [beta_P beta_I beta_A ]);
        
beta_A = double(S.beta_A); % syms to double
beta_I = double(S.beta_I);
beta_P = double(S.beta_P);  

%% Cost Function Parameters
 
PIL_2019_Italy = 1787664000000; % source: istat
italian_pop = 59641488;
PIL_ITA_pc = (PIL_2019_Italy/italian_pop)/365;  % Italian per capita GDP 
 
VSL = 20*PIL_ITA_pc; % value of statistical life
 
u = 0.2; % weight between economic losses and life losses
theta = 0.8; % effectiveness of the lockdown