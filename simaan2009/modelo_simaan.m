% glauberrleite
% Modelagem do sistema cardiovascular humano 2019.1
% Universidade Federal de Alagoas

clc;
clear all;

%% Model Parameters
R_s = 1.0000; % Systemic Vascular Resistance (SVR)
R_m = 0.0050; % Mitral Valve Resistance
R_a = 0.0010; % Aortic Valve Resistance
R_c = 0.0398; % Characteristic Resistance

C_ae = 4.4000; % Left Atrial Compliance
C_s = 1.3300; % Systemic Compliance
C_ao = 0.0800; % Aortic Compliance

L_s = 0.0005; % Inertance of blood in Aorta
V_0 = 10; % Reference volume (mL)

E_max = 2;
E_min = 0.05;
HeartRate = 75;
tc = 60/HeartRate; % Cardiac cycle interval
t_max = 0.2 + 0.15*tc;

%% Model functions and state space matrices
E_n = @(t_n) 1.55 ...
    * (((t_n/0.7)^1.9)/(1 + (t_n/0.7)^1.9)) ...
    * 1/(1 + (t_n/1.17)^21.9);

E = @(t) (E_max - E_min) * E_n(mod(t, tc)/t_max) ...
    + E_min;

P_ve = @(t, V_ve) E(t)*(V_ve - V_0);

% x = [V_ve, P_ae, Q_a, P_ao, P_s]
A = @(t, D_m, D_a) [-(D_m/R_m + D_a/R_a)*E(t), D_m/R_m, 0, D_a/R_a, 0 ;
    (D_m * E(t))/(R_m * C_ae), -(C_ae^-1)*(1/R_s + D_m/R_m), 0, 0, (R_s*C_ae)^-1;
    0, 0, -R_c/L_s, L_s^-1, -L_s^-1;
    (D_a * E(t))/(R_a*C_ao), 0, -C_ao^-1, -D_a/(R_a*C_ao), 0;
    0, (R_s*C_s)^-1, C_s^-1, 0, -(R_s*C_s)^-1];

p = @(t, D_m, D_a) [(D_m/R_m + D_a/R_a)*E(t)*V_0;
    -(D_m*E(t)*V_0)/(R_m*C_ae);
    0;
    -(D_a*E(t)*V_0)/(R_a*C_ao);
    0];

%% Preparing simulation
T_s = 0.0001;
t = 0;
t_end = 3;

data.x.V_ve = zeros(t_end/T_s, 1);
data.x.P_ae = zeros(t_end/T_s, 1);
data.x.Q_a = zeros(t_end/T_s, 1);
data.x.P_ao = zeros(t_end/T_s, 1);
data.x.P_s = zeros(t_end/T_s, 1);
data.P_ve = zeros(t_end/T_s, 1);
data.t = t:T_s:t_end;

%x = @(i) [data.x.V_ve(i); data.x.P_ae(i); data.x.Q_a(i); data.x.P_ao(i); data.x.P_s(i)];

% Initial conditions
data.x.V_ve(1) = 140;
data.x.P_ae(1) = 5;
data.x.Q_a(1) = 0;
data.x.P_ao(1) = 90;
data.x.P_s(1) = 90;
data.P_ve(1) = P_ve(t, data.x.V_ve(1));

D_a = 1;
D_m = 1;

%% Running simulation
n = 1;
while t < t_end
    t = t + T_s;
    
    % Diod values
    if (data.x.P_ae(n) > data.P_ve(n))
        D_m = 1;
        D_a = 0;
    elseif (data.P_ve(n) > data.x.P_ao(n))
        D_m = 0;
        D_a = 1;
    else
        D_m = 0;
        D_a = 0;
    end
    
    x = [data.x.V_ve(n); data.x.P_ae(n); data.x.Q_a(n); data.x.P_ao(n); data.x.P_s(n)];
    
    % RK4
    dx = A(t, D_m, D_a) * x +  p(t, D_m, D_a);
    kx1 = T_s * dx;
    x1 = x + 0.5*kx1;
    
    dx = A(t, D_m, D_a) * x1 +  p(t, D_m, D_a);
    kx2 = T_s * dx;
    x1 = x + 0.5*kx2;
    
    dx = A(t, D_m, D_a) * x1 +  p(t, D_m, D_a);
    kx3 = T_s * dx;
    x1 = x + kx3;
    
    dx = A(t, D_m, D_a) * x1 +  p(t, D_m, D_a);
    kx4 = T_s * dx;
    
    xf = x + (kx1 + 2 * kx2 + 2 * kx3 + kx4)/6;
    
    % Saving data
    data.x.V_ve(n+1) = xf(1);
    data.x.P_ae(n+1) = xf(2);
    data.x.Q_a(n+1) = xf(3);
    data.x.P_ao(n+1) = xf(4);
    data.x.P_s(n+1) = xf(5);
    
    data.P_ve(n+1) = P_ve(t, data.x.V_ve(n+1));
    
    n = n + 1;
end

%% Ploting results
figure(1);
subplot(3,1,1);
plot(data.t, data.x.P_ao, 'k', data.t, data.P_ve, ':k', data.t, data.x.P_ae, '--k');
title('Results - Glauber R. Leite');
legend('AoP', 'LVP', 'LAP');
ylabel('Pressures (mmHg)');
axis([0, 3, -5, 150]);
grid on;

subplot(3,1,2);
plot(data.t, data.x.Q_a, 'k');
ylabel('Aortic Flow (ml/s)');
axis([0, 3, -50, 700]);
grid on;

subplot(3,1,3);
plot(data.t, data.x.V_ve, 'k');
ylabel('LVV (ml)');
xlabel('time (s)');
axis([0, 3, 50, 160]);
grid on;