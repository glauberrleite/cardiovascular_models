% glauberrleite
% Modelagem do sistema cardiovascular humano 2019.1
% Universidade Federal de Alagoas
% Projeto 3 - Simaan + ECG + LVAD pulsátil (PVAD)

clc;
clear;

%% Model and simulation parameters
R_s = 0.8738; % Systemic Vascular Resistance (SVR)
R_m = 0.0050; % Mitral Valve Resistance
R_a = 0.010; % Aortic Valve Resistance
R_c = 0.0398; % Characteristic Resistance

C_ae = 4.0000; % Left Atrial Compliance
C_s = 2.896; % Systemic Compliance
C_ao = 0.0800; % Aortic Compliance

L_s = 0.001025; % Inertance of blood in Aorta
V_0 = 15; % Reference volume (mL)

R_i = 0.15; % Inlet Resistance of Cannulae
R_o = @(Q_o) 0.05 + 0.00015 * abs(Q_o); % Outlet Resistance of Cannulae
L_i = 0.0854; % Inlet Inertance of Cannulae
L_o = 0.0087; % Outlet Inertance of Cannulae
L_p = 0.0033; % Inertance of Thoratec® Bladder
R_p = 0.05; % Resistance of Thoratec® Bladder
C_d = 4; % Compliance of pneumatic drive line
R_d = 0.01; % Resistance of pneumatic drive line
alpha_PVAD = 0.15; % PVAD parameter (s/mL)
C_p = 2; % Compliance of Thoratec® Bladder
Vd_vad = 107; % Reference Volume

P_e = 180; % Ejection pressure
P_f = 0; % Filling pressure

E_max = .8;
E_min = 0.05;
HeartRate = 60;
tc = 60/HeartRate; % Cardiac cycle interval
t_max = 0.2 + 0.1555*tc;

ejection_time = tc / 2; % Pe will remain constant for this fixed period
ejection_delay = 0.25 * tc; % The amount of delay after E(t) is triggered to start ejection

T_s = 0.0001;
t = 0;
t_end = 3;
N = t_end/T_s;
data.t = t:T_s:t_end;

%% ECG simulation applying MaMeMi Filter
data.x.xecg = zeros(N, 1);
data.x.yecg = zeros(N, 1);
data.x.zecg = zeros(N, 1);
maxi        = zeros(1,N);
mini        = zeros(1,N);
h           = zeros(1,N);
a           = zeros(1,N);
n           = zeros(1,N);
g           = zeros(1,N);
r           = zeros(1,N);
v           = zeros(1,N);
w           = zeros(1,N);
R_dtct      = zeros(1,N);
z_0         = zeros(1,N);
deltaS      = zeros(1,N);
deltaPVAD   = zeros(1,N);
gammad      = zeros(1,N);


data.x.xecg(1) = -1;

deltaS(1) = 1;

ip = 0;
lg = 0;
e  = zeros(1,N);

% R detector variables
sigma = 2;
delta = 2;
beta = 15;
atv = 0;

A_ecg = [1.2, -5, 30, -7.5, 0.75];
B_ecg = [0.25, 0.1, 0.1, 0.1, 0.4];
TH = [-1/3, -1/12, 0, 1/12, 1/2]*pi;

for i = 1:N-1
    % RK4
    k = 0;
    w_ecg = 2*pi/tc;
    alpha = 1-sqrt(data.x.xecg(i)^2 + data.x.yecg(i)^2);
    th = atan2(data.x.yecg(i),data.x.xecg(i));
    for ii = 1:5
       k = k - (A_ecg(ii)*(th-TH(ii))*exp(-(th-TH(ii))^2/(2*B_ecg(ii)^2)));
    end

    A_i = [alpha, -w_ecg, 0; w_ecg, alpha, 0; 0, 0, -1];
    B_i = [0; 0; k + z_0(i)];

    dx = A_i * [data.x.xecg(i); data.x.yecg(i); data.x.zecg(i)] +  B_i;
    kx1 = T_s * dx;
    x1 = [data.x.xecg(i); data.x.yecg(i); data.x.zecg(i)] + 0.5 * kx1;
       
    k = 0;
    alpha = 1-sqrt(x1(1)^2 + x1(2)^2);
    th = atan2(x1(2),x1(1));
    for ii = 1:5
       k = k - (A_ecg(ii)*(th-TH(ii))*exp(-(th-TH(ii))^2/(2*B_ecg(ii)^2)));
    end

    A_i = [alpha, -w_ecg, 0; w_ecg, alpha, 0; 0, 0, -1];
    B_i = [0; 0; k + z_0(i)];
    
    dx = A_i * x1 +  B_i;
    kx2 = T_s * dx;
    x1 = [data.x.xecg(i); data.x.yecg(i); data.x.zecg(i)] + 0.5 * kx2;
    
    k = 0;
    alpha = 1-sqrt(x1(1)^2 + x1(2)^2);
    th = atan2(x1(2),x1(1));
    for ii = 1:5
       k = k - (A_ecg(ii)*(th-TH(ii))*exp(-(th-TH(ii))^2/(2*B_ecg(ii)^2)));
    end

    A_i = [alpha, -w_ecg, 0; w_ecg, alpha, 0; 0, 0, -1];
    B_i = [0; 0; k + z_0(i)];
    
    dx = A_i * x1 +  B_i;
    kx3 = T_s * dx;
    x1 = [data.x.xecg(i); data.x.yecg(i); data.x.zecg(i)] + kx3;
    
    k = 0;
    alpha = 1-sqrt(x1(1)^2 + x1(2)^2);
    th = atan2(x1(2),x1(1));
    for ii = 1:5
       k = k - (A_ecg(ii)*(th-TH(ii))*exp(-(th-TH(ii))^2/(2*B_ecg(ii)^2)));
    end

    A_i = [alpha, -w_ecg, 0; w_ecg, alpha, 0; 0, 0, -1];
    B_i = [0; 0; k + z_0(i)];
    
    dx = A_i * x1 +  B_i;
    kx4 = T_s * dx;
    
    xf = [data.x.xecg(i); data.x.yecg(i); data.x.zecg(i)] + (kx1 + 2 * kx2 + 2 * kx3 + kx4)/6;
 
    % Updating values
    z_0(i+1) = 0.15*sin(2*pi*(60/(12+randn))*data.t(i+1));
    HeartRate = 60 + 2*randn;
    tc = 60/HeartRate;

    data.x.xecg(i+1) = xf(1);
    data.x.yecg(i+1) = xf(2);
    data.x.zecg(i+1) = xf(3);

    % R-wave detection
    if data.x.zecg(i+1) > maxi(i) 
        maxi(i+1) = maxi(i) + sigma*delta;    
    elseif data.x.zecg(i+1) <= maxi(i)
        maxi(i+1) = maxi(i) - delta;
    end

    if data.x.zecg(i+1) < mini(i) 
        mini(i+1) = mini(i) - sigma*delta;    
    elseif data.x.zecg(i+1) >= mini(i)
        mini(i+1) = mini(i) + delta;
    end

    h(i+1) = data.x.zecg(i+1) - (maxi(i+1)+mini(i+1))/2;

    a(i+1) = maxi(i+1)-mini(i+1);

    if a(i+1) <= h(i+1)
        n(i+1) = sign(h(i+1)*(abs(h(i+1))-a(i+1)));
    else
        n(i+1) = 0;
    end

    if i > beta
        if (n(i)>0) && (n(i)>n(i-beta)) && (n(i)>n(i+beta))   
            g(i) = n(i) - max(n(i-beta),n(i+beta));
        elseif (n(i)<0) && (n(i)<n(i-beta)) && (n(i)<n(i+beta))   
            g(i) = n(i) + min(n(i-beta),n(i+beta));
        else
            g(i) = 0;
        end

        if(g(i)>g(i-1) && g(i)>g(i+1))
            r(i) = g(i);
        else
            r(i) = 0;
        end

        if(g(i)<g(i-1) && g(i)<g(i+1))
            v(i) = g(i);
        else
            v(i) = 0;
        end
    end
    if r(i)>0
        w(i) = r(i);
    end
    
    if v(i)<0
        w(i) = -v(i);
    end

    if(w(i) == 1)
        atv = 1;
    end

    if ((atv ==1) && (data.x.zecg(i)>=0.03))
        atv = 2;
    end

    R_dtct(i) = 0;
    if ((atv ==2) && (data.x.zecg(i)<=0.03))
        j = 1;
        aux_E = 1;
        atv = 0;
        R_dtct(i) = 1;
    end

    if (i > ip)
        fprintf('Executing ... \t%d %%\r',lg);
        lg = lg + 10;
        ip = ip + (N-1)/10;
    end
end

%plot(data.t(1:end-1),data.x.zecg,'b', data.t(1:end-1),max(data.x.zecg)*R_dtct,'k');
r_wave_times = find(max(data.x.zecg)*R_dtct) * T_s;

%% Model functions and state space matrices
E_n = @(t_n) 1.55 ...
    * (((t_n/0.7)^1.9)/(1 + (t_n/0.7)^1.9)) ...
    * 1/(1 + (t_n/1.17)^21.9);

E = @(t) (t > r_wave_times(1))*(E_max - E_min) ...
    * E_n(mod(t, ...
    r_wave_times(abs(t - r_wave_times) == min(abs(t - r_wave_times)))... % Technical replacement for argmin :)
    )/t_max) ...
    + E_min;

P_ve = @(E_t, V_ve) E_t * (V_ve - V_0);

% x = [V_ve, P_ae, Q_a, P_ao, P_s, Q_i, Q_o, P_d, V_c]
A = @(E_t, D_m, D_a, beta_i, beta_o, gamma, omega_i_t, omega_o_t) [-(D_m/R_m + D_a/R_a) * E_t, D_m/R_m, 0, D_a/R_a, 0, -1, 0, 0, 0;
    (D_m * E_t)/(R_m * C_ae), -(C_ae^-1) * (1/R_s + D_m/R_m), 0, 0, (R_s * C_ae)^-1, 0, 0, 0, 0;
    0, 0, -R_c/L_s, L_s^-1, -L_s^-1, 0, 0, 0, 0;
    (D_a * E_t)/(R_a * C_ao), 0, -C_ao^-1, -D_a/(R_a * C_ao), 0, 0, C_ao^-1, 0, 0;
    0, (R_s * C_s)^-1, C_s^-1, 0, -(R_s * C_s)^-1, 0, 0, 0, 0;
    beta_i * E_t / (1 - gamma * L_p), 0, 0, - gamma / (1 - gamma * L_p), 0, - (alpha_PVAD * (beta_i - gamma) + (beta_i * R_p + omega_i_t + gamma * R_p))/(1 - gamma * L_p), (alpha_PVAD * (beta_i - gamma) + (beta_i * R_p - gamma * R_p - beta_i * L_p * omega_o_t))/(1 - gamma * L_p), -(beta_i - gamma)/(1 - gamma * L_p), -(beta_i - gamma)/(C_p * (1 - gamma * L_p));
    0, 0, 0, -beta_o / (1 - gamma * L_p), 0, (alpha_PVAD * (beta_o - gamma) + (beta_o * R_p - beta_o * L_p * omega_i_t - gamma * R_p))/(1 - gamma * L_p), -(alpha_PVAD * (beta_o - gamma) + (beta_o * R_p + omega_o_t - gamma * R_p))/(1 - gamma * L_p), (beta_o - gamma) / (1 - gamma * L_p), (beta_o - gamma)/(C_p * (1 - gamma * L_p));
    0, 0, 0, 0, 0, 0, 0, -(R_d * C_d)^-1, 0;
    0, 0, 0, 0, 0, 1, -1, 0, 0];

p = @(E_t, D_m, D_a, beta_i, beta_o, gamma, P_x) [(D_m/R_m + D_a/R_a)*E_t*V_0;
    -(D_m*E_t*V_0)/(R_m*C_ae);
    0;
    -(D_a*E_t*V_0)/(R_a*C_ao);
    0;
    ((beta_i - gamma)/C_p) * Vd_vad - beta_i * E_t * V_0;
    ((beta_o - gamma)/C_p) * Vd_vad;
    P_x / (R_d * C_d);
    0];

%% Preparing main simulation
data.x.V_ve = zeros(N, 1);
data.x.P_ae = zeros(N, 1);
data.x.Q_a = zeros(N, 1);
data.x.P_ao = zeros(N, 1);
data.x.P_s = zeros(N, 1);
data.x.Q_i = zeros(N, 1);
data.x.Q_o = zeros(N, 1);
data.x.P_d = zeros(N, 1);
data.x.V_c = zeros(N, 1);
data.P_ve = zeros(N, 1);
data.E = zeros(N, 1);
data.omega_i = zeros(N, 1);
data.omega_o = zeros(N, 1);
data.P_x = zeros(N, 1);

% Initial conditions
data.x.V_ve(1) = 140;
data.x.P_ae(1) = 5;
data.x.Q_a(1) = 0;
data.x.P_ao(1) = 90;
data.x.P_s(1) = 90;
data.E(1) = E(0);
data.P_ve(1) = P_ve(data.E(1), data.x.V_ve(1));

D_a = 1;
D_m = 1;
D_i = 1;
D_o = 1;

ejection_counter = 0;

%% Running main simulation
n = 1;
t = 0;
while t < t_end
    
    % Calculating PVAD pressures
    if ejection_counter > 0
        P_x = P_e;
        ejection_counter = ejection_counter - 1;
    elseif ismember(n, (r_wave_times / T_s) + (ejection_delay / T_s))
        P_x = P_e;
        ejection_counter = ejection_time / T_s;
    else 
        P_x = P_f;
    end
    
    P_air = data.x.P_d(n) + alpha_PVAD * (data.x.Q_i(n) - data.x.Q_o(n));
    P_c = P_air;
    
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
    
    if (data.P_ve(n) > P_c)
        D_i = 1;
        D_o = 0.03;
    elseif (P_c > data.x.P_ao(n))
        D_i = 0.08;
        D_o = 1;
    else
        D_i = 0.08;
        D_o = 0.03;
    end
        
    % State vector
    x = [data.x.V_ve(n); data.x.P_ae(n); data.x.Q_a(n); data.x.P_ao(n); data.x.P_s(n); data.x.Q_i(n); data.x.Q_o(n); data.x.P_d(n); data.x.V_c(n)];
    
    % Updating PVAD values
    beta_i = D_i/(L_i + D_i * L_p);
    beta_o = D_o/(L_o + D_o * L_p);
    gamma = beta_i * L_p * beta_o;
    R_i_new = R_i + exp(-0.25*data.x.V_ve(n));
    omega_i_t = R_i_new * (L_i + D_i * L_p);
    omega_o_t = R_o(data.x.Q_o(n)) * (L_o + D_o * L_p);
        
    % RK4 integration
    A_i = A(data.E(n), D_m, D_a, beta_i, beta_o, gamma, omega_i_t, omega_o_t);
    p_i = p(data.E(n), D_m, D_a, beta_i, beta_o, gamma, P_x);
    
    dx =  A_i * x +  p_i;
    kx1 = T_s * dx;
    x1 = x + 0.5*kx1;
    
    dx = A_i * x1 +  p_i;
    kx2 = T_s * dx;
    x1 = x + 0.5*kx2;
    
    dx = A_i * x1 +  p_i;
    kx3 = T_s * dx;
    x1 = x + kx3;
    
    dx = A_i * x1 +  p_i;
    kx4 = T_s * dx;
    
    xf = x + (kx1 + 2 * kx2 + 2 * kx3 + kx4)/6;
    
    % Iterating
    t = t + T_s;
    n = n + 1;
    
    % Saving data
    data.x.V_ve(n) = xf(1);
    data.x.P_ae(n) = xf(2);
    data.x.Q_a(n) = xf(3);
    data.x.P_ao(n) = xf(4);
    data.x.P_s(n) = xf(5);
    data.x.Q_i(n) = xf(6);
    data.x.Q_o(n) = xf(7);
    data.x.P_d(n) = xf(8);
    data.x.V_c(n) = xf(9);
    
    data.P_x(n) = P_x;    
    data.E(n) = E(t);
    data.omega_i(n) = omega_i_t;
    data.omega_o(n) = omega_o_t;
    
    data.P_ve(n) = P_ve(data.E(n), data.x.V_ve(n));
    
end

%% Ploting results
figure(1);
subplot(3,1,1);
plot(data.t, data.x.P_ao, 'k', data.t, data.P_ve, ':k', data.t, data.x.P_ae, '--k');
title('Results - Glauber R. Leite');
legend('AoP', 'LVP', 'LAP');
ylabel('Pressures (mmHg)');
%axis([0, 3, -5, 150]);
grid on;

subplot(3,1,2);
plot(data.t, data.x.Q_a, 'k');
ylabel('Aortic Flow (ml/s)');
%axis([0, 3, -50, 700]);
grid on;

subplot(3,1,3);
plot(data.t, data.x.V_ve, 'k');
ylabel('LVV (ml)');
xlabel('time (s)');
%axis([0, 3, 50, 160]);
grid on;

figure(2);
subplot(3,1,1);
plot(data.t(1:end-1),data.x.zecg,'b', data.t(1:end-1),max(data.x.zecg)*R_dtct,'k');
title('MaMeMi performance');
legend('ECG', 'MaMeMi');
axis([0, 3, -0.03, 0.06]);
grid on;
 
subplot(3,1,2);
plot(data.t(1:end-1),data.E(1:end-1),'b', data.t(1:end-1),E_max * 20*max(data.x.zecg)*R_dtct,'k');
title('Elastance');
legend('E(t)', 'MaMeMi');
axis([0, 3, -0.2, E_max + 0.2]);
grid on;
 
subplot(3,1,3);
plot(data.t(1:end-1),data.P_x(1:end-1),'b', data.t(1:end-1), P_e * 20*max(data.x.zecg)*R_dtct,'k');
title('Pneumatic Signal');
legend('P_x', 'MaMeMi');
axis([0, 3, P_f - 20, P_e + 20]);
grid on;

%figure(3);