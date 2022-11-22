clear all 
clc
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nfig=0;
video_steps = 5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%simulations 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L = 400;                %length of network in m

gamma = 10;             %alignment strength

T = 10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Internal damping parameters
gamma_dsr = gamma;
beta2_dsr = 1;              %=1 for DSR
delta_t = 10^(-4); 
dt = delta_t;


%Viscous damping and undamped
beta2_m_1 = 0.999; beta1_m_1 = 0;%viscous damping

beta2_m_2 = 0.99999999; beta1_m_2 = 0; %undamped


%interagent distance and nominal wave speed for internal damping
a = 0.2; %0.0224; %sqrt(v^2*beta2*delta_t*2/gamma)
v = sqrt(gamma*a^2/(2*1*delta_t*beta2_dsr))        %wave velocity in m/s

tend = 0.95*L/v;           %s - duration of simulation time -- BEFORE REFLECTION 
% return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%constrained parameters 
D = 1; %number of spatial dimensions
n = 1*round(L/a) %number of agents


%return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%network definition
connection_struct = [-0.5 1 -0.5];
A = zeros(n,n);
A(1,1) = 1; A(1, 2) = -0.5; 
A(n,n) = 1; A(n,n-1) = -1;
for i=2:1:n-1
    A(i,i-1:1:i+1) = connection_struct;
end
A = A;
B = [0.5; zeros(n-1,1)];

lambda_A = eig(A);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta1_star = ((beta2_dsr+1) - gamma_dsr*delta_t/2 )/max(lambda_A);
%different choices of beta1
beta1_dsr = 0.9*beta1_star
beta1_dsr_2 = 0.1*beta1_star
beta1_dsr_3 = 0.01*beta1_star

%time variable
ts = 0:delta_t:T/2; %time period of initial BC creating pulse
t = [ts (ts(end)+delta_t):delta_t:tend]; %time vector
% return

%% creating theory vs simulation plots
% close all
% clc

%frequency vectors
w_mat = 2*pi*logspace(log(0.001)/log(10),log(10^2)/log(10), 100);

w_mat_Hz = w_mat/(2*pi);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %DSR with beta1 = 0.9 X beta1_star
    coeff_dsr = w_mat./( 1 + beta1_dsr^2*w_mat.^2/gamma_dsr^2 );
%     return
    realpart_dsr = ( w_mat - (1-beta2_dsr)*beta1_dsr*w_mat./(beta2_dsr*dt*gamma_dsr) );
    imagpart_dsr = -( beta1_dsr*w_mat.^2/gamma_dsr + (1-beta2_dsr)/(beta2_dsr*dt) );
    c2_k2_dsr = coeff_dsr.*complex(realpart_dsr, imagpart_dsr);
    c2_k2_dsr_approx = w_mat.^2;
%     return
    %taking square root
    a_real = coeff_dsr.*realpart_dsr; 
    b = coeff_dsr.*imagpart_dsr; 
    z_mod = sqrt(a_real.^2 + b.^2);
    c_k_dsr_positive = complex( sqrt( (z_mod+a_real)/2 ), (b./abs(b)).*sqrt( (z_mod-a_real)/2 ) ); 
    c_k_dsr_negative = -complex( sqrt( (z_mod+a_real)/2 ), (b./abs(b)).*sqrt( (z_mod-a_real)/2 ) );
    
    c = sqrt( (gamma*a^2)/(2*D*beta2_dsr*dt) );
    speed_of_travel_dsr  = ( -(c)*w_mat./real(c_k_dsr_negative) );
    speed_of_travel_dsr_lowfreq_approx = c*(1 + 3*beta1_dsr^2*w_mat.^2/(8*gamma_dsr^2));
    
    dissipation_sigma_dsr = - imag(c_k_dsr_negative)/c;
    dissipation_sigma_dsr_lowfreq_approx  = beta1_dsr*w_mat.^2/(2*gamma_dsr*c);
    dissipation_sigma_dsr_highfreq_approx = sqrt(gamma_dsr*w_mat/(2*beta1_dsr))/c;


    
    %comparing square root correctness

    sqrt_root_realpart_1 = real(-sqrt(c2_k2_dsr));
    sqrt_root_imagpart_1 = imag(-sqrt(c2_k2_dsr));

    sqrt_root_realpart_2 = real(c_k_dsr_negative);
    sqrt_root_imagpart_2 = imag(c_k_dsr_negative);

  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %DSR with beta1 = 0.1 X beta1_star
    coeff_dsr = w_mat./( 1 + beta1_dsr_2^2*w_mat.^2/gamma_dsr^2 );
%     return
    realpart_dsr = ( w_mat - (1-beta2_dsr)*beta1_dsr_2*w_mat./(beta2_dsr*dt*gamma_dsr) );
    imagpart_dsr = -( beta1_dsr_2*w_mat.^2/gamma_dsr + (1-beta2_dsr)/(beta2_dsr*dt) );
    c2_k2_dsr = coeff_dsr.*complex(realpart_dsr, imagpart_dsr);
    c2_k2_dsr_approx = w_mat.^2;
%     return
    %taking square root
    a_real = coeff_dsr.*realpart_dsr; 
    b = coeff_dsr.*imagpart_dsr; 
    z_mod = sqrt(a_real.^2 + b.^2);
    c_k_dsr_positive = complex( sqrt( (z_mod+a_real)/2 ), (b./abs(b)).*sqrt( (z_mod-a_real)/2 ) ); 
    c_k_dsr_negative = -complex( sqrt( (z_mod+a_real)/2 ), (b./abs(b)).*sqrt( (z_mod-a_real)/2 ) );
    
    c = sqrt( (gamma*a^2)/(2*D*beta2_dsr*dt) );
    speed_of_travel_dsr_2 = ( -(c)*w_mat./real(c_k_dsr_negative) );
    dissipation_sigma_dsr_2 = - imag(c_k_dsr_negative)/c;
%   return    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %DSR with beta1 = 0.01 X beta1_star
    coeff_dsr = w_mat./( 1 + beta1_dsr_3^2*w_mat.^2/gamma_dsr^2 );
%     return
    realpart_dsr = ( w_mat - (1-beta2_dsr)*beta1_dsr_3*w_mat./(beta2_dsr*dt*gamma_dsr) );
    imagpart_dsr = -( beta1_dsr_3*w_mat.^2/gamma_dsr + (1-beta2_dsr)/(beta2_dsr*dt) );
    c2_k2_dsr = coeff_dsr.*complex(realpart_dsr, imagpart_dsr);
    c2_k2_dsr_approx = w_mat.^2;
%     return
    %taking square root
    a_real = coeff_dsr.*realpart_dsr; 
    b = coeff_dsr.*imagpart_dsr; 
    z_mod = sqrt(a_real.^2 + b.^2);
    c_k_dsr_positive = complex( sqrt( (z_mod+a_real)/2 ), (b./abs(b)).*sqrt( (z_mod-a_real)/2 ) ); 
    c_k_dsr_negative = -complex( sqrt( (z_mod+a_real)/2 ), (b./abs(b)).*sqrt( (z_mod-a_real)/2 ) );
    
    c = sqrt( (gamma*a^2)/(2*D*beta2_dsr*dt) );
    speed_of_travel_dsr_3  = ( -(c)*w_mat./real(c_k_dsr_negative) );
    dissipation_sigma_dsr_3 = - imag(c_k_dsr_negative)/c;
%   return    


    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Viscous damping
    gamma_m = gamma;
    coeff_m = w_mat./( 1 + beta1_m_1^2*w_mat.^2/gamma_m^2 );
    realpart_m = ( w_mat - (1-beta2_m_1)*beta1_m_1*w_mat./(beta2_m_1*dt*gamma_m) );
    imagpart_m = -( beta1_m_1*w_mat.^2/gamma_m + (1-beta2_m_1)/(beta2_m_1*dt) );
    c2_k2_m = coeff_m.*complex(realpart_m, imagpart_m); 
    real_approx = (w_mat.^2)*( 1 - beta1_m_1*(1-beta2_m_1)/(gamma_m*beta2_m_1*dt) );
    imag_approx = - w_mat*(1-beta2_m_1)/(beta2_m_1*dt);
    c2_k2_m_approx = complex(real_approx, imag_approx);
    
    %taking square root
    a_real = coeff_m.*realpart_m; 
    b = coeff_m.*imagpart_m; 
    z_mod = sqrt(a_real.^2 + b.^2);
    c_k_m_positive = complex( sqrt( (z_mod+a_real)/2 ), (b./abs(b)).*sqrt( (z_mod-a_real)/2 ) ); 
    c_k_m_negative = -complex( sqrt( (z_mod+a_real)/2 ), (b./abs(b)).*sqrt( (z_mod-a_real)/2 ) ); 
    
    c = sqrt( (gamma*a^2)/(2*D*beta2_m_1*dt) );
    speed_of_travel_m_1  =  ( -(c)*w_mat./real(c_k_m_negative) );
    dissipation_sigma_m_1 = - imag(c_k_m_negative)/c;


 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Undamped
    gamma_m = gamma;
    coeff_m = w_mat./( 1 + beta1_m_2^2*w_mat.^2/gamma_m^2 );
    realpart_m = ( w_mat - (1-beta2_m_2)*beta1_m_2*w_mat./(beta2_m_2*dt*gamma_m) );
    imagpart_m = -( beta1_m_2*w_mat.^2/gamma_m + (1-beta2_m_2)/(beta2_m_2*dt) );
    c2_k2_m = coeff_m.*complex(realpart_m, imagpart_m); 
    real_approx = (w_mat.^2)*( 1 - beta1_m_2*(1-beta2_m_2)/(gamma_m*beta2_m_2*dt) );
    imag_approx = - w_mat*(1-beta2_m_2)/(beta2_m_2*dt);
    c2_k2_m_approx = complex(real_approx, imag_approx);
    
    %taking square root
    a_real = coeff_m.*realpart_m; 
    b = coeff_m.*imagpart_m; 
    z_mod = sqrt(a_real.^2 + b.^2);
    c_k_m_positive = complex( sqrt( (z_mod+a_real)/2 ), (b./abs(b)).*sqrt( (z_mod-a_real)/2 ) ); 
    c_k_m_negative = -complex( sqrt( (z_mod+a_real)/2 ), (b./abs(b)).*sqrt( (z_mod-a_real)/2 ) ); 
    
    c = sqrt( (gamma_m*a^2)/(2*D*beta2_m_2*dt) );
    speed_of_travel_m_2  =  ( -(c)*w_mat./real(c_k_m_negative) );
    dissipation_sigma_m_2 = - imag(c_k_m_negative)/c;





% 
% end

% return
% f_dsr = [0.01 0.1 0.3];
% simul_speed_dsr = [4.2222 4.5333 4.5333];
% simul_sigma_dsr = [];
% 
% f_m_1 = [0.005 0.01 0.1];
% simul_speed_m_1 = [4.4104 4.4403 4.4627];
% simul_sigma_m_1 = [];
% 
% f_m_2 = [];
% simul_speed_m_2 = [];
% simul_sigma_m_2 = [];
% 
% 
% f_n = [0.01];
% simul_speed_n = [3.7303];
% simul_sigma_n = [];

%Internal damping
f_dsr = [0.08 0.1 4 10];%4];
simul_speed_dsr = [44.7693 44.7898 89.01 98.04];%89.01];
simul_sigma_dsr = [-2.5378*10^(-4) -3.9625*10^(-4) -0.15408 -0.2723];%-0.3205]; %-0.15408]; 

%Viscous damping
f_m_1 = [0.08 0.1 4 10];%4];
simul_speed_m_1 = [14.9049 16.5373 43.6915 48];%43.6915];
simul_sigma_m_1 = [-0.029807 -0.033063 -0.10209 -0.1090];%-0.10209];

%Undamped
f_m_2 = [0.08 0.1 4 10];%4];
simul_speed_m_2 = v*ones(1,4);%[45.1405 44.5455 44.5455 44.5455];
simul_sigma_m_2 = [3.9017e-08 1.1265e-07 -3.9074e-05 -3.9074e-03];



 






vertical_line = [0:0.001:0.99];
nfig=nfig+1; ff=figure(nfig);
subplot(1,2,1)
loglog(w_mat_Hz, exp(dissipation_sigma_m_2), 'g-', w_mat_Hz, exp(dissipation_sigma_dsr), 'b--',  ...
    w_mat_Hz, exp(dissipation_sigma_m_1), 'r--',  ...
    f_m_2,  exp(simul_sigma_m_2), 'gd', ...
      f_dsr,  exp(simul_sigma_dsr), 'bd', f_m_1,  exp(simul_sigma_m_1), 'rd',  'LineWidth', 3, 'MarkerSize', 14);
hold on
loglog(0.5031, 0.99, 'k*', 0.5031*ones(size(vertical_line)), vertical_line, 'k--',  'LineWidth', 3, 'MarkerSize', 14)
% legend('Undamped', 'Internal damping', 'Viscous damping', 'Location', 'northhwest')
xlabel('Frequency (Hz)')
ylabel('Propagated Info. Amplitude (e^{-k_i})')
title('(a)')
vertical_line = [0:0.001:1.5*c];
set(gca, 'FontSize', 20);
ylim([0.5 1.6])
xticks([10^(-3) 10^(-2) 10^(-1) 10^(0) 10^(1) 10^(2) 10^(3)])
xticklabels({'0.001','0.01','0.1','1','10','100','1000'})
 grid on
subplot(1,2,2)
loglog(w_mat_Hz, speed_of_travel_m_2, 'g-',  ...
    w_mat_Hz, speed_of_travel_dsr, 'b--', ...
    w_mat_Hz, speed_of_travel_m_1, 'r--',  ...
    f_m_2, simul_speed_m_2, 'gd', ...
    f_dsr, simul_speed_dsr, 'bd', ...
    f_m_1, simul_speed_m_1, 'rd', ...
     'LineWidth', 3, 'MarkerSize', 14);
hold on
loglog(0.9134, 1.1*c, 'k*', 0.9134*ones(size(vertical_line)), vertical_line, 'k--',  'LineWidth', 3, 'MarkerSize', 14)
 legend('Undamped', 'Internal damping', 'Viscous damping', 'Location', 'southwest')
xlabel('Frequency (Hz)')
ylabel('Speed of Info. Propagation')
title('(b)')
 ylim([10^(-2) 10^(4)])
 xticks([10^(-3) 10^(-2) 10^(-1) 10^(0) 10^(1) 10^(2) 10^(3)])
xticklabels({'0.001','0.01','0.1','1','10','100','1000'})
grid on
set(gca, 'FontSize', 20);
ff.Position = [100 100 540*2.2 500];

saveas(gcf,'fig_theoreticalplots_v2','epsc')