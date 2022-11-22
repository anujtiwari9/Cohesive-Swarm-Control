clear all 
clc
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nfig=0;
video_steps = 5;
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



%% damping ratio as a function of length of the network

%for L <= Lstart
L_mat = [1:1:10^4]*a;
c = sqrt( (gamma_dsr*a^2)/(2*D*beta2_dsr*dt) );

Lstar_by_a = beta1_star*pi*c./(4*gamma_dsr*a)
L_by_a = beta1_dsr*pi*c./(4*gamma_dsr*a)


zeta_mat_0 = beta1_star*pi*c./(4*gamma_dsr*L_mat);
zeta_mat = beta1_dsr*pi*c./(4*gamma_dsr*L_mat);
zeta_mat_2 = beta1_dsr_2*pi*c./(4*gamma_dsr*L_mat);
zeta_mat_3 = beta1_dsr_3*pi*c./(4*gamma_dsr*L_mat);

vertical_line = 0.000001:0.1:3;

nfig=nfig+1; figure(nfig);
loglog(L_mat, zeta_mat_0, 'm:', 'LineWidth', 3);
hold on
loglog(L_mat, zeta_mat, 'b-', 'LineWidth', 3);
hold on
loglog(L_mat, zeta_mat_2, 'k-', 'LineWidth', 3);
hold on
loglog(L_mat, zeta_mat_3, 'r--', 'LineWidth', 3);
hold on
loglog(L_mat, ones(size(L_mat)), 'k:', 'LineWidth', 2);
hold on 
loglog(Lstar_by_a*a*ones(size(vertical_line)), vertical_line, 'k:', 'LineWidth', 3);

%text(15*a-1.5, vertical_line(end)+2, 'L=15 \times a', 'FontSize', 20);

%loglog(L*ones(size(vertical_line)), vertical_line, 'k:', 'LineWidth', 3);

%text(L-200, vertical_line(end)+2, ' L = 2000 \times a', 'FontSize', 20);

xlabel('Network length (L = N \times a)');
ylabel('Damping ratio \zeta (with \beta_2 = 1)');
%grid on
set(gca, 'FontSize', 24);
legend('\beta_1 = \beta_1^*', '\beta_1 = 0.9 \times \beta_1^*', '\beta_1 = 0.1 \times \beta_1^*', '\beta_1 = 0.01 \times \beta_1^*', 'Location','northeastoutside')



%for L > Lstar
%1 - beta_2 with both viscous and internal damping
one_minus_beta2_mat = zeros(size(L_mat));
Ts_mat = one_minus_beta2_mat;


for i=1:1:length(L_mat)

    L = L_mat(i);
    
    if (L <= Lstar_by_a*a)
        
        beta2_dsr = 1;

        c = sqrt( (gamma_dsr*a^2)/(2*D*beta2_dsr*dt) );

        omega_0 = pi*c/(2*L);

        Ts = 6/(1*omega_0);

        Ts_mat(i) = Ts;
        

    else

        L = L_mat(i);

        beta2_dsr = (pi^2*a^2/(4*L^2))*( sqrt( gamma_dsr*dt/(2*D) + 4*(L^2)/(pi^2*a^2) + beta1_dsr/(2*D) ) - sqrt( gamma_dsr*dt/(2*D)) )^2;
        
        c = sqrt( (gamma_dsr*a^2)/(2*D*beta2_dsr*dt) );
        
        zeta_dsr = (1-beta2_dsr)*L/(pi*c*beta2_dsr*dt) + beta1_dsr*pi*c/(4*gamma_dsr*L);

        one_minus_beta2_mat(i) = 1 - beta2_dsr;

        omega_0 = pi*c/(2*L);

        Ts = 6/(zeta_dsr*omega_0);

        Ts_mat(i) = Ts;

    end


end


vertical_line = 0.000001:0.0001:1.2*10^(-3);

% nfig=nfig+1; figure(nfig);
% yyaxis left
% semilogx(L_mat, one_minus_beta2_mat, '-', 'LineWidth', 3);
% hold on;
% semilogx(Lstar_by_a*a*ones(size(vertical_line)), vertical_line, 'k:', 'LineWidth',3)
% hold on;
% semilogx(500*a*ones(size(vertical_line)), vertical_line, 'r:', 'LineWidth',3)
% ylabel('Viscous damping (1-\beta_2)')
% %hold on
% yyaxis right
% loglog(L_mat, Ts_mat, '--', 'LineWidth', 3);
% ylabel('Settling time (T_s)')
% xlabel('Network length (L = N \times a)')
% set(gca, 'FontSize', 24)



nfig=nfig+1; figure(nfig);
semilogx(L_mat, one_minus_beta2_mat, '-', 'LineWidth', 4);
hold on;
semilogx(Lstar_by_a*a*ones(size(vertical_line)), vertical_line, 'k:', 'LineWidth',4)
hold on;
% semilogx(1000*a*ones(size(vertical_line)), vertical_line, 'k:', 'LineWidth',4)
% hold on;
semilogx(2000*a*ones(size(vertical_line)), vertical_line, 'k:', 'LineWidth',4)
ylabel('Viscous damping (1-\beta_2)')
%hold on
xlabel('Network length (L = N \times a)')
set(gca, 'FontSize', 20)
saveas(gcf,'Fig_viscous_damping_vs_L','epsc')