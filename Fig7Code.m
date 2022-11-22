clear all 
clc
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nfig=0;
video_steps = 100; %200
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%simulations 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L = 3;     %length of network in m
a = 0.2; 
khat = 1;    %stiffness between agents
f1 = 0.08; f2 = 0.1; f3 = 10;   %Hz source frequency
disp(['Simulating frequencies: ' num2str(f1) ' Hz, ' num2str(f2) ' Hz, ' num2str(f3) ' Hz.'])
%wave properties

T1 = 1/f1;                  %time period of propagated in seconds
T2 = 1/f2;
T3 = 1/f3;

cases_to_run = [1 1 0 1];   %[DSR Mom-small-damping Mom-large-damping Nesterov]

%topological radius
threshold_dist = 2*a;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DSR parameters
gamma_dsr = 10;
beta2_dsr = 1;   %=1 for DSR
delta_t = 10^(-4); 
dt = delta_t;

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


if (L <= 17.55*a)
        %%%selecting stable beta
    
    % % % 
%        beta2_dsr = (pi^2*a^2/(4*L^2))*( sqrt( gamma_dsr*dt/(2*D) + 4*(L^2)/(pi^2*a^2) + beta1_dsr/(2*D) ) - sqrt( gamma_dsr*dt/(2*D)) )^2
     c = sqrt(gamma_dsr*a^2/(2*D*delta_t*beta2_dsr));
    
       beta1_dsr = 1; %4*gamma_dsr*L/(pi*c)

else    

    %%%selecting stable beta1
    beta1_dsr = 0.9*((beta2_dsr+1) - gamma_dsr*delta_t/2 )/max(lambda_A);

    
    
    
    % % % 
       beta2_dsr = 1;%(pi^2*a^2/(4*L^2))*( sqrt( gamma_dsr*dt/(2*D) + 4*(L^2)/(pi^2*a^2) + beta1_dsr/(2*D) ) - sqrt( gamma_dsr*dt/(2*D)) )^2
     c = sqrt(gamma_dsr*a^2/(2*D*delta_t*beta2_dsr));
    
    %   beta1_dsr = 4*gamma_dsr*L/(pi*c)

    video_name = ['Oct_27_2Dsimuls_longnetwork3.mp4'];



end

zeta_dsr = (1-beta2_dsr)*L/(pi*c*beta2_dsr*dt) + beta1_dsr*pi*c/(4*gamma_dsr*L)
omega_0 = pi*c/(2*L);
predicted_settling_time = 6/(zeta_dsr*omega_0)


%Momentum parameter
gamma_m = gamma_dsr;
beta2_m_1 = 0.999;
beta2_m_2 = 1;
beta1_m = 0;

%Nesterov parameter
gamma_n = gamma_dsr;
beta2_n = 0.999;
beta1_n = gamma_n*dt*beta2_n*(1+beta2_n);
beta2_n_2 = 1;
beta1_n_2 = gamma_n*dt*beta2_n_2*(1+beta2_n_2);


v = sqrt(gamma_dsr*a^2/(2*1*delta_t*beta2_dsr))        %wave velocity in m/s

tend = max(0.9*L/v, 20) %2.8835e+03 %20*L/v;           %s - duration of simulation time -- BEFORE REFLECTION 
% return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% 
% beta2_dsr = (pi^2*a^2/(4*L^2))*( sqrt( gamma_dsr*dt/(2*D) + 4*(L^2)/(pi^2*a^2) + beta1_dsr/(2*D) ) - sqrt( gamma_dsr*dt/(2*D)) )^2
%  c = sqrt(gamma_dsr*a^2/(2*D*delta_t*beta2_dsr));
%  zeta_dsr = (1-beta2_dsr)*L/(pi*c*beta2_dsr*dt) + beta1_dsr*pi*c/(4*gamma_dsr*L)
% omega_0 = pi*c/(2*L);
%  predicted_settling_time = 6/(zeta_dsr*omega_0)

predicted_settling_time = tend;
 
% 
%   return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%time variable
ts = 0:delta_t:T2/2; %time period of initial BC creating pulse
t = [ts (ts(end)+delta_t):delta_t:tend]; %time vector

noise_amp_mat = 0:0.1:1;

velocity_alignment_mat = zeros(3,length(noise_amp_mat));

for kkk=1:1:length(noise_amp_mat) 

    kkk/length(noise_amp_mat)

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %define input source
            input1 = sin(2*pi*ts/T2);
            input2 = noise_amp_mat(kkk)*sin(2*pi*ts/T3);
            input = input1 + input2;
            initialzeros = 10000; %2
            Is1 = [zeros(1,initialzeros) input1 zeros(1, length(t)-length(ts)-initialzeros)];
            Is2 = [zeros(1,initialzeros) input2 zeros(1, length(t)-length(ts)-initialzeros)];
            Is = [zeros(1,initialzeros) input zeros(1, length(t)-length(ts)-initialzeros)];
            % Is_1 = ones(size(t)); %sin(2*pi*f1*t);
            % Is_2 = sin(2*pi*f2*t);
            % Is_3 = sin(2*pi*f3*t);
            
            %  return




            %%%%%%%%METHOD 2: DSR %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %number of discrete x points
            num_of_x = n; %equal to number of agents in 1d network
            
            %number of time steps 
            num_of_t = length(t);
            
            Idsr1 = zeros(num_of_x, 3);  
            Im_1 = zeros(num_of_x,3);
            Im_2 = zeros(num_of_x, 3);
            Im_3 = zeros(num_of_x, 3);
            
            %DSR (Internal damping)
            thetadsr = 0*pi/2*ones(num_of_x, 3); %initially the angle is pi/2
            xdsr = linspace(a, num_of_x*a, num_of_x);
            ydsr = 0*xdsr;
            vel = 1; 
             leader_xytheta_dsr = [xdsr(1) ydsr(1) thetadsr(1,1);
                                   xdsr(1) ydsr(1) thetadsr(1,2)];
             follower_xytheta_dsr = [xdsr(end) ydsr(end) thetadsr(end,1);
                                   xdsr(end) ydsr(end) thetadsr(end,2)];
             time_vec = [t(1:2)];

             %M1 (viscous damping)
             thetam1 = thetadsr; xm1 = xdsr; ym1 = ydsr;
             leader_xytheta_m1 = leader_xytheta_dsr;
             follower_xytheta_m1 = follower_xytheta_dsr;

             %M2 (undamped)
             thetam2 = thetadsr; xm2 = xdsr; ym2 = ydsr;
             leader_xytheta_m2 = leader_xytheta_dsr;
             follower_xytheta_m2 = follower_xytheta_dsr;




       %%%%%%%%%%%%%%
       theta_dsr_mat = [];
       theta_m1_mat = [];
       theta_m2_mat = [];

            
            for k = [3:1:num_of_t] 
                %t(k)/t(end)

                iii=3;

                %frequency f1
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

               bigA = get_topological_A_matrix(xdsr, ydsr, threshold_dist, n);

               Delta_kminus1 = bigA*[Idsr1(:,iii-1); Is(k-1)];
               Delta_kminus2 = bigA*[Idsr1(:,iii-2); Is(k-2)];

               Idsr1(:,iii) = Idsr1(:,iii-1) - gamma_dsr*delta_t*Delta_kminus1 - beta1_dsr*(Delta_kminus1-Delta_kminus2) ...
                                    + beta2_dsr*(Idsr1(:,iii-1)-Idsr1(:,iii-2));% + gamma_dsr*delta_t*B(i)*Is(k-1); 


                
                       thetadsr = Idsr1*pi/2; %thetadsr + pi/2*Idsr1;
                       xdsr = xdsr + vel*delta_t*cos(thetadsr(:,3)');
                       ydsr = ydsr + vel*delta_t*sin(thetadsr(:,3)');



                Idsr = Idsr1;
                f=f2;

                Idsr1(:,1:2) = Idsr1(:,2:3);


                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
               bigA = get_topological_A_matrix(xm1, ym1, threshold_dist, n);

               Delta_kminus1 = bigA*[Im_1(:,iii-1); Is(k-1)];
               Delta_kminus2 = bigA*[Im_1(:,iii-2); Is(k-2)];
                
               Im_1(:,iii) = Im_1(:,iii-1) - gamma_m*delta_t*Delta_kminus1 - beta1_m*(Delta_kminus1-Delta_kminus2) ...
                                + beta2_m_1*(Im_1(:,iii-1)-Im_1(:,iii-2));




                   diag_A_m1 = diag(bigA(1:n,1:n));

                   for i=1:1:n

                       if (diag_A_m1(i) == 0)

                        xm1(i) = xm1(i);% + vel*delta_t*cos(thetam2(:,3)');
                        ym1(i) = ym1(i);% + vel*delta_t*sin(thetam2(:,3)');
    
                       else

                         thetam1(i,:) = Im_1(i,:)*pi/2; %thetadsr + pi/2*Idsr1;
                         xm1(i) = xm1(i) + vel*delta_t*cos(thetam1(i,3)');
                         ym1(i) = ym1(i) + vel*delta_t*sin(thetam1(i,3)');
    
                       end

                   end



                Im_1(:,1:2) = Im_1(:,2:3);



               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
               bigA = get_topological_A_matrix(xm2, ym2, threshold_dist, n);

               Delta_kminus1 = bigA*[Im_2(:,iii-1); Is(k-1)];
               Delta_kminus2 = bigA*[Im_2(:,iii-2); Is(k-2)];
                
            Im_2(:,iii) = Im_2(:,iii-1) - gamma_m*delta_t*Delta_kminus1 - beta1_m*(Delta_kminus1-Delta_kminus2) ...
                    + beta2_m_2*(Im_2(:,iii-1)-Im_2(:,iii-2));% + gamma_m*delta_t*B(i)*Is(k-1); 

                  
                 
                   diag_A_m2 = diag(bigA(1:n,1:n));

                   for i=1:1:n

                       if (diag_A_m2(i) == 0)

                        xm2(i) = xm2(i);% + vel*delta_t*cos(thetam2(:,3)');
                        ym2(i) = ym2(i);% + vel*delta_t*sin(thetam2(:,3)');
    
                       else

                         thetam2(i,:) = Im_2(i,:)*pi/2; %thetadsr + pi/2*Idsr1;
                         xm2(i) = xm2(i) + vel*delta_t*cos(thetam2(i,3)');
                         ym2(i) = ym2(i) + vel*delta_t*sin(thetam2(i,3)');
    
                       end

                   end


                Im_2(:,1:2) = Im_2(:,2:3);              


            
                
            
                if (mod(k,video_steps) == 0)


                        theta_dsr_mat = [theta_dsr_mat thetadsr(:,3)];
                        theta_m1_mat = [theta_m1_mat thetam1(:,3)];
                        theta_m2_mat = [theta_m2_mat thetam2(:,3)];

                          %%2D simul
                          leader_xytheta_dsr = [leader_xytheta_dsr; xdsr(1) ydsr(1) thetadsr(1,3)];
                          follower_xytheta_dsr = [follower_xytheta_dsr; xdsr(end) ydsr(end) thetadsr(end,3)];

                          leader_xytheta_m1 = [leader_xytheta_m1; xm1(1) ym1(1) thetam1(1,3)];
                          follower_xytheta_m1 = [follower_xytheta_m1; xm1(end) ym1(end) thetam1(end,3)];

                          leader_xytheta_m2 = [leader_xytheta_m2; xm2(1) ym2(1) thetam2(1,3)];
                          follower_xytheta_m2 = [follower_xytheta_m2; xm2(end) ym2(end) thetam2(end,3)];

                          time_vec = [time_vec t(k)];
                          
                          



                end

           


    
                        
            end



   
            %%%%%%%%DSR
            cos_sum_dsr = sum( cos(theta_dsr_mat) );
            sin_sum_dsr = sum( sin(theta_dsr_mat) );

            velocity_alignment_mat(1,kkk) = min(sqrt( cos_sum_dsr.^2 + sin_sum_dsr.^2 )/n);

            %%%%%%%m1
            cos_sum_m1 = sum( cos(theta_m1_mat) );
            sin_sum_m1 = sum( sin(theta_m1_mat) );

            velocity_alignment_mat(2,kkk) = min(sqrt( cos_sum_m1.^2 + sin_sum_m1.^2 )/n);

            %%%%%%m2
            cos_sum_m2 = sum( cos(theta_m2_mat) );
            sin_sum_m2 = sum( sin(theta_m2_mat) );

            velocity_alignment_mat(3,kkk) = min(sqrt( cos_sum_m2.^2 + sin_sum_m2.^2 )/n);





end

%%

nfig=nfig+1; figure(nfig);
semilogx(noise_amp_mat, velocity_alignment_mat(1,:), 'd--', 'MarkerSize', 12, 'MarkerFaceColor', [0 0 1], 'LineWidth', 3, 'Color',[0 0 1]);
hold on
semilogy(noise_amp_mat, velocity_alignment_mat(2,:), 'd--', 'MarkerSize', 12, 'MarkerFaceColor', [1 0 0], 'LineWidth', 3, 'Color',[1 0 0]);
hold on
semilogy(noise_amp_mat, velocity_alignment_mat(3,:), 'd--', 'MarkerSize', 12, 'MarkerFaceColor', [0 1 0], 'LineWidth', 3, 'Color',[0 1 0]);
hold on
xlabel('Noise amplitude')
ylabel('Velocity alignment (v_a)')
 yticks([0 0.1 0.2:0.2:1.2])
yticklabels({'0','0.1','0.2','0.4','0.6','0.8', '1', '1.2'})
 xticks(0.0:0.2:1.2)
xticklabels({'0','0.2','0.4','0.6','0.8','1','1.2'})
legend('Internal damping', 'Viscous damping', 'Undamped', 'Location','southwest')
set(gca, 'FontSize', 24);
saveas(gcf,'Fig_velocity_alignment_vs_noise','epsc')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bigA = get_topological_A_matrix(xdsr, ydsr, radius, N)
%returns graph lapacian (N rows and columns where last N+1 column is for the source node
%and connected only to the leader agent) based on current x y location of agents with
%topological dist of given radius 
    bigA = zeros(N, N+1);

    pos_of_agents = sqrt(xdsr.^2+ydsr.^2);

    for i=1:1:N

            for j=1:1:N

                if ( j ~= i && abs(pos_of_agents(j)-pos_of_agents(i)) < radius )

                    bigA(i,j) = bigA(i,j) - 1;
%                     bigA(i,i) = bigA(i,i) + 1;
                        
                end    

            end

    end
    
    %leader agent
%     bigA(1,1) = bigA(1,1) + 1;
    bigA(1,N+1) = bigA(1, N+1) - 1;

    for i=1:1:N
            
            sum_of_row = - sum(bigA(i,:));
            
            if (sum_of_row > 0)
                bigA(i,:) = bigA(i,:)/sum_of_row;
            
                bigA(i,i) = 1;

            end

    end




end


