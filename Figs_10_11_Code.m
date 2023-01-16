clear all 
clc
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nfig=0; 
video_steps = 500; %200
take_video = 0; %1 to save video of individual cases
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

% beta1_dsr = 0*4/(max(lambda_A)*(gamma_dsr*delta_t+2));

video_name = ['Nov17-MovieS2-white.mp4'];

if (L <= 17.55*a)
        %%%selecting stable beta1
    
    
    
    % % % 
%        beta2_dsr = (pi^2*a^2/(4*L^2))*( sqrt( gamma_dsr*dt/(2*D) + 4*(L^2)/(pi^2*a^2) + beta1_dsr/(2*D) ) - sqrt( gamma_dsr*dt/(2*D)) )^2
     c = sqrt(gamma_dsr*a^2/(2*D*delta_t*beta2_dsr));
    
       beta1_dsr = 0.9*((beta2_dsr+1) - gamma_dsr*delta_t/2 )/max(lambda_A) %4*gamma_dsr*L/(pi*c)

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

tend = max(0.9*L/v, 6) %2.8835e+03 %20*L/v;           %s - duration of simulation time -- BEFORE REFLECTION 
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


%topological radius
threshold_dist = 2*a;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%define input source
input1 = sin(2*pi*ts/T2);
input2 = 0*sin(2*pi*ts/T3);
input = input1 + input2;
initialzeros = 500; %2
Is1 = [zeros(1,initialzeros) input1 zeros(1, length(t)-length(ts)-initialzeros)];
Is2 = [zeros(1,initialzeros) input2 zeros(1, length(t)-length(ts)-initialzeros)];
Is = [zeros(1,initialzeros) input zeros(1, length(t)-length(ts)-initialzeros)];
% Is_1 = ones(size(t)); %sin(2*pi*f1*t);
% Is_2 = sin(2*pi*f2*t);
% Is_3 = sin(2*pi*f3*t);

Is_unitstep = ones(size(t));

% nfig=nfig+1; figure(nfig);
% plot(t, Is_1, '-', t, Is_2, '-', t, Is_3, '-', 'LineWidth', 3);
% xlabel('Time (s)');
% ylabel('Is');
% title('Source Trajectory');

nfig=nfig+1; figure(nfig);
plot(t, Is1, '-','LineWidth', 3);
xlabel('Time (s)');
ylabel('Source Signal (0.1 Hz)');
grid on
set(gca, 'FontSize', 24);
ylim([-1 2])

nfig=nfig+1; figure(nfig);
plot(t, Is2, '-','LineWidth', 3);
xlabel('Time (s)');
ylabel('Source Signal (10 Hz)');
grid on
set(gca, 'FontSize', 24);
ylim([-1 2]);

nfig=nfig+1; figure(nfig);
plot(t, Is, '-','LineWidth', 3);
xlabel('Time (s)');
ylabel('Source Signal (I_s)');
grid on
set(gca, 'FontSize', 24);
ylim([-1 2])

%  return



if (cases_to_run(1)==1)


            %%%%%%%%METHOD 2: DSR %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %number of discrete x points
            num_of_x = n; %equal to number of agents in 1d network
            
            %number of time steps 
            num_of_t = length(t);
            
            Idsr1 = zeros(num_of_x, 3); 
            Idsr2 = zeros(num_of_x, 3);
            
            %DSR (Internal damping) with tracking (DAC)
            thetadsr = 0*pi/2*ones(num_of_x, 3); %initially the angle is pi/2
            xdsr = linspace(a, num_of_x*a, num_of_x);
            ydsr = 0*xdsr;
            vel = 1; 
             leader_xytheta_dsr = [xdsr(1) ydsr(1) thetadsr(1,1);
                                   xdsr(1) ydsr(1) thetadsr(1,2)];
             follower_xytheta_dsr = [xdsr(end) ydsr(end) thetadsr(end,1);
                                   xdsr(end) ydsr(end) thetadsr(end,2)];
             time_vec = [t(1:2)];

             %DSR (Internal damping) without tracking (DAC)
             thetadsr2 = thetadsr; xdsr2 = xdsr; ydsr2 = ydsr;
             leader_xytheta_dsr2 = leader_xytheta_dsr;
             follower_xytheta_dsr2 = follower_xytheta_dsr;

             theta_desired = [Is(1)*pi/2; Is(2)*pi/2];

     
             


        if (take_video == 1)
             vidObj = VideoWriter(video_name, 'MPEG-4'); open(vidObj);
        end

            
            for k = [3:1:num_of_t] 
                t(k)/t(end)

                iii=3;

                %DSR Internal damping with tracking (DAC)
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

               bigA = get_topological_A_matrix(xdsr, ydsr, threshold_dist, n);

               Delta_kminus1 = bigA*[Idsr1(:,iii-1); Is(k-1)];
               Delta_kminus2 = bigA*[Idsr1(:,iii-2); Is(k-2)];

               Idsr1(:,iii) = Idsr1(:,iii-1) - gamma_dsr*delta_t*Delta_kminus1 - beta1_dsr*(Delta_kminus1-Delta_kminus2) ...
                                    + beta2_dsr*(Idsr1(:,iii-1)-Idsr1(:,iii-2));% + gamma_dsr*delta_t*B(i)*Is(k-1); 


                 diag_A_dsr = diag(bigA(1:n,1:n));


                 for i=1:1:n

                    %                      if (diag_A_dsr(i)==0)
                    % 
                    %                         xdsr(i) = xdsr(i);
                    %                         ydsr(i) = ydsr(i);
                    % 
                    %                      else

                       thetadsr(i,:) = Idsr1(i,:)*pi/2; %thetadsr + pi/2*Idsr1;
                       xdsr(i) = xdsr(i) + vel*delta_t*cos(thetadsr(i,3)');
                       ydsr(i) = ydsr(i) + vel*delta_t*sin(thetadsr(i,3)');

                    %                      end


                 end      



                Idsr = Idsr1;
                f=f2;

                Idsr1(:,1:2) = Idsr1(:,2:3);


                %DSR Internal damping without tracking (DAC)
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

               bigA = get_topological_A_matrix(xdsr2, ydsr2, threshold_dist, n);

               Delta_kminus1 = bigA*[Idsr2(:,iii-1); Is(k-1)];
               Delta_kminus2 = bigA*[Idsr2(:,iii-2); Is(k-1)];

               Idsr2(:,iii) = Idsr2(:,iii-1) - gamma_dsr*delta_t*Delta_kminus1 - beta1_dsr*(Delta_kminus1-Delta_kminus2) ...
                                    + beta2_dsr*(Idsr2(:,iii-1)-Idsr2(:,iii-2));


                 diag_A_dsr2 = diag(bigA(1:n,1:n));


                 for i=1:1:n

                    %                      if (diag_A_dsr(i)==0)
                    % 
                    %                         xdsr(i) = xdsr(i);
                    %                         ydsr(i) = ydsr(i);
                    % 
                    %                      else

                       thetadsr2(i,:) = Idsr2(i,:)*pi/2; %thetadsr + pi/2*Idsr1;
                       xdsr2(i) = xdsr2(i) + vel*delta_t*cos(thetadsr2(i,3)');
                       ydsr2(i) = ydsr2(i) + vel*delta_t*sin(thetadsr2(i,3)');

                    %                      end


                 end      

                Idsr2(:,1:2) = Idsr2(:,2:3);




            
                
            
                if (mod(k,video_steps) == 0)


                          %%2D simul
                          leader_xytheta_dsr = [leader_xytheta_dsr; xdsr(1) ydsr(1) thetadsr(1,3)];
                          follower_xytheta_dsr = [follower_xytheta_dsr; xdsr(end) ydsr(end) thetadsr(end,3)];

                          leader_xytheta_dsr2 = [leader_xytheta_dsr2; xdsr2(1) ydsr2(1) thetadsr2(1,3)];
                          follower_xytheta_dsr2 = [follower_xytheta_dsr2; xdsr2(end) ydsr2(end) thetadsr2(end,3)];
                          
                          time_vec = [time_vec t(k)];

                          theta_desired = [theta_desired; Is(k)*pi/2];
                          
                          
                          ff = figure(100);
                          subplot(1,2,1);
                          plot(leader_xytheta_dsr(:,1), leader_xytheta_dsr(:,2), '-', 'LineWidth',1.5, 'Color', [0 0.7 0.7]);
                          hold on;
                          plot(follower_xytheta_dsr(:,1), follower_xytheta_dsr(:,2), 'b-', 'LineWidth',1.5);
                          plot(xdsr, ydsr, 'ko', 'MarkerSize', 6)%, 'MarkerFaceColor','k');
                          plot(xdsr(1), ydsr(1), 'bo', 'MarkerSize', 6, 'MarkerFaceColor', [0 0.7 0.7])
                          plot(xdsr(end), ydsr(end), 'bo', 'MarkerSize', 6, 'MarkerFaceColor', [0 0 1])
%                           quiver(xdsr(1), ydsr(1), cos(thetadsr(1,3)'), sin(thetadsr(1,3)'), 0.5, 'Color',[0 0 1])
%                           quiver(xdsr(end), ydsr(end), cos(thetadsr(end,3)'), sin(thetadsr(end,3)'), 0.5, 'Color',[1 0 0])
%                           quiver(xdsr(2:end-1), ydsr(2:end-1), cos(thetadsr(2:end-1,3)'), sin(thetadsr(2:end-1,3)'), 0.5, 'Color',[0 0 0])
                          quiver(xdsr, ydsr, cos(thetadsr(:,3)'), sin(thetadsr(:,3)'), 0.5, 'Color',[0 0 0])
           
                          hold off
                          xlabel('X')
                          ylabel('Y')
                          xlim([0 6.5]);
                          ylim([-0.5 6]);
                          title('(a) Internal damping with DAC')
                          set(gca, 'FontSize', 16)

                          subplot(1,2,2);
                          plot(leader_xytheta_dsr2(:,1), leader_xytheta_dsr2(:,2), '-', 'LineWidth',1.5, 'Color', [0 0.7 0.7]);
                          hold on;
                          plot(follower_xytheta_dsr2(:,1), follower_xytheta_dsr2(:,2), 'b-', 'LineWidth',1.5);
                          plot(xdsr2, ydsr2, 'ko', 'MarkerSize', 6)%, 'MarkerFaceColor','k');
                          plot(xdsr2(1), ydsr2(1), 'bo', 'MarkerSize', 6, 'MarkerFaceColor', [0 0.7 0.7])
                          plot(xdsr2(end), ydsr2(end), 'bo', 'MarkerSize', 6, 'MarkerFaceColor', [0 0 1])
%                           quiver(xdsr(1), ydsr(1), cos(thetadsr(1,3)'), sin(thetadsr(1,3)'), 0.5, 'Color',[0 0 1])
%                           quiver(xdsr(end), ydsr(end), cos(thetadsr(end,3)'), sin(thetadsr(end,3)'), 0.5, 'Color',[1 0 0])
%                           quiver(xdsr(2:end-1), ydsr(2:end-1), cos(thetadsr(2:end-1,3)'), sin(thetadsr(2:end-1,3)'), 0.5, 'Color',[0 0 0])
                          quiver(xdsr2, ydsr2, cos(thetadsr2(:,3)'), sin(thetadsr2(:,3)'), 0.5, 'Color',[0 0 0])
                          
                          hold off
                          xlabel('X')
                          ylabel('Y')
                          xlim([0 6.5]);
                          ylim([-0.5 6]);
                          title('(b) Internal damping without DAC')
                          set(gca, 'FontSize', 16)

          
                          ff.Position = [100 100 540*2 400];
                          set(gcf,'color','w');

                          pause(0.001)
                           if (take_video == 1)
                               F = getframe(gcf);
                               writeVideo(vidObj,F);
                           end


                end

          
    
                        
            end

            if (take_video == 1)
%              close(vidObj0);
%              close(vidObj1);
%              close(vidObj2);
              close(vidObj);
            end



%             save('Oct_27_2Dsimuls_L3_withnoise.mat', "time_vec", "leader_xytheta_dsr", "follower_xytheta_dsr", ...
%                 "leader_xytheta_m1", "follower_xytheta_m1", ...
%                 "leader_xytheta_m2", "follower_xytheta_m2");

%             save('Oct_27_2Dsimuls_L3_withoutnoise_v2.mat', "time_vec", "leader_xytheta_dsr", "follower_xytheta_dsr", ...
%                 "leader_xytheta_m1", "follower_xytheta_m1", ...
%                 "leader_xytheta_m2", "follower_xytheta_m2");



%             save('Oct_27_2Dsimuls_L10_withnoise.mat', "time_vec", "leader_xytheta_dsr", "follower_xytheta_dsr", ...
%                 "leader_xytheta_m1", "follower_xytheta_m1", ...
%                 "leader_xytheta_m2", "follower_xytheta_m2");

%             save('Oct_27_2Dsimuls_L10_withoutnoise_v2.mat', "time_vec", "leader_xytheta_dsr", "follower_xytheta_dsr", ...
%                 "leader_xytheta_m1", "follower_xytheta_m1", ...
%                 "leader_xytheta_m2", "follower_xytheta_m2");


%             save('Oct_27_2Dsimuls_L200_withoutnoise.mat', "time_vec", "leader_xytheta_dsr", "follower_xytheta_dsr", ...
%                 "leader_xytheta_m1", "follower_xytheta_m1", ...
%                 "leader_xytheta_m2", "follower_xytheta_m2");

%             save('Oct_27_2Dsimuls_L200_withoutnoise_v2.mat', "time_vec", "leader_xytheta_dsr", "follower_xytheta_dsr", ...
%                 "leader_xytheta_m1", "follower_xytheta_m1", ...
%                 "leader_xytheta_m2", "follower_xytheta_m2");



end

%%
nfig=nfig+1; figure(nfig);
plot(time_vec, abs(follower_xytheta_dsr(:,3)-theta_desired), 'b:', ...
    time_vec, abs(follower_xytheta_dsr2(:,3)-theta_desired), 'r:', 'LineWidth', 3);
legend('With DAC', 'Without DAC');
ylabel('Follower (i=N) tracking error (e_N)');
xlabel('Time (s)')
set(gca, 'FontSize', 16)

nfig=nfig+1; figure(nfig);
plot(time_vec, theta_desired, 'g-', ...
    time_vec, follower_xytheta_dsr(:,3), 'b--', ...
    time_vec, follower_xytheta_dsr2(:,3), 'r:', 'LineWidth', 3);
legend('I_s', 'With Dynamic Consensus', 'Without Dynamic Consensus');
ylabel('Follower (i=N) orientation response (I_N)');
xlabel('Time (s)')
set(gca, 'FontSize', 16)


return

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vf = filtered_derivative(t, x, wf, n)
% FILTERED_dERIVATIVE  take derivative of a filtered signal.
%   vf = filtered_derivative(t, x, wf, n) filters x to xf and takes its
%   derivative vf

    A = -wf*eye(n); B = wf*eye(n);
    C = -wf*eye(n); D = wf*eye(n);
    
    sys = ss(A, B, C, D);
    
    vf = lsim(sys, x, t);

end



function [t, I, Idot] = second_order_simulation(t, Is, n, c_s, c_sd, cd)


        %number of discrete x points
            num_of_x = n; %equal to number of agents in 1d network

            %number of time steps 
            num_of_t = length(t);

            I = zeros(num_of_x, num_of_t);
            Idot = zeros(num_of_x, num_of_t);

            %adding the initial condition on left end with sin input
            Idsr(1,:) = Is;
            % I(2,:) = Is;


            for k=3:1:num_of_t

               for i=1:1:(num_of_x) 

                   if i==1 
                       Deltai_kminus1 =  Idsr(i,k-1) - 0.5*Idsr(i+1,k-1) -0.5*Is(k-1) ; 
                       Deltai_kminus2 =  Idsr(i,k-2) - 0.5*Idsr(i+1,k-2) -0.5*Is(k-2) ; 
                   elseif i == num_of_x        
                       Deltai_kminus1 = -1*Idsr(i-1,k-1) + Idsr(i,k-1); 
                       Deltai_kminus2 = -1*Idsr(i-1,k-2) + Idsr(i,k-2);         
                   else
                       Deltai_kminus1 = -0.5*Idsr(i-1,k-1) + Idsr(i,k-1) - 0.5*Idsr(i+1,k-1);  
                       Deltai_kminus2 = -0.5*Idsr(i-1,k-2) + Idsr(i,k-2) - 0.5*Idsr(i+1,k-2); 
                   end

        
                       Idsr(i,k) = Idsr(i,k-1) - gamma_dsr*delta_t*Deltai_kminus1 - beta1_dsr*(Deltai_kminus1-Deltai_kminus2) ...
                                + beta2_dsr*(Idsr(i,k-1)-Idsr(i,k-2));



               end

            end

end


%return