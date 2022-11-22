clear all 
clc
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nfig=0;
video_steps = 250; %200
take_video = 1; %1 to save video of individual cases
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%simulations 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L = 400;     %length of network in m
khat = 1;   %stiffness between agents
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


a = 0.2; %0.0224; %sqrt(v^2*beta2*delta_t*2/gamma)
v = sqrt(gamma_dsr*a^2/(2*1*delta_t*beta2_dsr))        %wave velocity in m/s

tend = max(0.9*L/v) %2.8835e+03 %20*L/v;           %s - duration of simulation time -- BEFORE REFLECTION 
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

% beta1_dsr = 0*4/(max(lambda_A)*(gamma_dsr*delta_t+2));
beta1_dsr = 0.9*((beta2_dsr+1) - gamma_dsr*delta_t/2 )/max(lambda_A);


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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%define input source
input1 = sin(2*pi*ts/T2);
input2 = 1*sin(2*pi*ts/T3);
input = input1 + input2;
initialzeros = 10000; %2
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
ylim([-1 2])

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
            Im_1 = zeros(num_of_x,3);
            Im_2 = zeros(num_of_x, 3);
            Im_3 = zeros(num_of_x, 3);

%             thetadsr = pi/2*ones(num_of_x, 3); %initially the angle is pi/2
%             xdsr = linspace(a, num_of_x*a, num_of_x);
%             ydsr = 0*xdsr;
%             vel = L/4; 
%              leader_x_history_dsr = [];
%              leader_y_history_dsr = [];
%              follower_x_history_dsr = [];
%              follower_y_history_dsr = [];

            %defining variables to take snaps of wave at two instants
                snap1_taken = 0;
                snap2_taken = 0;
                snap3_taken = 0;
                snap4_taken = 0;
                snap5_taken = 0;
                snap6_taken = 0; 

        if (take_video == 1)
%             vidObj0 = VideoWriter('Oct_24_Case0.mp4', 'MPEG-4'); open(vidObj0);
%             vidObj1 = VideoWriter('Oct_24_Case1.mp4', 'MPEG-4'); open(vidObj1);
%             vidObj2 = VideoWriter('Oct_24_Case2.mp4', 'MPEG-4'); open(vidObj2);
            vidObj = VideoWriter('MovieS1.mp4', 'MPEG-4'); open(vidObj);
        end

            
            for k = [3:1:num_of_t] 
                t(k)/t(end)

                iii=3;

            %frequency f1
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
               for i=1:1:(num_of_x) 
                   
                   if i==1 
                    Deltai_kminus1 =  Idsr1(i,iii-1) - 0.5*Idsr1(i+1,iii-1) -0.5*Is(k-1) ; 
                    Deltai_kminus2 =  Idsr1(i,iii-2) - 0.5*Idsr1(i+1,iii-2) -0.5*Is(k-2) ; 
                   elseif i == num_of_x        
                    Deltai_kminus1 = -1*Idsr1(i-1,iii-1) + Idsr1(i,iii-1); 
                    Deltai_kminus2 = -1*Idsr1(i-1,iii-2) + Idsr1(i,iii-2);         
                   else
                    Deltai_kminus1 = -0.5*Idsr1(i-1,iii-1) + Idsr1(i,iii-1) - 0.5*Idsr1(i+1,iii-1);  
                    Deltai_kminus2 = -0.5*Idsr1(i-1,iii-2) + Idsr1(i,iii-2) - 0.5*Idsr1(i+1,iii-2); 
                   end
                
                   Idsr1(i,iii) = Idsr1(i,iii-1) - gamma_dsr*delta_t*Deltai_kminus1 - beta1_dsr*(Deltai_kminus1-Deltai_kminus2) ...
                                + beta2_dsr*(Idsr1(i,iii-1)-Idsr1(i,iii-2));% + gamma_dsr*delta_t*B(i)*Is(k-1); 
               end

%                for i=1:1:(num_of_x) 
%                    
%                    if i==1 
%                     Deltai_kminus1 =  Idsr1(i,iii-1) - 0.5*Idsr1(i+1,iii-1);% -0.5*Is(k-1) ; 
%                     Deltai_kminus2 =  Idsr1(i,iii-2) - 0.5*Idsr1(i+1,iii-2);% -0.5*Is(k-2) ; 
%                    elseif i == num_of_x        
%                     Deltai_kminus1 = -1*Idsr1(i-1,iii-1) + Idsr1(i,iii-1); 
%                     Deltai_kminus2 = -1*Idsr1(i-1,iii-2) + Idsr1(i,iii-2);         
%                    else
%                     Deltai_kminus1 = -0.5*Idsr1(i-1,iii-1) + Idsr1(i,iii-1) - 0.5*Idsr1(i+1,iii-1);  
%                     Deltai_kminus2 = -0.5*Idsr1(i-1,iii-2) + Idsr1(i,iii-2) - 0.5*Idsr1(i+1,iii-2); 
%                    end
%                 
%                    Idsr1(i,iii) = Idsr1(i,iii-1) - gamma_dsr*delta_t*Deltai_kminus1 - beta1_dsr*(Deltai_kminus1-Deltai_kminus2) ...
%                                 + beta2_dsr*(Idsr1(i,iii-1)-Idsr1(i,iii-2)) + gamma_dsr*delta_t*B(i)*Is(k-1); 
%                end
                

%                thetadsr = Idsr1*pi/2; %thetadsr + pi/2*Idsr1;
%                xdsr = xdsr + delta_t*cos(thetadsr(:,3)');
%                ydsr = ydsr + delta_t*sin(thetadsr(:,3)');

                Idsr = Idsr1;
                f=f2;

                Idsr1(:,1:2) = Idsr1(:,2:3);


                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
               for i=1:1:(num_of_x) 
                   
                   if i==1 
                    Deltai_kminus1 =  Im_1(i,iii-1) - 0.5*Im_1(i+1,iii-1) -0.5*Is(k-1) ; 
                    Deltai_kminus2 =  Im_1(i,iii-2) - 0.5*Im_1(i+1,iii-2) -0.5*Is(k-2) ; 
                   elseif i == num_of_x        
                    Deltai_kminus1 = -1*Im_1(i-1,iii-1) + Im_1(i,iii-1); 
                    Deltai_kminus2 = -1*Im_1(i-1,iii-2) + Im_1(i,iii-2);         
                   else
                    Deltai_kminus1 = -0.5*Im_1(i-1,iii-1) + Im_1(i,iii-1) - 0.5*Im_1(i+1,iii-1);  
                    Deltai_kminus2 = -0.5*Im_1(i-1,iii-2) + Im_1(i,iii-2) - 0.5*Im_1(i+1,iii-2); 
                   end
                
                   Im_1(i,iii) = Im_1(i,iii-1) - gamma_m*delta_t*Deltai_kminus1 - beta1_m*(Deltai_kminus1-Deltai_kminus2) ...
                                + beta2_m_1*(Im_1(i,iii-1)-Im_1(i,iii-2));% + gamma_m*delta_t*B(i)*Is(k-1); 
               end


             %                 Idsr = Idsr1;
                %                 f=f2;

                Im_1(:,1:2) = Im_1(:,2:3);



                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
               for i=1:1:(num_of_x) 
                   
                   if i==1 
                    Deltai_kminus1 =  Im_2(i,iii-1) - 0.5*Im_2(i+1,iii-1) -0.5*Is(k-1) ; 
                    Deltai_kminus2 =  Im_2(i,iii-2) - 0.5*Im_2(i+1,iii-2) -0.5*Is(k-2) ; 
                   elseif i == num_of_x        
                    Deltai_kminus1 = -1*Im_2(i-1,iii-1) + Im_2(i,iii-1); 
                    Deltai_kminus2 = -1*Im_2(i-1,iii-2) + Im_2(i,iii-2);         
                   else
                    Deltai_kminus1 = -0.5*Im_2(i-1,iii-1) + Im_2(i,iii-1) - 0.5*Im_2(i+1,iii-1);  
                    Deltai_kminus2 = -0.5*Im_2(i-1,iii-2) + Im_2(i,iii-2) - 0.5*Im_2(i+1,iii-2); 
                   end
                
                   Im_2(i,iii) = Im_2(i,iii-1) - gamma_m*delta_t*Deltai_kminus1 - beta1_m*(Deltai_kminus1-Deltai_kminus2) ...
                                + beta2_m_2*(Im_2(i,iii-1)-Im_2(i,iii-2));% + gamma_m*delta_t*B(i)*Is(k-1); 
               end


             %                 Idsr = Idsr1;
                %                 f=f2;

                Im_2(:,1:2) = Im_2(:,2:3);         





               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
               for i=1:1:(num_of_x) 
                   
                   if i==1 
                    Deltai_kminus1 =  Im_3(i,iii-1) - 0.5*Im_3(i+1,iii-1) -0.5*Is1(k-1) ; 
                    Deltai_kminus2 =  Im_3(i,iii-2) - 0.5*Im_3(i+1,iii-2) -0.5*Is1(k-2) ; 
                   elseif i == num_of_x        
                    Deltai_kminus1 = -1*Im_3(i-1,iii-1) + Im_3(i,iii-1); 
                    Deltai_kminus2 = -1*Im_3(i-1,iii-2) + Im_3(i,iii-2);         
                   else
                    Deltai_kminus1 = -0.5*Im_3(i-1,iii-1) + Im_3(i,iii-1) - 0.5*Im_3(i+1,iii-1);  
                    Deltai_kminus2 = -0.5*Im_3(i-1,iii-2) + Im_3(i,iii-2) - 0.5*Im_3(i+1,iii-2); 
                   end
                
                   Im_3(i,iii) = Im_3(i,iii-1) - gamma_m*delta_t*Deltai_kminus1 - beta1_m*(Deltai_kminus1-Deltai_kminus2) ...
                                + beta2_m_2*(Im_3(i,iii-1)-Im_3(i,iii-2));% + gamma_m*delta_t*B(i)*Is1(k-1); 
               end


             %                 Idsr = Idsr1;
                %                 f=f2;

                Im_3(:,1:2) = Im_3(:,2:3);  


            
                
            
                if (mod(k,video_steps) == 0)
 
                         ff = figure(100);
                           subplot(1,3,1)
                          plot(linspace(a, num_of_x*a, num_of_x), Im_2(:,iii-1), '.', 'MarkerSize',6);
                          axis([min(1*a) max(num_of_x*a) -1 2]);
                          grid on;
%                           xlabel('X-location of Agents','fontSize',24);
                          ylabel('Agent Responses','fontSize',20);              
                          %titlestring = ['DSR, \beta_2=' num2str(beta2_dsr) ', \omega_s=' num2str(f) ' Hz, Time = ',num2str(t(k)), ' s']; %['TIME STEP = ',num2str(j), ' TIME = ',num2str(t(j)), 'second'];
                          titlestring = ['Case 0: Undamped'];
                          title(titlestring ,'fontsize',20);  
                          set(gca, 'FontSize', 18); 
%                           pause(0.001)
%                           if (take_video == 1)
%                             F = getframe(gcf);
%                             writeVideo(vidObj0,F);
%                           end


%                            figure(200);
                             subplot(1,3,2)
                          plot(linspace(a, num_of_x*a, num_of_x), Im_1(:,iii-1), '.', 'MarkerSize',6);
                          axis([min(1*a) max(num_of_x*a) -1 2]);
                          grid on;
                          xlabel('X-location of Agents','fontSize',20);
%                           ylabel('Agent Responses','fontSize',24);              
                          %titlestring = ['DSR, \beta_2=' num2str(beta2_dsr) ', \omega_s=' num2str(f) ' Hz, Time = ',num2str(t(k)), ' s']; %['TIME STEP = ',num2str(j), ' TIME = ',num2str(t(j)), 'second'];
                          titlestring = ['Case 1: Viscous damping'];
                          title(titlestring ,'fontsize',20);  
                          set(gca, 'FontSize', 18);
%                          pause(0.001)
%                          if (take_video == 1)
%                             F = getframe(gcf);
%                             writeVideo(vidObj1,F);
%                          end
                            
%                          figure(300);
                          subplot(1,3,3)
                          plot(linspace(a, num_of_x*a, num_of_x), Idsr1(:,iii-1), '.', 'MarkerSize',6);
                          axis([min(1*a) max(num_of_x*a) -1 2]);
                          grid on;
%                           xlabel('X-location of Agents','fontSize',24);
%                           ylabel('Agent Responses','fontSize',24);              
                          %titlestring = ['DSR, \beta_2=' num2str(beta2_dsr) ', \omega_s=' num2str(f) ' Hz, Time = ',num2str(t(k)), ' s']; %['TIME STEP = ',num2str(j), ' TIME = ',num2str(t(j)), 'second'];
                          titlestring = ['Case 2: Internal damping'];
                          title(titlestring ,'fontsize',20);  
                          set(gca, 'FontSize', 18);
                          ff.Position = [100 100 540*2 400];
                          pause(0.001)
%                           if (take_video == 1)
%                                F = getframe(gcf);
%                                writeVideo(vidObj2,F);
%                           end
                           if (take_video == 1)
                               F = getframe(gcf);
                               writeVideo(vidObj,F);
                           end
    



%                           %%2D simul
%                           leader_x_history_dsr = [leader_x_history_dsr xdsr(1)];
%                           leader_y_history_dsr = [leader_y_history_dsr ydsr(1)];
%                           follower_x_history_dsr = [follower_x_history_dsr xdsr(end)];
%                           follower_y_history_dsr = [follower_y_history_dsr ydsr(end)];
%                           
%                           
%                           figure(100);
%                           plot(leader_x_history_dsr, leader_y_history_dsr, 'b-', 'LineWidth',1);
%                           hold on;
%                           plot(follower_x_history_dsr, follower_y_history_dsr, 'b--', 'LineWidth',1);
%                           plot(xdsr, ydsr, 'ko', 'MarkerSize', 4, 'MarkerFaceColor','k');
% %                          quiver(xdsr, ydsr, cos(thetadsr(:,3)'), sin(thetadsr(:,3)'), 'k')
%                           hold off
%                           xlabel('X')
%                           ylabel('Y')
% %                           xlim([0 10*L/3]);
% %                           ylim([-1 5*L/3]);
% 
%                           pause(0.001)


                end

           

              %if (t(k) >= predicted_settling_time/1.5 && snap6_taken == 0)
              if (t(k) >= 7 && snap6_taken == 0)
                 
                  k_6 = k;
            
                  snap6_taken = 1;
                    

                  signal_6_undamped = Im_2(:,iii-1);
                  signal_6_viscous = Im_1(:,iii-1);
                  signal_6_internal = Idsr1(:,iii-1);

                  signal_6_ideal = Im_3(:,iii-1);
            
              end

    
                        
            end

            if (take_video == 1)
%              close(vidObj0);
%              close(vidObj1);
%              close(vidObj2);
             close(vidObj);
            end




end

%%
nfig=nfig+1;
ff=figure(nfig);
subplot(2,3,1);
plot(t, Is1, 'k-', 'LineWidth', 2)
xlabel('Time (s)')
ylabel('I_{s,s}');
ylim([-1, 2]);
grid on;
title('(a)')
set(gca, 'FontSize', 20);
subplot(2,3,2);
plot(t, Is2, 'k-', 'LineWidth', 0.25)
xlabel('Time (s)')
ylabel('I_{s,l}');
ylim([-1, 2]);
grid on;
title('(b)')
set(gca, 'FontSize', 20);
subplot(2,3,3);
plot(t, Is, 'k-', 'LineWidth', 0.5)
xlabel('Time (s)')
ylabel('Source signal (I_s = I_{s,s}+I_{s,l})');
ylim([-1, 2]);
grid on;
title('(c)')
set(gca, 'FontSize', 20);
subplot(2,3,4);
plot(linspace(a, num_of_x*a, num_of_x), signal_6_undamped, 'g-', 'LineWidth', 2)
hold on
plot(linspace(a, num_of_x*a, num_of_x), signal_6_ideal, 'k:', 'LineWidth', 3)
xlabel('X-location')
ylabel('Agent responses');
legend('Undamped', 'Ideal propagation', 'Location','southeast')
ylim([-1, 2]);
grid on;
title('(d)')
set(gca, 'FontSize', 20);
subplot(2,3,5);
plot(linspace(a, num_of_x*a, num_of_x), signal_6_viscous, 'r-', 'LineWidth',2)
hold on
plot(linspace(a, num_of_x*a, num_of_x), signal_6_ideal, 'k:', 'LineWidth', 3)
xlabel('X-location')
ylabel('Agent responses');
legend('Viscous damping', 'Ideal propagation', 'Location','southeast')
ylim([-1, 2]);
grid on;
title('(e)')
set(gca, 'FontSize', 20);
subplot(2,3,6);
plot(linspace(a, num_of_x*a, num_of_x), signal_6_internal, 'b-', 'LineWidth', 2)
hold on
plot(linspace(a, num_of_x*a, num_of_x), signal_6_ideal, 'k:', 'LineWidth', 3)
xlabel('X-location')
ylabel('Agent responses');
legend('Internal damping', 'Ideal propagation', 'Location','southeast')
ylim([-1, 2]);
grid on;
title('(f)')
set(gca, 'FontSize', 20);
ff.Position = [100 100 540*2.2 650];

saveas(gcf,'fig_noiseeffects_v3','epsc')
