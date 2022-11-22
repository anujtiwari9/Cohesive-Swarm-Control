clear all 
clc
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nfig=0;
video_steps = 200; %200
take_video = 0; %1 to save video of individual cases
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%simulations 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L = 400;     %length of network in m
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

video_name = ['Oct_27_2Dsimuls_shortnetwork_withnoise_point1.mp4'];

if (L <= 17.55*a)
        %%%selecting stable beta1
    
    
    
    % % % 
%        beta2_dsr = (pi^2*a^2/(4*L^2))*( sqrt( gamma_dsr*dt/(2*D) + 4*(L^2)/(pi^2*a^2) + beta1_dsr/(2*D) ) - sqrt( gamma_dsr*dt/(2*D)) )^2
     c = sqrt(gamma_dsr*a^2/(2*D*delta_t*beta2_dsr));
    
       beta1_dsr = 4*gamma_dsr*L/(pi*c)

else    

    %%%selecting stable beta1
    beta1_dsr = 0.9*((beta2_dsr+1) - gamma_dsr*delta_t/2 )/max(lambda_A);

    
    
    
    % % % 
       beta2_dsr = (pi^2*a^2/(4*L^2))*( sqrt( gamma_dsr*dt/(2*D) + 4*(L^2)/(pi^2*a^2) + beta1_dsr/(2*D) ) - sqrt( gamma_dsr*dt/(2*D)) )^2
     c = sqrt(gamma_dsr*a^2/(2*D*delta_t*beta2_dsr));
    
    %   beta1_dsr = 4*gamma_dsr*L/(pi*c)

    video_name = ['Oct_27_2Dsimuls_longnetwork3.mp4'];



end

% return

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

tend = max(0.9*L/v, 40) %2.8835e+03 %20*L/v;           %s - duration of simulation time -- BEFORE REFLECTION 
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%define input source
input1 = sin(2*pi*ts/T2);
input2 = 0*sin(2*pi*ts/T3);
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
             vidObj = VideoWriter(video_name, 'MPEG-4'); open(vidObj);
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
                

               thetadsr = Idsr1*pi/2; %thetadsr + pi/2*Idsr1;
               xdsr = xdsr + delta_t*cos(thetadsr(:,3)');
               ydsr = ydsr + delta_t*sin(thetadsr(:,3)');

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


               thetam1 = Im_1*pi/2; %thetadsr + pi/2*Idsr1;
               xm1 = xm1 + delta_t*cos(thetam1(:,3)');
               ym1 = ym1 + delta_t*sin(thetam1(:,3)');


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

               thetam2 = Im_2*pi/2; %thetadsr + pi/2*Idsr1;
               xm2 = xm2 + delta_t*cos(thetam2(:,3)');
               ym2 = ym2 + delta_t*sin(thetam2(:,3)');



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
 
%                          ff = figure(100);
%                            subplot(1,3,1)
%                           plot(linspace(a, num_of_x*a, num_of_x), Im_2(:,iii-1), '.', 'MarkerSize',6);
%                           axis([min(1*a) max(num_of_x*a) -1 2]);
%                           grid on;
% %                           xlabel('X-location of Agents','fontSize',24);
%                           ylabel('Agent Responses','fontSize',20);              
%                           %titlestring = ['DSR, \beta_2=' num2str(beta2_dsr) ', \omega_s=' num2str(f) ' Hz, Time = ',num2str(t(k)), ' s']; %['TIME STEP = ',num2str(j), ' TIME = ',num2str(t(j)), 'second'];
%                           titlestring = ['Case 0: Undamped'];
%                           title(titlestring ,'fontsize',20);  
%                           set(gca, 'FontSize', 20); 
% %                           pause(0.001)
% %                           if (take_video == 1)
% %                             F = getframe(gcf);
% %                             writeVideo(vidObj0,F);
% %                           end
% 
% 
% %                            figure(200);
%                              subplot(1,3,2)
%                           plot(linspace(a, num_of_x*a, num_of_x), Im_1(:,iii-1), '.', 'MarkerSize',6);
%                           axis([min(1*a) max(num_of_x*a) -1 2]);
%                           grid on;
%                           xlabel('X-location of Agents','fontSize',20);
% %                           ylabel('Agent Responses','fontSize',24);              
%                           %titlestring = ['DSR, \beta_2=' num2str(beta2_dsr) ', \omega_s=' num2str(f) ' Hz, Time = ',num2str(t(k)), ' s']; %['TIME STEP = ',num2str(j), ' TIME = ',num2str(t(j)), 'second'];
%                           titlestring = ['Case 1: Viscous damping'];
%                           title(titlestring ,'fontsize',20);  
%                           set(gca, 'FontSize', 20);
% %                          pause(0.001)
% %                          if (take_video == 1)
% %                             F = getframe(gcf);
% %                             writeVideo(vidObj1,F);
% %                          end
%                             
% %                          figure(300);
%                           subplot(1,3,3)
%                           plot(linspace(a, num_of_x*a, num_of_x), Idsr1(:,iii-1), '.', 'MarkerSize',6);
%                           axis([min(1*a) max(num_of_x*a) -1 2]);
%                           grid on;
% %                           xlabel('X-location of Agents','fontSize',24);
% %                           ylabel('Agent Responses','fontSize',24);              
%                           %titlestring = ['DSR, \beta_2=' num2str(beta2_dsr) ', \omega_s=' num2str(f) ' Hz, Time = ',num2str(t(k)), ' s']; %['TIME STEP = ',num2str(j), ' TIME = ',num2str(t(j)), 'second'];
%                           titlestring = ['Case 2: Internal damping'];
%                           title(titlestring ,'fontsize',20);  
%                           set(gca, 'FontSize', 20);
%                           ff.Position = [100 100 540*2 400];
%                           pause(0.001)
% %                           if (take_video == 1)
% %                                F = getframe(gcf);
% %                                writeVideo(vidObj2,F);
% %                           end
% %                            if (take_video == 1)
% %                                F = getframe(gcf);
% %                                writeVideo(vidObj,F);
% %                            end
    



                          %%2D simul
                          leader_xytheta_dsr = [leader_xytheta_dsr; xdsr(1) ydsr(1) thetadsr(1,3)];
                          follower_xytheta_dsr = [follower_xytheta_dsr; xdsr(end) ydsr(end) thetadsr(end,3)];

                          leader_xytheta_m1 = [leader_xytheta_m1; xm1(1) ym1(1) thetam1(1,3)];
                          follower_xytheta_m1 = [follower_xytheta_m1; xm1(end) ym1(end) thetam1(end,3)];

                          leader_xytheta_m2 = [leader_xytheta_m2; xm2(1) ym2(1) thetam2(1,3)];
                          follower_xytheta_m2 = [follower_xytheta_m2; xm2(end) ym2(end) thetam2(end,3)];

                          time_vec = [time_vec t(k)];
                          
                          
%                           ff = figure(100);
%                           subplot(1,3,3);
%                           plot(leader_xytheta_dsr(:,1), leader_xytheta_dsr(:,2), 'b-', 'LineWidth',1);
%                           hold on;
%                           plot(follower_xytheta_dsr(:,1), follower_xytheta_dsr(:,2), 'b--', 'LineWidth',1);
%                           plot(xdsr, ydsr, 'ko', 'MarkerSize', 4, 'MarkerFaceColor','k');
% %                           quiver(xdsr(1), ydsr(1), cos(thetadsr(1,3)'), sin(thetadsr(1,3)'), 0.5, 'Color',[0 0 1])
% %                           quiver(xdsr(end), ydsr(end), cos(thetadsr(end,3)'), sin(thetadsr(end,3)'), 0.5, 'Color',[1 0 0])
% %                           quiver(xdsr(2:end-1), ydsr(2:end-1), cos(thetadsr(2:end-1,3)'), sin(thetadsr(2:end-1,3)'), 0.5, 'Color',[0 0 0])
%                           quiver(xdsr, ydsr, cos(thetadsr(:,3)'), sin(thetadsr(:,3)'), 0.5, 'Color',[0 0 0])
%                           hold off
%                           xlabel('X')
%                           ylabel('Y')
%                           xlim([0 10*L/3]);
%                           ylim([-1 5]);
%                           title('(c) Internal damping')
%                           set(gca, 'FontSize', 24)
% 
%                           subplot(1,3,2);
%                           plot(leader_xytheta_m1(:,1), leader_xytheta_m1(:,2), 'r-', 'LineWidth',1);
%                           hold on;
%                           plot(follower_xytheta_m1(:,1), follower_xytheta_m1(:,2), 'r--', 'LineWidth',1);
%                           plot(xm1, ym1, 'ko', 'MarkerSize', 4, 'MarkerFaceColor','k');
% %                           quiver(xdsr(1), ydsr(1), cos(thetadsr(1,3)'), sin(thetadsr(1,3)'), 0.5, 'Color',[0 0 1])
% %                           quiver(xdsr(end), ydsr(end), cos(thetadsr(end,3)'), sin(thetadsr(end,3)'), 0.5, 'Color',[1 0 0])
% %                           quiver(xdsr(2:end-1), ydsr(2:end-1), cos(thetadsr(2:end-1,3)'), sin(thetadsr(2:end-1,3)'), 0.5, 'Color',[0 0 0])
%                           quiver(xm1, ym1, cos(thetam1(:,3)'), sin(thetam1(:,3)'), 0.5, 'Color',[0 0 0])
%                           hold off
%                           xlabel('X')
%                           ylabel('Y')
%                           xlim([0 10*L/3]);
%                           ylim([-1 5]);
%                           title('(b) Viscous damping')
%                           set(gca, 'FontSize', 24)
% 
%                           subplot(1,3,1);
%                           plot(leader_xytheta_m2(:,1), leader_xytheta_m2(:,2), 'g-', 'LineWidth',1);
%                           hold on;
%                           plot(follower_xytheta_m2(:,1), follower_xytheta_m2(:,2), 'g--', 'LineWidth',1);
%                           plot(xm2, ym2, 'ko', 'MarkerSize', 4, 'MarkerFaceColor','k');
% %                           quiver(xdsr(1), ydsr(1), cos(thetadsr(1,3)'), sin(thetadsr(1,3)'), 0.5, 'Color',[0 0 1])
% %                           quiver(xdsr(end), ydsr(end), cos(thetadsr(end,3)'), sin(thetadsr(end,3)'), 0.5, 'Color',[1 0 0])
% %                           quiver(xdsr(2:end-1), ydsr(2:end-1), cos(thetadsr(2:end-1,3)'), sin(thetadsr(2:end-1,3)'), 0.5, 'Color',[0 0 0])
%                           quiver(xm2, ym2, cos(thetam2(:,3)'), sin(thetam2(:,3)'), 0.5, 'Color',[0 0 0])
%                           hold off
%                           xlabel('X')
%                           ylabel('Y')
%                           xlim([0 10*(L/3)]);
%                           ylim([-1 5]);
%                           title('(a) Undamped')
%                           set(gca, 'FontSize', 24)
%                           ff.Position = [100 100 540*2 400];
% 
% 
%                           pause(0.001)
%                            if (take_video == 1)
%                                F = getframe(gcf);
%                                writeVideo(vidObj,F);
%                            end


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


%             save('Oct_27_2Dsimuls_L400_withoutnoise.mat', "time_vec", "leader_xytheta_dsr", "follower_xytheta_dsr", ...
%                 "leader_xytheta_m1", "follower_xytheta_m1", ...
%                 "leader_xytheta_m2", "follower_xytheta_m2");

%             save('Oct_27_2Dsimuls_L400_withoutnoise_v2.mat', "time_vec", "leader_xytheta_dsr", "follower_xytheta_dsr", ...
%                 "leader_xytheta_m1", "follower_xytheta_m1", ...
%                 "leader_xytheta_m2", "follower_xytheta_m2");



end

% return

%%

nfig=nfig+1;
ff=figure(nfig);
subplot(1,2,1);
load('Oct_27_2Dsimuls_L3_withoutnoise.mat');
plot(time_vec, leader_xytheta_dsr(:,3), 'b-', time_vec, follower_xytheta_dsr(:,3), 'r--', 'LineWidth', 4);
xlabel('Time (s)');
ylabel('Information signal (I_i(t))');
% legend('Leader (i=1)', 'Follower (i=N)');
set(gca, 'FontSize', 24);
grid on;
title('(a) N= 15')
set(gca, 'FontSize', 24);
ylim([-1 3]);


% subplot(1,3,2);
% load('Oct_27_2Dsimuls_L200_withoutnoise.mat');
% plot(time_vec, leader_xytheta_dsr(:,3), 'b-', time_vec, follower_xytheta_dsr(:,3), 'r--', 'LineWidth', 4);
% hold on
% load('Oct_27_2Dsimuls_L200_withoutnoise_v2.mat');
% plot(time_vec, follower_xytheta_dsr(:,3), 'g:', 'LineWidth', 4);
% xlabel('Time (s)');
% ylabel('Information signal (I_i(t))');
% set(gca, 'FontSize', 24);
% % legend('Leader (i=1)', 'Follower (i=N, no viscous damping)', 'Follower (i=N, with viscous damping)');
% grid on;
% title('(b) N=1000')
% set(gca, 'FontSize', 24);
% ylim([-3 5]);



subplot(1,2,2);
load('Oct_27_2Dsimuls_L400_withoutnoise.mat');
plot(time_vec, leader_xytheta_dsr(:,3), 'b-', time_vec, follower_xytheta_dsr(:,3), 'r--', 'LineWidth', 4);
hold on
load('Oct_27_2Dsimuls_L400_withoutnoise_v2.mat');
plot(time_vec, follower_xytheta_dsr(:,3), 'g:', 'LineWidth', 4);
xlabel('Time (s)');
ylabel('Information signal (I_i(t))');
set(gca, 'FontSize', 24);
legend('Leader (i=1)', 'Follower (i=N, no viscous damping)', 'Follower (i=N, with viscous damping)');
grid on;
title('(b) N=2000')
set(gca, 'FontSize', 24);
ylim([-3 5]);

ff.Position = [100 100 540*2.3 400];
saveas(gcf,'Fig_impact_of_length_4','epsc')