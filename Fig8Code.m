clear all 
clc
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nfig=0; 
video_steps = 50; %200
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
    
       beta1_dsr = 0.9*((beta2_dsr+1) - gamma_dsr*delta_t/2 )/max(lambda_A)

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
input2 = 0.5*sin(2*pi*ts/T3);
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
             follower_xytheta_m2 = follower_xytheta_dsr
            

             %detect formation breaking
             formationbroken_dsr = 0;
             formationbroken_m1 = 0;
             formationbroken_m2 = 0;
             



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

               bigA = get_topological_A_matrix(xdsr, ydsr, threshold_dist, n);

               Delta_kminus1 = bigA*[Idsr1(:,iii-1); Is(k-1)];
               Delta_kminus2 = bigA*[Idsr1(:,iii-2); Is(k-2)];

               Idsr1(:,iii) = Idsr1(:,iii-1) - gamma_dsr*delta_t*Delta_kminus1 - beta1_dsr*(Delta_kminus1-Delta_kminus2) ...
                                    + beta2_dsr*(Idsr1(:,iii-1)-Idsr1(:,iii-2));% + gamma_dsr*delta_t*B(i)*Is(k-1); 

%                for i=1:1:(num_of_x) 
% 
% 
%                    if i==1 
%                     Deltai_kminus1 =  Idsr1(i,iii-1) - 0.5*Idsr1(i+1,iii-1) -0.5*Is(k-1) ; 
%                     Deltai_kminus2 =  Idsr1(i,iii-2) - 0.5*Idsr1(i+1,iii-2) -0.5*Is(k-2) ; 
%                    elseif i == num_of_x        
%                     Deltai_kminus1 = -1*Idsr1(i-1,iii-1) + Idsr1(i,iii-1);   
%                     Deltai_kminus2 = -1*Idsr1(i-1,iii-2) + Idsr1(i,iii-2);         
%                    else
%                     Deltai_kminus1 = -0.5*Idsr1(i-1,iii-1) + Idsr1(i,iii-1) - 0.5*Idsr1(i+1,iii-1);
%                     Deltai_kminus2 = -0.5*Idsr1(i-1,iii-2) + Idsr1(i,iii-2) - 0.5*Idsr1(i+1,iii-2); 
%                    end  
% 
%                    
%                    xdiff = xdsr(2:end) - xdsr(1:end-1); ydiff = ydsr(2:end)-ydsr(1:end-1);
%                    abs_diff = sqrt(xdiff.^2 + ydiff.^2);
% 
%                    if (max(abs_diff) > threshold_dist)
%                                     formationbroken_dsr = 1;
%                    end  
% 
%                    if (formationbroken_dsr == 0)
%                        Idsr1(i,iii) = Idsr1(i,iii-1) - gamma_dsr*delta_t*Deltai_kminus1 - beta1_dsr*(Deltai_kminus1-Deltai_kminus2) ...
%                                     + beta2_dsr*(Idsr1(i,iii-1)-Idsr1(i,iii-2));% + gamma_dsr*delta_t*B(i)*Is(k-1); 
%                    else
% 
%                         Idsr1(i,iii) = Idsr1(i,iii-1);
%                    end
% 
%                end
                
%                if (formationbroken_dsr == 0)
                       thetadsr = Idsr1*pi/2; %thetadsr + pi/2*Idsr1;
                       xdsr = xdsr + vel*delta_t*cos(thetadsr(:,3)');
                       ydsr = ydsr + vel*delta_t*sin(thetadsr(:,3)');
%                else 
%                        thetadsr = Idsr1*pi/2; %thetadsr + pi/2*Idsr1;
%                        xdsr = xdsr;
%                        ydsr = ydsr;
%                end



                Idsr = Idsr1;
                f=f2;

                Idsr1(:,1:2) = Idsr1(:,2:3);


                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
               bigA = get_topological_A_matrix(xm1, ym1, threshold_dist, n);

               Delta_kminus1 = bigA*[Im_1(:,iii-1); Is(k-1)];
               Delta_kminus2 = bigA*[Im_1(:,iii-2); Is(k-2)];
                
               Im_1(:,iii) = Im_1(:,iii-1) - gamma_m*delta_t*Delta_kminus1 - beta1_m*(Delta_kminus1-Delta_kminus2) ...
                                + beta2_m_1*(Im_1(:,iii-1)-Im_1(:,iii-2));


%                for i=1:1:(num_of_x) 
%                    
%                    if i==1 
%                     Deltai_kminus1 =  Im_1(i,iii-1) - 0.5*Im_1(i+1,iii-1) -0.5*Is(k-1) ; 
%                     Deltai_kminus2 =  Im_1(i,iii-2) - 0.5*Im_1(i+1,iii-2) -0.5*Is(k-2) ; 
%                    elseif i == num_of_x        
%                     Deltai_kminus1 = -1*Im_1(i-1,iii-1) + Im_1(i,iii-1); 
%                     Deltai_kminus2 = -1*Im_1(i-1,iii-2) + Im_1(i,iii-2);         
%                    else
%                     Deltai_kminus1 = -0.5*Im_1(i-1,iii-1) + Im_1(i,iii-1) - 0.5*Im_1(i+1,iii-1);  
%                     Deltai_kminus2 = -0.5*Im_1(i-1,iii-2) + Im_1(i,iii-2) - 0.5*Im_1(i+1,iii-2); 
%                    end
% 
%                     xdiff = xm1(2:end) - xm1(1:end-1); ydiff = ym1(2:end)-ym1(1:end-1);
%                     abs_diff = sqrt(xdiff.^2 + ydiff.^2);
% 
%                    if (max(abs_diff) > threshold_dist)
%                                     formationbroken_m1 = 1;
%                                     
%                    end  
%                 
%                    if (formationbroken_m1 == 0)
%                         Im_1(i,iii) = Im_1(i,iii-1) - gamma_m*delta_t*Deltai_kminus1 - beta1_m*(Deltai_kminus1-Deltai_kminus2) ...
%                                 + beta2_m_1*(Im_1(i,iii-1)-Im_1(i,iii-2));% + gamma_m*delta_t*B(i)*Is(k-1); 
%                    else
%                         Im_1(i,iii) = Im_1(i,iii-1);
% 
%                    end
%                end


%                if (formationbroken_m1 == 0)

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

%                else
%                    thetam1 = Im_1*pi/2; %thetadsr + pi/2*Idsr1;
%                    xm1 = xm1;% + delta_t*cos(thetam1(:,3)');
%                    ym1 = ym1;% + delta_t*sin(thetam1(:,3)');
%                end


             %                 Idsr = Idsr1;
                %                 f=f2;

                Im_1(:,1:2) = Im_1(:,2:3);



               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
               bigA = get_topological_A_matrix(xm2, ym2, threshold_dist, n);

               Delta_kminus1 = bigA*[Im_2(:,iii-1); Is(k-1)];
               Delta_kminus2 = bigA*[Im_2(:,iii-2); Is(k-2)];
                
            Im_2(:,iii) = Im_2(:,iii-1) - gamma_m*delta_t*Delta_kminus1 - beta1_m*(Delta_kminus1-Delta_kminus2) ...
                    + beta2_m_2*(Im_2(:,iii-1)-Im_2(:,iii-2));% + gamma_m*delta_t*B(i)*Is(k-1); 

%                for i=1:1:(num_of_x) 
%                    
%                    if i==1 
%                     Deltai_kminus1 =  Im_2(i,iii-1) - 0.5*Im_2(i+1,iii-1) -0.5*Is(k-1) ; 
%                     Deltai_kminus2 =  Im_2(i,iii-2) - 0.5*Im_2(i+1,iii-2) -0.5*Is(k-2) ; 
%                    elseif i == num_of_x        
%                     Deltai_kminus1 = -1*Im_2(i-1,iii-1) + Im_2(i,iii-1);
%                     Deltai_kminus2 = -1*Im_2(i-1,iii-2) + Im_2(i,iii-2);         
%                    else
%                     Deltai_kminus1 = -0.5*Im_2(i-1,iii-1) + Im_2(i,iii-1) - 0.5*Im_2(i+1,iii-1);  
%                     Deltai_kminus2 = -0.5*Im_2(i-1,iii-2) + Im_2(i,iii-2) - 0.5*Im_2(i+1,iii-2); 
%                    end
% 
%                      xdiff = xm2(2:end) - xm2(1:end-1); ydiff = ym2(2:end)-ym2(1:end-1);
%                     abs_diff = sqrt(xdiff.^2 + ydiff.^2);
% 
%                    if (max(abs_diff) > threshold_dist)
%                                     formationbroken_m2 = 1;
%                                     
%                    end  
% 
%                    
%                   if (formationbroken_m2 == 0)
%                         Im_2(i,iii) = Im_2(i,iii-1) - gamma_m*delta_t*Deltai_kminus1 - beta1_m*(Deltai_kminus1-Deltai_kminus2) ...
%                                 + beta2_m_2*(Im_2(i,iii-1)-Im_2(i,iii-2));% + gamma_m*delta_t*B(i)*Is(k-1); 
%                   else
%                         Im_2(i,iii) = Im_2(i,iii-1);
%                   end
%                end

%                if (formationbroken_m2 == 0)

                  
                 
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
% 
%                else
% 
%                    thetam2 = Im_2*pi/2; %thetadsr + pi/2*Idsr1;
%                    xm2 = xm2;% + delta_t*cos(thetam1(:,3)');
%                    ym2 = ym2;% + delta_t*sin(thetam1(:,3)');
% 
%                end



             %                 Idsr = Idsr1;
                %                 f=f2;

                Im_2(:,1:2) = Im_2(:,2:3);         





            
                
            
                if (mod(k,video_steps) == 0)
 

    



                          %%2D simul
                          leader_xytheta_dsr = [leader_xytheta_dsr; xdsr(1) ydsr(1) thetadsr(1,3)];
                          follower_xytheta_dsr = [follower_xytheta_dsr; xdsr(end) ydsr(end) thetadsr(end,3)];

                          leader_xytheta_m1 = [leader_xytheta_m1; xm1(1) ym1(1) thetam1(1,3)];
                          follower_xytheta_m1 = [follower_xytheta_m1; xm1(end) ym1(end) thetam1(end,3)];

                          leader_xytheta_m2 = [leader_xytheta_m2; xm2(1) ym2(1) thetam2(1,3)];
                          follower_xytheta_m2 = [follower_xytheta_m2; xm2(end) ym2(end) thetam2(end,3)];

                          time_vec = [time_vec t(k)];
                          
                          
                          ff = figure(100);
                          subplot(1,3,3);
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
                          if (formationbroken_dsr==1)
                              text(1,4,'Formation broken', 'FontSize', 20, 'Color',[1 0 0]);             
                          end
                          hold off
                          xlabel('X')
                          ylabel('Y')
                          xlim([0 6.5]);
                          ylim([-0.5 6]);
                          title('(c) Case 2: Internal damping')
                          set(gca, 'FontSize', 16)

                          subplot(1,3,2);
                          plot(leader_xytheta_m1(:,1), leader_xytheta_m1(:,2), '-', 'LineWidth',1.5, 'Color', [0 0.7 0.7]);
                          hold on;
                          plot(follower_xytheta_m1(:,1), follower_xytheta_m1(:,2), 'r-', 'LineWidth',1.5);
                          plot(xm1, ym1, 'ko', 'MarkerSize', 6)%, 'MarkerFaceColor','k');
                          plot(xm1(1), ym1(1), 'bo', 'MarkerSize', 6, 'MarkerFaceColor', [0 0.7 0.7])
                          plot(xm1(end), ym1(end), 'ko', 'MarkerSize', 6, 'MarkerFaceColor', [1 0 0])
%                           quiver(xdsr(1), ydsr(1), cos(thetadsr(1,3)'), sin(thetadsr(1,3)'), 0.5, 'Color',[0 0 1])
%                           quiver(xdsr(end), ydsr(end), cos(thetadsr(end,3)'), sin(thetadsr(end,3)'), 0.5, 'Color',[1 0 0])
%                           quiver(xdsr(2:end-1), ydsr(2:end-1), cos(thetadsr(2:end-1,3)'), sin(thetadsr(2:end-1,3)'), 0.5, 'Color',[0 0 0])
                          quiver(xm1, ym1, cos(thetam1(:,3)'), sin(thetam1(:,3)'), 0.5, 'Color',[0 0 0])
                          if (formationbroken_m1==1)
                              text(1,4,'Formation broken', 'FontSize', 20, 'Color',[1 0 0]);        
                          end
                          hold off
                          xlabel('X')
                          ylabel('Y')
                          xlim([0 6.5]);
                          ylim([-0.5 6]);
                          title('(b) Case 1: Viscous damping')
                          set(gca, 'FontSize', 16)

                          subplot(1,3,1);
                          plot(leader_xytheta_m2(:,1), leader_xytheta_m2(:,2), '-', 'LineWidth',1.5, 'Color', [0 0.7 0.7]);
                          hold on;
                          plot(follower_xytheta_m2(:,1), follower_xytheta_m2(:,2), 'g-', 'LineWidth',1.5);
                          plot(xm2, ym2, 'ko', 'MarkerSize', 6)%, 'MarkerFaceColor','k');
                          plot(xm2(1), ym2(1), 'bo', 'MarkerSize', 6, 'MarkerFaceColor', [0 0.7 0.7])
                          plot(xm2(end), ym2(end), 'ko', 'MarkerSize', 6, 'MarkerFaceColor', [0 1 0])
%                           quiver(xdsr(1), ydsr(1), cos(thetadsr(1,3)'), sin(thetadsr(1,3)'), 0.5, 'Color',[0 0 1])
%                           quiver(xdsr(end), ydsr(end), cos(thetadsr(end,3)'), sin(thetadsr(end,3)'), 0.5, 'Color',[1 0 0])
%                           quiver(xdsr(2:end-1), ydsr(2:end-1), cos(thetadsr(2:end-1,3)'), sin(thetadsr(2:end-1,3)'), 0.5, 'Color',[0 0 0])
                          quiver(xm2, ym2, cos(thetam2(:,3)'), sin(thetam2(:,3)'), 0.5, 'Color',[0 0 0])
                          if (formationbroken_m2==1)
                              text(1,4,'Formation broken', 'FontSize', 20, 'Color',[1 0 0]);            
                          end
                          hold off
                          xlabel('X')
                          ylabel('Y')
                          xlim([0 6.5]);
                          ylim([-0.5 6]);
                          title('(a) Case 0: Without damping')
                          set(gca, 'FontSize', 16)
                          ff.Position = [100 100 540*2 400];
                          set(gcf,'color','w');

                          pause(0.001)
                           if (take_video == 1)
                               F = getframe(gcf);
                               writeVideo(vidObj,F);
                           end


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

saveas(gcf,'Fig_15agents_topological','epsc')

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