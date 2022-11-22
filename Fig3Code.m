clear all
clc
close all


load('f1_f2_DSR_m_data_v2.mat')  ;
nfig=0;
dt = 10^(-4);
            %shift in time from the first leaders peak             

%%
            nfig=nfig+1; figure(nfig);
            shift_t = 1/dt; %4 secs
subplot(1,3,1)
            xd_1 = 0.2; yd_1 = 0.998927;
            xd_2 = 49.4; yd_2 = 0.97872;  
%             plot(linspace(a, num_of_x*a, num_of_x), Idsr1(:,indx_DSR_L_peak_1), 'b-', ...
%                 linspace(a, num_of_x*a, num_of_x), Idsr2(:,indx_DSR_L_peak_1), 'k-', ... 
%                 linspace(a, num_of_x*a, num_of_x), Idsr1(:,indx_DSR_L_peak_1+fix(shift_t*1.14)), 'b--', ...
%                 linspace(a, num_of_x*a, num_of_x), Idsr2(:,indx_DSR_L_peak_1+fix(shift_t*1.14)), 'k--', 'LineWidth', 3);
            plot(linspace(a, num_of_x*a, num_of_x), Idsr2(:,indx_DSR_L_peak_2), 'b-', ...
                linspace(a, num_of_x*a, num_of_x), Idsr2(:,indx_DSR_L_peak_2+fix(shift_t*1.14)), 'b:', 'LineWidth', 3);
            hold on
            plot(xd_1, yd_1, 'bo', 'MarkerSize', 10, 'MarkerFaceColor','blue');
            plot(xd_2, yd_2, 'bo', 'MarkerSize', 10, 'LineWidth', 4);
            ylabel('Agent Responses (Internal damping)');
            xlabel('X-location of agents');
            grid on;
            set(gca, 'FontSize', 20);
            ylim([0 1.5])
            xlim([0 200])
            title('(a) 1 Hz', 'FontSize', 18)
subplot(1,3,2)
            xm_1 = 0.2; ym_1 = 0.993968; 
            xm_2 = 18.4; ym_2 = 0.422115;
            shift_t = 2/dt; %4 secs
%             plot(linspace(a, num_of_x*a, num_of_x), Im_1(:,indx_m_L_peak_1), 'r-', ...
%                 linspace(a, num_of_x*a, num_of_x), Im_2(:,indx_m_L_peak_1), 'k-', ... 
%                 linspace(a, num_of_x*a, num_of_x), Im_1(:,indx_m_L_peak_1+fix(shift_t*1.14)), 'r--', ...
%                 linspace(a, num_of_x*a, num_of_x), Im_2(:,indx_m_L_peak_1+fix(shift_t*1.14)), 'k--', 'LineWidth', 3);
           plot(linspace(a, num_of_x*a, num_of_x), Im_2(:,indx_m_L_peak_2), 'r-', ...
                linspace(a, num_of_x*a, num_of_x), Im_2(:,indx_m_L_peak_2+fix(shift_t*1.14)), 'r:', 'LineWidth', 3);
            hold on
            plot(xm_1, ym_1, 'ro', 'MarkerSize', 10, 'MarkerFaceColor','red');
            plot(xm_2, ym_2, 'ro', 'MarkerSize', 10, 'LineWidth', 4); 
           ylabel('Agent Responses (Viscous damping)');
            xlabel('X-location of agents');
            grid on;
            set(gca, 'FontSize', 20);
            ylim([0 1.5])
            xlim([0 200])
            title('(b) 1 Hz', 'FontSize', 18)

            indx = 0.0025/dt;
            shift_t = 0.07/dt;

%             nfig=nfig+1; figure(nfig);
subplot(1,3,3)
%             xm_1 = 1; ym_1 = 0.877572;
%             xm_2 = 5.4; ym_2 = 0.571359; 
% 
%             xd_1 = 0.2; yd_1 = 0.492328;
%             xd_2 = 5.6; yd_2 = 0.247916; 
%             plot(linspace(a, num_of_x*a, num_of_x), Im_3(:,indx+fix(shift_t/2)), 'r-', ...
%                 linspace(a, num_of_x*a, num_of_x), Im_3(:,indx+fix(shift_t*1.14)), 'r:', ...
%                 linspace(a, num_of_x*a, num_of_x), Idsr3(:,indx+fix(shift_t/2)), 'b-', ...
%                 linspace(a, num_of_x*a, num_of_x), Idsr3(:,indx+fix(shift_t*1.14)), 'b:', ...
%                 'LineWidth', 3);
            xm_1 = 0.2; ym_1 = 0.982412;
            xm_2 = 2.6; ym_2 = 0.755849; 

            xd_1 = 1.8; yd_1 = 0.148954; %xd_1 = 0.2; yd_1 = 0.262958;
            xd_2 = 2.8; yd_2 = 0.113442; 
              plot(linspace(a, num_of_x*a, num_of_x), Im_3(:,indx_m_L_peak_3), 'r-', ...
                linspace(a, num_of_x*a, num_of_x), Im_3(:,indx+fix(shift_t*1.14)), 'r:', ...
                linspace(a, num_of_x*a, num_of_x), Idsr3(:,indx_DSR_L_peak_3+250), 'b-', ...
                linspace(a, num_of_x*a, num_of_x), Idsr3(:,indx+fix(shift_t*1.14)), 'b:', ...
                'LineWidth', 3);
             hold on
            plot(xd_1, yd_1, 'bo', 'MarkerSize', 10, 'MarkerFaceColor','blue');
            plot(xd_2, yd_2, 'bo', 'MarkerSize', 10, 'LineWidth', 4);
             hold on
            plot(xm_1, ym_1, 'ro', 'MarkerSize', 10, 'MarkerFaceColor','red');
            plot(xm_2, ym_2, 'ro', 'MarkerSize', 10, 'LineWidth', 4);
            ylabel('Agent Responses');
            xlabel('X-location of agents');
             grid on;
            set(gca, 'FontSize', 20);
%             legend('Viscous damping (t=t_0)', 'Viscous damping (t=t_0+\Delta t)', ...
%                 'Internal damping (t=t_0)', 'Internal damping (t=t_0+\Delta t)', 'FontSize', 16)
             ylim([0 1.5])
            xlim([0 10])
            title('(c) 10 Hz', 'FontSize', 18)

            ff.Position = [100 100 540*1.2 1500];

saveas(gcf,'fig_threefreqs_plot','epsc')
