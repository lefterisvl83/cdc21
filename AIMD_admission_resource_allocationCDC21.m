clear all
close all
clc
%
%%%%%%%% ENTER INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tol = .01;                  % tolerance used below in the big for loop %
num_vm = 4;                 % number of nodes                          %                  
uk0 = [0 5 10 20]';       % initial admission rates                  %        
wk0 = [5 10 15 20]'*1.5;    % initial queued requests                  %        
lambda = 100;               % workload                                 %
alfa = [1 2 3 4]'*5;        % alfa_i                                   %
veta = [.5 .5 .5 .5]';      % beta_i                                   %
events = 25;                % number of events                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
uk = []; % admission
uk(:,1) = uk0;
wk = []; % queued requests
wk(:,1) = wk0;
Tk = []; % cycle period
%
for i=1:events
    Tk(i) = 2*(lambda - veta'*uk(:,i))/(sum(alfa));
    uk(:,i+1) = diag(veta)*uk(:,i) + alfa*Tk(i);
end
figure(1) % Maximum admission rates ui(k)
for i = 1:4
    plot(1:events, uk(i,1:events), 'LineWidth', 2)
    xlabel('Events','FontSize', 18)
    ylabel('$u_{i}(k)$','Interpreter','latex', 'FontSize', 18)
    title('Maximum admission rates (end of cycle)','FontSize', 18)
    leg_text{i}=['node-', num2str(i)]; 
    legend(leg_text)
    grid on
    hold on
end
figure(2) % Cycle Period T(k)
plot(1:events, Tk(1:events), 'LineWidth', 2)
xlabel('Events','FontSize', 18)
ylabel('$T(k)$', 'Interpreter', 'latex', 'FontSize', 18)
title('Cycle Period','FontSize', 18)
grid on
%
%
Tkcum = []; % actual time instants at which events occur
Tkcum(1) = Tk(1);
for i = 1:(events-1)
    Tkcum(i+1) = Tkcum(i) + Tk(i+1);
end
%
Tkcum = [0 Tkcum];
t_cont = 0:.01:Tkcum(end);
interv = length(t_cont);
u_cont = []; 
u_cont(:,1) = uk(:,1);
w_cont = zeros(4,interv);
w_cont(:,1) = wk(:,1);
gkamak = [];
gkamak(:,1) = diag(veta)*uk(:,1) + sqrt(2*alfa).*sqrt(wk(:,1)); %%%
j = 1;
%
for i = 1:interv
    if abs(t_cont(i) - Tkcum(j)) > tol
        u_cont(:,i) = diag(veta)*uk(:,j-1) + alfa*(t_cont(i)-Tkcum(j-1));
        w_cont(:,i) = floor(wk(:,j-1) + diag(veta)*uk(:,j-1).*(t_cont(i)-Tkcum(j-1))...
           + .5*alfa*(t_cont(i)-Tkcum(j-1))*(t_cont(i)-Tkcum(j-1)) ...
            - gkamak(:,j-1)*(t_cont(i)-Tkcum(j-1))); % USE ceil() for realistic!!!
    else
        if j == 1
            wk(:,j) = wk0;
        else
            wk(:,j) = w_cont(:,i-1);
        end
        j = j + 1;
        gkamak(:,j-1) = diag(veta)*uk(:,j-1) + sqrt(2*alfa).*sqrt(wk(:,j-1));
         u_cont(:,i) = diag(veta)*uk(:,j-1) + alfa*(t_cont(i)-Tkcum(j-1));
        w_cont(:,i) = (wk(:,j-1) + diag(veta)*uk(:,j-1).*(t_cont(i)-Tkcum(j-1))...
           + .5*alfa*(t_cont(i)-Tkcum(j-1))*(t_cont(i)-Tkcum(j-1)) ...
            - gkamak(:,j-1)*(t_cont(i)-Tkcum(j-1))); % USE ceil() for realistic!!!
    end
end
%
% Calculation of queueing time: .5(wi(k) + wi(k+1) / ( ui(k) + .5*ai*T(k) )
Tqueueing = [];
for i = 1:(events-1)
    for j = 1:4
    Tqueueing(j,i) =  .5*(wk(j,i) + wk(j,i+1))/(uk(j,i) + .5*alfa(j)*Tk(i));
    end
end
% 
% Calculation of average admission rates
uk_av = [];
for i = 1:(events-1)
    uk_av(:,i) = .5*(diag(veta)*uk(:,i)+uk(:,i+1));
end
%
figure(3) % number of queued requests at node-i
for i = 1:4
    plot(t_cont, w_cont(i,:), 'LineWidth', 1.5)
    xlabel('Time [sec]','FontSize', 18)
    ylabel('$w_{i}(t)$','Interpreter','latex', 'FontSize', 18)
    title('Queued Requests','FontSize', 18)
    leg_text{i}=['node-', num2str(i)]; 
    legend(leg_text)
    grid on
    hold on
end
hold off
%
figure(4) % continuous-time admission rates
for i = 1:4
    plot(t_cont, u_cont(i,:), 'LineWidth', 1.5)
    xlabel('Time [sec]','FontSize', 18)
    ylabel('$u_{i}(t)$','Interpreter','latex', 'FontSize', 18)
    title('Admission rates','FontSize', 18)
    leg_text{i}=['node-', num2str(i)]; 
    legend(leg_text)
    grid on
    hold on
end
hold off
%
j = 1;
C = {'b','r','g','m',[.5 .6 .7],[.8 .2 .6]}; % Cell array of colors.
figure(5) % resource allocation gamma(k) (constant throughout cycles)
for i = 1:4
    plot(1:events, gkamak(i,1:events),'-', 'LineWidth', 1.5, 'color', C{i})
    legtxt{j} = ['node-', num2str(i)]; 
    hold on
    j = j + 1;
    plot( 1:(events-1), uk_av(i,1:(events-1)),'--', 'LineWidth', 1.5, 'color', C{i})
    legtxt{j} = ['node-', num2str(i)];
    j = j + 1;
    xlabel('Events','FontSize', 18)
    ylabel('$\gamma_{i}(k), u_{i}^{\mathbf{av}}(k)$','Interpreter','latex', 'FontSize', 18)
    title('Service and average admission rates','FontSize', 18)
%     leg_text{i}=['node-', num2str(i)]; 
    legend(legtxt)
    grid on
end
hold off
%
figure(6) % Queueing time associated with each node
for i = 1:4
    plot(1:(events-1), Tqueueing(i,:), 'LineWidth', 1.5)
    xlabel('Events','FontSize', 18)
    ylabel('$T_{i}(k)$','Interpreter','latex', 'FontSize', 18)
    title('Queueing time','FontSize', 18)
    leg_text{i}=['node-', num2str(i)]; 
    legend(leg_text)
    grid on
    hold on
end
hold off



    
    
    
    
    
    
    
    
    
    
    
    
    
    

