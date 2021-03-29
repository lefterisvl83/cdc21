clear all
close all
clc
%
%%%%%%%% ENTER INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tol = .01;                 % tolerance used below in the big for loop  %             
% uk0 = [20 20 30 40]';    % initial admission rates                   % 
uk0 = rand(30,1)*1;        % >> >>                                     %
% wk0 = [10 20 30 40]*1';  % initial queued requests                   % 
wk0 = rand(30,1)*10;       % >> >>                                     %
lambda = 2000;             % MEAN arrival rate                         %
% alfa = [1 2 3 4]'*5;     % alfa_i                                    %
alfa = rand(30,1)*1;       % >> >>                                     %
% veta = [1 1 1 1]'*.5;    % beta_i                                    %
veta = rand(30,1);         % >> >>                                     %
sim_time = 10000;          % continuous time (actual time)             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
n = length(veta); % number of nodes
uk = []; % admission
uk(:,1) = uk0;
wk = []; % queued requests
wk(:,1) = wk0;
Tk = []; % cycle period
t_event = 0; % event time instant (scalar changes every iteration)        
eventss = 1; % number of events
%
t_cont = 0:0.01:sim_time;
for i = 1:length(t_cont)
y(i) = exprnd(1/100); % random interarrival times (according to exponential distribution)
end
% %
yy = [];
yy(1) = 0;
for j = 1:length(y)
    yy(j+1) = yy(j) + y(j); % actual time at which an arrival occurs
end
%
admission = zeros(1,length(t_cont)); % total admissions at each simulation step
tot_admitted = 0;
%
for i = 1:length(t_cont)
    admission(i) = veta'*uk(:,eventss)*(yy(i) - t_event) + .5*sum(alfa)*...
        (yy(i) - t_event)*(yy(i) - t_event);
    if admission(i) >= (i - tot_admitted) && yy(i) - t_event > .2 % second condition is for smoothing Tk a little bit....
        tot_admitted = tot_admitted + admission(i);
        Tk(eventss) = yy(i) - t_event;
        t_event = yy(i);
        eventss = eventss + 1;
        uk(:,eventss) = diag(veta)*uk(:,eventss-1) + alfa*Tk(eventss-1);
    end
end
% 
% the remaining code is identical to the deterministic case, 
% refer to AIMD_admission_resource_allocationCDC21.m 
%
figure(1) % Maximum admission rates ui(k)
for i = 1:n
    plot(1:eventss, uk(i,1:eventss), 'LineWidth', 2)
    xlabel('Events')
    ylabel('$u_{i}(k)$','Interpreter','latex', 'FontSize', 18)
    title('Maximum admission rates (end of cycle)')
    leg_text{i}=['node-', num2str(i)]; 
    legend(leg_text)
    grid on
    hold on
end
%
figure(2) % Cycle Period T(k)
plot(1:eventss-1, Tk(1:eventss-1), 'LineWidth', 2)
xlabel('Events')
ylabel('$T(k)$', 'Interpreter', 'latex', 'FontSize', 18)
title('Cycle Period')
grid on
%
%
Tkcum = []; % actual time instants at which events occur
Tkcum(1) = Tk(1);
for i = 1:(eventss-2)
    Tkcum(i+1) = Tkcum(i) + Tk(i+1);
end
%
Tkcum = [0 Tkcum];
t_cont = 0:.01:Tkcum(end);
interv = length(t_cont);
u_cont = []; 
u_cont(:,1) = uk(:,1);
w_cont = zeros(n,interv);
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
            - gkamak(:,j-1)*(t_cont(i)-Tkcum(j-1))); % USE floor() for realistic queue length!!!
    else
        if j == 1
            wk(:,j) = wk0;
        else
            wk(:,j) = w_cont(:,i-1);
        end
        j = j + 1;
        gkamak(:,j-1) = diag(veta)*uk(:,j-1) + sqrt(2*alfa).*sqrt(wk(:,j-1));
         u_cont(:,i) = diag(veta)*uk(:,j-1) + alfa*(t_cont(i)-Tkcum(j-1));
        w_cont(:,i) = floor(wk(:,j-1) + diag(veta)*uk(:,j-1).*(t_cont(i)-Tkcum(j-1))...
           + .5*alfa*(t_cont(i)-Tkcum(j-1))*(t_cont(i)-Tkcum(j-1)) ...
            - gkamak(:,j-1)*(t_cont(i)-Tkcum(j-1))); % USE floor() for realistic queue length!!!
    end
end
%
% Calculation of queueing time: .5(wi(k) + wi(k+1) / ( ui(k) + .5*ai*T(k) )
Tqueueing = [];
for i = 1:(eventss-2)
    for j = 1:n
    Tqueueing(j,i) =  .5*(wk(j,i) + wk(j,i+1))/(uk(j,i) + .5*alfa(j)*Tk(i));
    end
end
% 
% Calculation of average admission rates
uk_av = [];
for i = 1:(eventss-1)
    uk_av(:,i) = .5*(diag(veta)*uk(:,i)+uk(:,i+1));
end
%
figure(3) % number of queued requests at node-i
for i = 1:n
    plot(t_cont, w_cont(i,:), 'LineWidth', 1.5)
    xlabel('Time [sec]')
    ylabel('$w_{i}(t)$','Interpreter','latex', 'FontSize', 18)
    title('Queued Requests')
    leg_text{i}=['node-', num2str(i)]; 
    legend(leg_text)
    grid on
    hold on
end
hold off
%
figure(4) % continuous-time admission rates
for i = 1:n
    plot(t_cont, u_cont(i,:), 'LineWidth', 1.5)
    xlabel('Time [sec]')
    ylabel('$u_{i}(t)$','Interpreter','latex', 'FontSize', 18)
    title('Admission rates')
    leg_text{i}=['node-', num2str(i)]; 
    legend(leg_text)
    grid on
    hold on
end
hold off
%
j = 1;
col = gray(n+1);
figure(5) % resource allocation gamma(k) (constant throughout cycles)
for i = 1:n
    plot(1:eventss, gkamak(i,1:eventss),'-', 'LineWidth', 1.5, 'color', col(i,:))
    legtxt{j} = ['node-', num2str(i)]; 
    hold on
    j = j + 1;
    plot( 1:(eventss-1), uk_av(i,1:(eventss-1)),'--', 'LineWidth', 1.5, 'color', col(i,:))
    legtxt{j} = ['node-', num2str(i)];
    j = j + 1;
    xlabel('Events')
    ylabel('$\gamma_{i}(k), u_{i}^{\mathbf{av}}(k)$','Interpreter','latex', 'FontSize', 18)
    title('Service (solid) and average admission (dashed) rates')
    legend(legtxt)
    grid on
end
hold off
%
figure(6) % Queueing time associated with each node
for i = 1:n
    plot(1:(eventss-2), Tqueueing(i,:), 'LineWidth', 1.5)
    xlabel('Events')
    ylabel('$T_{\mathbf{queueing}}(k)$','Interpreter','latex', 'FontSize', 18)
    title('Queueing/waiting time (delay)')
    leg_text{i}=['node-', num2str(i)]; 
    legend(leg_text)
    grid on
    hold on
end
hold off
%


    
    
    
    
    
    
    
    
    
    
    
    
    
    





