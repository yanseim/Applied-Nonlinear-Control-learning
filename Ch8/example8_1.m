% Applied Nonlinear Control. P317. example 8.1:MRAC control of an unknown
% mass. 20210731
clc
close all
%% parameters
m = 2;
lambda1 = 10;lambda2 = 25;lambda = 6;
gamma = 0.5;
dt = 0.001;
t = 0:dt:3;

func_xmdot2 = @(xmdot2, r)[xmdot2(2); -lambda1*xmdot2(2)-lambda2*xmdot2(1)+lambda2*r];% xmdot2 = [xm xmdot]
func_xdot2 = @(xdot2, u)[xdot2(2); u/m];% xdot2 = [x xdot]

%% (1) r=0
m_hat = 0;
% initial position for r=0
xdot2 = [0.5;0];
xmdot2 = [0.5;0];
log_r = []; log_u = []; log_x = []; log_xm = []; log_mhat = [];

for i = 1:length(t)
    % x_m
    r = 0;
    xmddot2 = func_xmdot2(xmdot2,r);
    xmdot2 = xmdot2 + xmddot2 * dt;
    % u
    xtilda = xdot2(1)-xmdot2(1);
    xtildadot = xdot2(2)-xmdot2(2);
    nu = xmddot2(2)-2*lambda*xtildadot-lambda^2*xtilda; %这里第一项容易犯错，是二阶导！如果写成一阶导，mhat就会一直上升
    u = m_hat*nu;
    % actual state
    xdot2 = xdot2 + func_xdot2(xdot2,u) * dt;
    % update m_hat
    s = xtildadot + lambda* xtilda;
    m_hat = m_hat - gamma * nu * s * dt;
    % log
    log_x = [log_x, xdot2];
    log_xm = [log_xm, xmdot2];
    log_u = [log_u, u];
    log_r = [log_r r];
    log_mhat = [log_mhat, m_hat]; 
end

%% plot r=0
title_text = 'results for $r=0$';
figure(1)
plot(t,log_x(1,:), 'LineWidth',2)
hold on
plot(t,log_xm(1,:),'--','LineWidth',2)
ylabel('tracking performance')
xlabel('time(sec)')
legend('actual x','reference model')
handle = title(title_text);
set(handle,'Interpreter','latex','FontSize',12);

figure(2)
plot(t, log_mhat, 'LineWidth',2);
hold on
plot(t, m*ones(1,length(t)),'--','LineWidth',2);
ylabel('parameter estimation');
xlabel('time(sec)')
legend('estimated','actual')
handle = title(title_text);
set(handle,'Interpreter','latex','FontSize',12);

%% (2) r = sin(4t), initial positions also change
m_hat = 0;
% intial position for r=sin(4t)
xdot2 = [0;0];
xmdot2 = [0;0];
log_r = []; log_u = []; log_x = []; log_xm = []; log_mhat = [];

for i = 1:length(t)
    % x_m
    r = sin(4*t(i));
    xmddot2 = func_xmdot2(xmdot2,r);
    xmdot2 = xmdot2 + xmddot2 * dt;
    % u
    xtilda = xdot2(1)-xmdot2(1);
    xtildadot = xdot2(2)-xmdot2(2);
    nu = xmddot2(2)-2*lambda*xtildadot-lambda^2*xtilda; %这里第一项容易犯错，是二阶导！如果写成一阶导，mhat就会一直上升
    u = m_hat*nu;
    % actual state
    xdot2 = xdot2 + func_xdot2(xdot2,u) * dt;
    % update m_hat
    s = xtildadot + lambda* xtilda;
    m_hat = m_hat - gamma * nu * s * dt;
    % log
    log_x = [log_x, xdot2];
    log_xm = [log_xm, xmdot2];
    log_u = [log_u, u];
    log_r = [log_r r];
    log_mhat = [log_mhat, m_hat]; 
end

%% plot r=sin(4t)
title_text = 'results for $r=sin(4t)$';
figure(3)
plot(t,log_x(1,:), 'LineWidth',2)
hold on
plot(t,log_xm(1,:),'--','LineWidth',2)
ylabel('tracking performance')
xlabel('time(sec)')
legend('actual x','reference model')
handle = title(title_text);
set(handle,'Interpreter','latex','FontSize',12);

figure(4)
plot(t, log_mhat, 'LineWidth',2);
hold on
plot(t, m*ones(1,length(t)),'--','LineWidth',2);
ylabel('parameter estimation');
xlabel('time(sec)')
legend('estimated','actual')
handle = title(title_text);
set(handle,'Interpreter','latex','FontSize',12);