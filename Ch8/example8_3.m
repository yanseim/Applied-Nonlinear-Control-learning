% Applied Nonlinear Control. P329. example 8.3: a first order plant
% 20210731
clc
close all
%%
dt = 0.01;
t = 0:dt:10;

gamma = 2;

a_p = -1; b_p = 3;
a_m = 4; b_m = 4;

func_ydot = @(y, u)-a_p*y+b_p*u;
func_xmdot = @(x_m, r)-a_m*x_m+b_m*r;
%% r=4
ar_hat = 0; ay_hat = 0;% intial values of both parameters of the controllers are 0
y = 0; x_m = 0;% initial conditions of the plant and the model are both zero
log_r = []; log_u = []; log_y = []; log_xm = []; log_ahat = []; log_bhat = [];
for i = 1:length(t)
    % x_m
    r = 4;
%     r = 4*sin(3*t(i));
    x_m = x_m + func_xmdot(x_m, r)* dt;
    % u
    u = ar_hat * r+ ay_hat * y;
    % actual state
    y = y + func_ydot(y, u) * dt;
    % update ar_hat and ay_hat
    e = y-x_m;
    ar_hat = ar_hat - gamma * e * r * dt;
    ay_hat = ay_hat - gamma * e * y * dt;
    % log
    log_y = [log_y, y];
    log_xm = [log_xm, x_m];
    log_u = [log_u, u];
    log_r = [log_r r];
    log_ahat = [log_ahat, ar_hat]; 
    log_bhat = [log_bhat, ay_hat]; 
end

% plot
figure(1)
plot(t,log_y(1,:), 'LineWidth',2)
hold on
plot(t,log_xm(1,:),'--','LineWidth',2)
ylabel('tracking performance')
xlabel('time(sec)')
legend('actual x','reference model')

figure(2)
plot(t, log_ahat, 'LineWidth',2);
hold on
plot(t, log_bhat, 'LineWidth',2);
plot(t, b_m/b_p*ones(1,length(t)),'--','LineWidth',2);% 这个真值不是简单的a_p b_p，是要计算过的
plot(t, (a_p-a_m)/b_p*ones(1,length(t)),'--','LineWidth',2);
ylabel('parameter estimation');
xlabel('time(sec)')
legend('estimated ar','estimated ay','actual ar','actual ay')

%% r = 4sin(3t)

ar_hat = 0; ay_hat = 0;% intial values of both parameters of the controllers are 0
y = 0; x_m = 0;% initial conditions of the plant and the model are both zero
log_r = []; log_u = []; log_y = []; log_xm = []; log_ahat = []; log_bhat = [];
for i = 1:length(t)
    % x_m
%     r = 4;
    r = 4*sin(3*t(i));
    x_m = x_m + func_xmdot(x_m, r)* dt;
    % u
    u = ar_hat * r+ ay_hat * y;
    % actual state
    y = y + func_ydot(y, u) * dt;
    % update ar_hat and ay_hat
    e = y-x_m;
    ar_hat = ar_hat - gamma * e * r * dt;
    ay_hat = ay_hat - gamma * e * y * dt;
    % log
    log_y = [log_y, y];
    log_xm = [log_xm, x_m];
    log_u = [log_u, u];
    log_r = [log_r r];
    log_ahat = [log_ahat, ar_hat]; 
    log_bhat = [log_bhat, ay_hat]; 
end

% plot
figure(3)
plot(t,log_y(1,:), 'LineWidth',2)
hold on
plot(t,log_xm(1,:),'--','LineWidth',2)
ylabel('tracking performance')
xlabel('time(sec)')
legend('actual x','reference model')

figure(4)
plot(t, log_ahat, 'LineWidth',2);
hold on
plot(t, log_bhat, 'LineWidth',2);
plot(t, b_m/b_p*ones(1,length(t)),'--','LineWidth',2);
plot(t, (a_p-a_m)/b_p*ones(1,length(t)),'--','LineWidth',2);
ylabel('parameter estimation');
xlabel('time(sec)')
legend('estimated ar','estimated ay','actual ar','actual ay')