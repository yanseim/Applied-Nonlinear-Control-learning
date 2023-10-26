% example 9.1
clc
close all
%% parameters
m1 = 1;
l1 = 1;
me = 2;
deltae = pi/6;
I1 = 0.12;
lc1 = 0.5;
Ie = 0.25;
lce = 0.6;

global a1 a2 a3 a4
a1 = I1+m1*lc1^2+Ie+me*lce^2+me*l1^2;
a2 = Ie+me*lce^2;
a3 = me*l1*lce*cos(deltae);
a4 = me*l1*lce*sin(deltae);

Kd = 100*eye(2);
Kp = 2000*eye(2);

log_qddot = [];log_qdot = [];log_q = [];log_tau = [];

%% initial
q = [0;0];qdot = [0;0];
qd = [pi/3;pi/2];

%% main loop
dt = 0.001;T = 1;
i = 1;
for t = 0:dt:T
    tau = -Kp*(q-qd)-Kd*qdot;
    qddot = dynamics(tau,q,qdot);
    qdot = qdot+qddot*dt;
    q = q+qdot*dt;
    log_qddot = [log_qddot qddot];
    log_q = [log_q q];
    log_qdot = [log_qdot qdot];
    log_tau = [log_tau tau];
    i = i+1;
end

%% plot
t = 0:dt:T;
figure(1);
subplot(221)
plot(t,(log_q(1,:)-qd(1))*180/pi);
xlabel('time/s')
ylabel('deg')
title('position error 1');
% set(gca, 'Position', [0.08 0.05 0.88 0.36])
subplot(222)
plot(t,(log_q(2,:)-qd(2))*180/pi);
axis tight
ylabel('deg')
xlabel('time/s')
title('position error 2')
subplot(223)
plot(t,log_tau(1,:));
axis tight
xlabel('time/s')
title('control torque 1')
subplot(224)
plot(t,log_tau(2,:));
axis tight
xlabel('time/s')
title('control torque 2')




