% example 9.2
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
Lambda = 20*eye(2);

log_qddot = [];log_qdot = [];log_q = [];log_tau = [];
log_qd = []; log_ahat = [];

%% estimated parameters
acc = 0.8;
m1hat = 1*acc;
mehat = 2*acc;

a1hat = I1+m1hat*lc1^2+Ie+mehat*lce^2+mehat*l1^2;
a2hat = Ie+mehat*lce^2;
a3hat = mehat*l1*lce*cos(deltae);
a4hat = mehat*l1*lce*sin(deltae);

%% initial
q = [0;0];qdot = [0;0];

%% main loop
dt = 0.001;T = 3;
i = 1;
Phi = 0.05;
k = [100*Phi;100*Phi];
Yita = diag([0.03,0.05,0.1,0.3]);
ahat = [0;0;0;0];
for t = 0:dt:T
    H11hat = a1hat+2*a3hat+cos(q(2))+2*a4hat*sin(q(2));
    H12hat = a2hat+a3hat*cos(q(2))+a4hat*sin(q(2));
    H21hat = H12hat;
    H22hat = a2hat;
    hhat = a3hat*sin(q(2))-a4hat*cos(q(2));
    Hhat = [H11hat H12hat;H21hat H22hat];
    Chat = [-hhat*qdot(2) -hhat*(qdot(1)+qdot(2));hhat*qdot(1) 0];
  
    qd = [pi/6*(1-cos(2*pi*t));pi/4*(1-cos(2*pi*t))];
    qd_dot = [pi/6*2*pi*sin(2*pi*t);pi/4*2*pi*sin(2*pi*t)];
    qd_ddot = [pi/6*(2*pi)^2*cos(2*pi*t);pi/4*(2*pi)^2*cos(2*pi*t)];
    
    s = qdot-qd_dot+Lambda*(q-qd);
    qrdot = qd_dot-Lambda*(q-qd);
    qrddot = qd_ddot-Lambda*(qdot-qd_dot);

    Y11 = qrddot(1);
    Y12 = qrddot(2);
    Y21 = 0;
    Y22 = qrddot(1)+qrddot(2);
    Y13 = (2*qrddot(1)+qrddot(2))*cos(q(2))-(qdot(2)*qrdot(1)+qdot(1)*qrdot(2)+qdot(2)*qrdot(2))*sin(q(2));
    Y14 = (2*qrddot(1)+qrddot(2))*sin(q(2))-(qdot(2)*qrdot(1)+qdot(1)*qrdot(2)+qdot(2)*qrdot(2))*cos(q(2));
    Y23 = qrddot(1)*cos(q(2))+qdot(1)*qrdot(1)*sin(q(2));
    Y24 = qrddot(1)*sin(q(2))+qdot(1)*qrdot(1)*cos(q(2));
    Y = [Y11 Y12 Y13 Y14; Y21 Y22 Y23 Y24];
    
    tau = Y*ahat-Kd*s;
    ahat = ahat-Yita*Y'*s*dt;
    
%     tau = -Kp*(q-qd)-Kd*qdot;
    qddot = dynamics(tau,q,qdot);
    qdot = qdot+qddot*dt;
    q = q+qdot*dt;
    log_qd = [log_qd qd];
    log_qddot = [log_qddot qddot];
    log_q = [log_q q];
    log_qdot = [log_qdot qdot];
    log_tau = [log_tau tau];
    log_ahat = [log_ahat ahat];
    i = i+1;
end

%% plot
t = 0:dt:T;
figure(1);
subplot(221)
plot(t,log_qd(1,:)*180/pi);
xlabel('time/s')
ylabel('deg')
title('desired q 1');
% set(gca, 'Position', [0.08 0.05 0.88 0.36])
subplot(222)
plot(t,log_qd(2,:)*180/pi);
axis tight
ylabel('deg')
xlabel('time/s')
title('desired q 2')

figure(2);
subplot(221)
plot(t,(log_q(1,:)-log_qd(1,:))*180/pi);
xlabel('time/s')
ylabel('deg')
title('position error 1');
% set(gca, 'Position', [0.08 0.05 0.88 0.36])
subplot(222)
plot(t,(log_q(2,:)-log_qd(2,:))*180/pi);
axis tight
ylabel('deg')
xlabel('time/s')
title('position error 2')
subplot(223)
plot(t,log_tau(1,:));
xlabel('time/s')
title('control torque 1')
subplot(224)
plot(t,log_tau(2,:));
xlabel('time/s')
title('control torque 1')

figure(3);
subplot(221)
plot(t,log_ahat(1,:));
hold on
plot(t,a1*ones(1,length(t)));
hold off
xlabel('time/s')
handle = title('$\hat{a}_1$');
set(handle,'Interpreter','latex','FontSize',12);
% set(gca, 'Position', [0.08 0.05 0.88 0.36])
subplot(222)
plot(t,log_ahat(2,:));
hold on
plot(t,a2*ones(1,length(t)));
hold off
xlabel('time/s')
handle = title('$\hat{a}_2$');
set(handle,'Interpreter','latex','FontSize',12);
subplot(223)
plot(t,log_ahat(3,:));
hold on
plot(t,a3*ones(1,length(t)));
hold off
xlabel('time/s')
handle = title('$\hat{a}_3$');
set(handle,'Interpreter','latex','FontSize',12);
subplot(224)
plot(t,log_ahat(4,:));
hold on
plot(t,a4*ones(1,length(t)));
hold off
xlabel('time/s')
handle = title('$\hat{a}_4$');
set(handle,'Interpreter','latex','FontSize',12);
