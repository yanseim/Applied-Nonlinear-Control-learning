%% dynamics in example9.1,9.2,9.3

function qddot = dynamics(tau,q,qdot)
    global a1 a2 a3 a4
    H11 = a1+2*a3+cos(q(2))+2*a4*sin(q(2));
    H12 = a2+a3*cos(q(2))+a4*sin(q(2));
    H21 = H12;
    H22 = a2;
    h = a3*sin(q(2))-a4*cos(q(2));
    H = [H11 H12;H21 H22];
    
    C = [-h*qdot(2) -h*(qdot(1)+qdot(2));h*qdot(1) 0];
    qddot = H\(tau-C*qdot);
%     disp(tau-C*qdot)
%     disp(qddot)
end