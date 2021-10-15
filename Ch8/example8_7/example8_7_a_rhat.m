function [sys,x0,str,ts] = spacemodel(t,x,u,flag)
    switch flag
    case 0
        [sys,x0,str,ts]=mdlInitializeSizes;
    case 1
        sys=mdlDerivatives(t,x,u);
    case 3
        sys=mdlOutputs(t,x,u);
    case {2,4,9}
        sys=[];
    otherwise
        error(['Unhandled flag = ',num2str(flag)]);
    end
end

function [sys,x0,str,ts]=mdlInitializeSizes
    sizes = simsizes;
    sizes.NumContStates  = 1;
    sizes.NumDiscStates  = 0;
    sizes.NumOutputs     = 1; 
    sizes.NumInputs      = 3;% r y e
    sizes.DirFeedthrough = 0;
    sizes.NumSampleTimes = 0;
    sys = simsizes(sizes);
    x0  = 0;
    str = [];
    ts  = [];
end

function sys=mdlDerivatives(t,x,u)
    gamma = 1;
    sys = -gamma*u(3)*u(1);
end

function sys=mdlOutputs(t,x,u)
    if isnan(u)
       u = [0;0;0]; 
    end
    sys = x*u(1)
end