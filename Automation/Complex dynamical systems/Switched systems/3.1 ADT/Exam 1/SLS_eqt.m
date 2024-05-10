%% System Dynamic Equation (Switched Linear System)

function dx = SLS_eqt(swi_t,x)
% Input arguments: switching signal，system state
% Output arguments: derivative of the state

A1=[-2 1;-2 -2]; % system matrices
A2=[-2 2;-1 -2];

switch swi_t   
    case 1
        A=A1;
    case 2
        A=A2;
end
dx=A*x; % system dynamic equation

