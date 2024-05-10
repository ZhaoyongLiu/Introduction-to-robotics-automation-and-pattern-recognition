%% This MATLAB program checks the feasibility of LMIs from Theorem 6
% Link: https://zhuanlan.zhihu.com/p/607642466
function [P,S,Q,R,X]=LMI_Th6(A,Ad,hm,hM,dm,dM)
% The program uses YALMIP parser (http://users.isy.liu.se/johanl/yalmip/)
% and SDPT3 solver (https://github.com/SQLP/SDPT3)

% Input: 
% A,Ad          - the parameters of the system; 
% hm,hM         - the minimum and maximum sampling interval; 
% dm,dM         - the lower and upper bound of the derivative of h(t); 

% Output: 
% P,S,Q,R,X - the matrix variables. 

n=size(A,1); 

G2=[eye(n),-eye(n),zeros(n),zeros(n),zeros(n)]; % [I -I 0 0 0]
G3=[eye(n),eye(n),zeros(n),-2*eye(n),zeros(n)]; % [I I 0 -2I 0]
G4=[zeros(n),eye(n),-eye(n),zeros(n),zeros(n)]; % [0 I -I 0 0]
G5=[zeros(n),eye(n),eye(n),zeros(n),-2*eye(n)]; % [0 I I 0 -2I]
Gamma=[G2; G3; G4; G5];
%% Decision variables 
P=sdpvar(3*n); 
S=sdpvar(n); 
Q=sdpvar(n); 
R=sdpvar(n); 
X=sdpvar(2*n,2*n,'f'); 

G0hpm=[A Ad zeros(n) zeros(n) zeros(n);...
      eye(n) -(1-dm)*eye(n) zeros(n) zeros(n) zeros(n);...
      zeros(n) (1-dm)*eye(n) -eye(n) zeros(n) zeros(n)];  % G0hp
G0hpM=[A Ad zeros(n) zeros(n) zeros(n);...
      eye(n) -(1-dM)*eye(n) zeros(n) zeros(n) zeros(n);...
      zeros(n) (1-dM)*eye(n) -eye(n) zeros(n) zeros(n)];
G1hm= [eye(n) zeros(n) zeros(n) zeros(n) zeros(n);...
      zeros(n) zeros(n) zeros(n) hm*eye(n) zeros(n);...
      zeros(n) zeros(n) zeros(n) zeros(n) (hM-hm)*eye(n)]; % G1h
G1hM= [eye(n) zeros(n) zeros(n) zeros(n) zeros(n);...
      zeros(n) zeros(n) zeros(n) hM*eye(n) zeros(n);...
      zeros(n) zeros(n) zeros(n) zeros(n) (hM-hM)*eye(n)];
  
Shat=blkvar; 
Shat(1,1)=S; 
Shat(2,2)=zeros(n); 
Shat(3,3)=-S; 
Shat(4,4)=eye(2*n); 
Shat=sdpvar(Shat); % Shat=diag{S,0,-S,0_{2n}}

Qhatm=blkvar; 
Qhatm(1,1)=Q; 
Qhatm(2,2)=-(1-dm)*Q; 
Qhatm(3,3)=eye(3*n); 
Qhatm=sdpvar(Qhatm); % Qhat=diag{Q,-(1-hp)Q,0_{3n}}

QhatM=blkvar; 
QhatM(1,1)=Q; 
QhatM(2,2)=-(1-dM)*Q; 
QhatM(3,3)=eye(3*n); 
QhatM=sdpvar(QhatM);

Rhat=blkvar;
Rhat(1,1)=R; 
Rhat(2,2)=zeros(2*n); 
Rhat=sdpvar(Rhat); % Rhat=diag{R,0_{2n}}

Rt=blkvar;
Rt(1,1)=R;
Rt(2,2)=3*R;
Rt=sdpvar(Rt); % Rt=diag{R,3R}

Phi2=blkvar;
Phi2(1,1)=Rt;
Phi2(1,2)=X;
Phi2(2,2)=Rt;
Phi2=sdpvar(Phi2); % Phi2=[Rt X;* Rt];
%% The LMI for hm,dm
Phi0mm=G1hm'*P*G0hpm+G0hpm'*P*G1hm+Shat+Qhatm+hM*G0hpm'*Rhat*G0hpm;
Phi1mm=Phi0mm-1/hM*Gamma'*Phi2*Gamma;

%% The LMI for hm,dM
Phi0mM=G1hm'*P*G0hpM+G0hpM'*P*G1hm+Shat+QhatM+hM*G0hpM'*Rhat*G0hpM;
Phi1mM=Phi0mM-1/hM*Gamma'*Phi2*Gamma;

%% The LMI for hM,dm
Phi0Mm=G1hM'*P*G0hpm+G0hpm'*P*G1hM+Shat+QhatM+hM*G0hpm'*Rhat*G0hpm;
Phi1Mm=Phi0Mm-1/hM*Gamma'*Phi2*Gamma;

%% The LMI for hM,dM
Phi0MM=G1hM'*P*G0hpM+G0hpM'*P*G1hM+Shat+QhatM+hM*G0hpM'*Rhat*G0hpM;
Phi1MM=Phi0MM-1/hM*Gamma'*Phi2*Gamma;

%% Solution of LMIs
LMIs=[P>=0, S>=0, Q>=0, R>=0, Phi2>=0, Phi1mm<=0, Phi1mM<=0, Phi1Mm<=0,Phi1MM<=0]; 
options=sdpsettings('solver','sdpt3','verbose',0);
sol=optimize(LMIs,[],options); 

if sol.problem == 0
    [primal,~]=check(LMIs); 
    if min(primal)>=0 
        disp('***All LMIs are negative definite!***'); 
    end
else
    yalmiperror(sol.problem) 
    disp('***Infeasible!***'); 
end
