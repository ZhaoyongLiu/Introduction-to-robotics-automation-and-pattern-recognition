%% YAMLIP version
% Link: https://zhuanlan.zhihu.com/p/536058938
% Author: Liu Zhaoyong
% Date: 2023.1.20  Version: 1.0

%%
current_dir=pwd;
cd('D:\MATLAB\R2021a\bin\SDPT3-4.0\SDPT3-4.0');
run('Installmex.m')
run('startup.m')
cd(current_dir);

%% Proposition 1
clc;clear;

beta=0.899;
A=[-2 0;0 -0.9];
A1=beta*[-1 0;-1 -1];
n=size(A,1); 

%% Decision variables 
P=sdpvar(n); 
Q=sdpvar(n);  

W=blkvar; 
W(1,1)=A'*P+P*A+Q;
W(1,2)=P*A1;
W(2,2)=-Q;
W=sdpvar(W); 

%% Solution of LMIs
LMIs=[P>=0,Q>=0,W<=0]; 
options=sdpsettings('solver','sdpt3','verbose',0);
sol=optimize(LMIs,[],options); 

if sol.problem == 0
    [primal,~]=check(LMIs); 
    if min(primal)>=0
        disp('***All LMIs are negative definite!!!***'); 
    else
        disp('*** Infeasible! ***');
    end
else
    yalmiperror(sol.problem) 
    disp('***Infeasible!!!***'); 
end
