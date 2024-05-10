%% LMI for Theorem 6
% Link: https://zhuanlan.zhihu.com/p/607642466
% Author: Liu Zhaoyong
% Date: 2022.12.26  Version: 1.0

%%
current_dir=pwd;
cd('D:\MATLAB\R2021a\bin\SDPT3-4.0');
run('Installmex.m')
run('startup.m')
cd(current_dir);

%% Parameters from Example 2
clc;clear;close all;

A=[-2 0; 0 -0.9]; Ad=[-1 0;-1 -1]; 

% Theorem 7
hm=0; hM=2.1280; % the upper and lower bounds of the time delay
dM=1; dm=-dM;    % the upper and lower bounds of the derivatives of the time delay
[P,S,Q,R,X]=LMI_Th6(A,Ad,hm,hM,dm,dM);
P=value(P);
S=value(S);
Q=value(Q);
R=value(R);
X=value(X);
