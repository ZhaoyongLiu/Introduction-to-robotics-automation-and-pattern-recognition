%% Lyapunov-Krasovskii Theorem
% Link: https://zhuanlan.zhihu.com/p/536058938
% Author: Liu Zhaoyong
% Date: 2023.2.20  Version: 1.0

%%
clc;clearvars 
LMIs=3; % the number of LMIs
A=[-2 0;0 -0.9];A1=0.899*[-1 0;-1 -1]; % system matrices
dimn=length(A);

%%
setlmis([]) 
p = lmivar(1,[dimn 1]); % define variable P
q = lmivar(1,[dimn 1]); 

%%
lmiterm([1 1 1 p],1,A,'s'); % LMI #1
lmiterm([1 1 1 q],1,1);
lmiterm([1 1 2 p],1,A1); 
lmiterm([1 2 2 q],-1,1);
%%
lmiterm([-2 1 1 p],1,1); % LMI #2
%%
lmiterm([-3 1 1 q],1,1); % LMI #3
%% Solving LMIs
lmis = getlmis;
[tmin,xfeas] = feasp(lmis,[0,0,10,0,0],0); % solving with function "feasp"
% [tmin,xfeas] = feasp(lmis);
P = dec2mat(lmis,xfeas,p); % exporting feasible value of matrix variable from "xfeas" using "dec2mat"
Q = dec2mat(lmis,xfeas,q);
disp('P = ');disp(P);
disp('Q = ');disp(Q);

%% Test whether the obtained solution satisfies the constraint of LMIs ******
evlmi=evallmi(lmis,xfeas); % bring the solved "xfeas" into LMIs
Lhs=cell(1,LMIs);Rhs=cell(1,LMIs);
Lhs_Rhs=cell(1,LMIs);
%----------------------------------verify the negative determinacy of each LMI
sum=0;label=zeros(1,LMIs);
store=cell(1,LMIs); % record the difference between the left and right eigenvalues of each LMI
for i=1:LMIs
    [lhs,rhs]=showlmi(evlmi,i);
    Lhs{1,i}=lhs;Rhs{1,i}=rhs;
    Lhs_Rhs{1,i}=[lhs,rhs];
    store{1,i}=eig(lhs-rhs);
    if(eig(lhs-rhs)<0)
        sum=sum+1;
        label(i)=1;
    end
end
if(sum==LMIs)
    disp('****** All LMIs are negative definite! ******');
else
    disp('****** Error! ******');
    notsat_num=find(label==0);
    for j=1:length(notsat_num)
        fprintf('The %dth LMI is infeasible!\n',notsat_num(j));
    end
end