%% LMI for Theorem 3
% Link: https://zhuanlan.zhihu.com/p/607642466
% Author: Liu Zhaoyong
% Date: 2023.2.20  Version: 1.0

%%
clc;clearvars 
LMIs=5; 
h=1.3403;   % the upper bound of tau(t)
d=0.8;      % the upper bound of the derivative of tau(t)
A=[-2 0;0 -0.9]; A1=[-1 0;-1 -1]; 
dimn=length(A);

%%
setlmis([]) 
P = lmivar(1,[dimn 1]); 
S = lmivar(1,[dimn 1]);
R = lmivar(1,[dimn 1]);
Q = lmivar(1,[dimn 1]);

%% LMI #1
lmiterm([1 1 1 P],1,A,'s'); 
lmiterm([1 1 1 S],1,1);
lmiterm([1 1 1 Q],1,1);
lmiterm([1 1 1 R],-1,1);

lmiterm([1 2 2 S],-1,1); 
lmiterm([1 2 2 R],-1,1); 

lmiterm([1 1 3 P],1,A1);
lmiterm([1 2 3 R],1,1);
lmiterm([1 3 3 Q],-(1-d),1);
lmiterm([1 3 3 R],-2,1);

lmiterm([1 1 4 R],h*A',1);
lmiterm([1 3 4 R],h*A1',1);
lmiterm([1 4 4 R],-1,1);

%%
lmiterm([-2 1 1 P],1,1); 
%%
lmiterm([-3 1 1 S],1,1); 
%%
lmiterm([-4 1 1 R],1,1); 
%%
lmiterm([-5 1 1 Q],1,1); 
%% LMI求解
lmis = getlmis;
[tmin,xfeas] = feasp(lmis,[0,0,10,0,0],0); 
% [tmin,xfeas] = feasp(lmis);
P = dec2mat(lmis,xfeas,P); 
S = dec2mat(lmis,xfeas,S);
R = dec2mat(lmis,xfeas,R);
Q = dec2mat(lmis,xfeas,Q);
disp('P = ');disp(P);

%% Test whether the obtained solution satisfies the constraint of LMIs ******
evlmi=evallmi(lmis,xfeas); 
Lhs=cell(1,LMIs);Rhs=cell(1,LMIs);
Lhs_Rhs=cell(1,LMIs);
%----------------------------------verify the negative determinacy of each LMI
sum=0;label=zeros(1,LMIs);
store=cell(1,LMIs);
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