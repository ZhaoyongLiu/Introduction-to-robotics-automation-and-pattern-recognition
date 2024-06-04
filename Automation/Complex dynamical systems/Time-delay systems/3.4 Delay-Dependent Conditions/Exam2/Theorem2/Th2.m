%% LMI for Theorem 2
% Link: https://zhuanlan.zhihu.com/p/607642466
% Author: Liu Zhaoyong
% Date: 2023.2.20  Version: 1.0

%%
clc;clearvars 
LMIs=6; % the number of LMIs
h=0.9998;
mu=1;
A=[-2 0;0 -0.9]; Ad=[-1 0;-1 -1]; % system matrices
dimn=length(A);

%%
setlmis([]) 
P = lmivar(1,[dimn 1]); % define variable P
Q = lmivar(1,[dimn 1]);
Z = lmivar(1,[dimn 1]);

X11 = lmivar(2,[dimn dimn]);
X12 = lmivar(2,[dimn dimn]);
X22 = lmivar(2,[dimn dimn]);
X13 = lmivar(2,[dimn dimn]);
X23 = lmivar(2,[dimn dimn]);
X33 = lmivar(2,[dimn dimn]);

N1 = lmivar(2,[dimn dimn]);
N2 = lmivar(2,[dimn dimn]);
N3 = lmivar(2,[dimn dimn]);

T1 = lmivar(2,[dimn dimn]);
T2 = lmivar(2,[dimn dimn]);

%% LMI #1
lmiterm([1 1 1 Q],1,1); 
lmiterm([1 1 1 N1],1,1,'s');
lmiterm([1 1 1 T1],-1,A,'s');
lmiterm([1 1 1 X11],h,1);

lmiterm([1 1 2 P],1,1); 
lmiterm([1 1 2 -N2],1,1);
lmiterm([1 1 2 T1],1,1);
lmiterm([1 1 2 -T2],-A',1);
lmiterm([1 1 2 X12],h,1);

lmiterm([1 2 2 Z],h,1);
lmiterm([1 2 2 T2],1,1,'s');
lmiterm([1 2 2 X22],h,1);

lmiterm([1 1 3 -N3],1,1);
lmiterm([1 1 3 N1],-1,1);
lmiterm([1 1 3 T1],-1,Ad);
lmiterm([1 1 3 X13],h,1);

lmiterm([1 2 3 N2],-1,1);
lmiterm([1 2 3 T2],-1,Ad); 
lmiterm([1 2 3 X23],h,1);

lmiterm([1 3 3 Q],-(1-mu),1);
lmiterm([1 3 3 N3],-1,1,'s');
lmiterm([1 3 3 X33],h,1);

%%
lmiterm([-2 1 1 X11],1,1); 
lmiterm([-2 1 2 X12],1,1); 
lmiterm([-2 2 2 X22],1,1); 
lmiterm([-2 1 3 X13],1,1); 
lmiterm([-2 2 3 X23],1,1); 
lmiterm([-2 3 3 X33],1,1); 
lmiterm([-2 1 4 N1],1,1); 
lmiterm([-2 2 4 N2],1,1); 
lmiterm([-2 3 4 N3],1,1); 
lmiterm([-2 4 4 Z],1,1); 

%%
lmiterm([-3 1 1 P],1,1); 
%%
lmiterm([-4 1 1 Q],1,1); 
%%
lmiterm([-5 1 1 Z],1,1); 
%%
lmiterm([-6 1 1 X11],1,1); 
lmiterm([-6 1 2 X12],1,1); 
lmiterm([-6 2 2 X22],1,1); 
lmiterm([-6 1 3 X13],1,1); 
lmiterm([-6 2 3 X23],1,1); 
lmiterm([-6 3 3 X33],1,1); 

%% Solving LMIs
lmis = getlmis;
[tmin,xfeas] = feasp(lmis,[0,0,10,0,0],0); % solving with function "feasp"
% [tmin,xfeas] = feasp(lmis);
P = dec2mat(lmis,xfeas,P); % exporting feasible value of matrix variable from "xfeas" using "dec2mat"
Q = dec2mat(lmis,xfeas,Q);
Z = dec2mat(lmis,xfeas,Z);
X11 = dec2mat(lmis,xfeas,X11);
X12 = dec2mat(lmis,xfeas,X12);
X22 = dec2mat(lmis,xfeas,X22);
X13 = dec2mat(lmis,xfeas,X13);
X23 = dec2mat(lmis,xfeas,X23);
X33 = dec2mat(lmis,xfeas,X33);
N1 = dec2mat(lmis,xfeas,N1);
N2 = dec2mat(lmis,xfeas,N2);
N3 = dec2mat(lmis,xfeas,N3);
disp('P = ');disp(P);

%% Test whether the obtained solution satisfies the constraints of LMIs
evlmi = evallmi(lmis, xfeas);  % bring the solved "xfeas" into LMIs
Lhs = cell(1, LMIs);Rhs=cell(1, LMIs);
Lhs_Rhs = cell(1, LMIs);
%-------------------------------------------------------------------------%
%         verify the negative determinacy of each LMI
%-------------------------------------------------------------------------%
sum = 0;label = zeros(1, LMIs);
store = cell(1, LMIs); % record the difference between the left and right eigenvalues of each LMI
for i=1:LMIs
    [lhs,rhs] = showlmi(evlmi,i);
    Lhs{1,i} = lhs; Rhs{1,i} = rhs;
    Lhs_Rhs{1,i} = [lhs,rhs];
    store{1,i} = eig(lhs-rhs);
    if(eig(lhs-rhs)<0)
        sum = sum+1;
        label(i) = 1; 
    end
end
if(sum == LMIs)
    disp('****** All LMIs are negative definite! ******');
else
    disp('****** Error! ******');
    notsat_num = find(label == 0);
    for j = 1:length(notsat_num)
        if mod(notsat_num(j),10) == 1 && notsat_num(j) ~= 11
             fprintf('The %dst LMI is infeasible!\n',notsat_num(j));
        elseif mod(notsat_num(j),10) == 2 && notsat_num(j) ~= 12
             fprintf('The %dnd LMI is infeasible!\n',notsat_num(j));
        elseif mod(notsat_num(j),10) == 3 && notsat_num(j) ~= 13
             fprintf('The %drd LMI is infeasible!\n',notsat_num(j));
        else
             fprintf('The %dth LMI is infeasible!\n',notsat_num(j));
        end
    end
end