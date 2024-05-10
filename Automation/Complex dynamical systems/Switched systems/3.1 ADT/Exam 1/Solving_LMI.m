%% Solving LMIs of example 1
% Link: https://zhuanlan.zhihu.com/p/539299349
% Author          Date          Version     Modification
% Zhaoyong Liu    Nov-3-2023    1.0    

%%
clc; clearvars 
LMIs=6; % the number of LMIs
A1=[-2 1;-2 -2]; A2=[-2 2;-1 -2]; % system matrices
dimn=length(A1);
lambda=3.999; mu=1.9986;  % parameters

%%
setlmis([]) 
p1 = lmivar(1,[dimn 1]); % define variable P1
p2 = lmivar(1,[dimn 1]); 
%%
lmiterm([1 1 1 p1],1,A1,'s'); % LMI #1
lmiterm([1 1 1 p1],lambda,1);
%%
lmiterm([2 1 1 p2],1,A2,'s'); % LMI #2
lmiterm([2 1 1 p2],lambda,1);
%% 
lmiterm([3 1 1 p1],1,1);  % LMI #3
lmiterm([-3 1 1 p2],mu,1);
%% 
lmiterm([4 1 1 p2],1,1);  % LMI #4
lmiterm([-4 1 1 p1],mu,1);
%%
lmiterm([-5 1 1 p1],1,1); % LMI #5
%%
lmiterm([-6 1 1 p2],1,1); % LMI #6

%% Solving LMIs
lmis = getlmis;
[tmin,xfeas] = feasp(lmis,[0,0,10,0,0],0); % solving with function "feasp"
% [tmin,xfeas] = feasp(lmis);
P1 = dec2mat(lmis,xfeas,p1) % exporting feasible value of matrix variable from "xfeas" using "dec2mat"
P2 = dec2mat(lmis,xfeas,p2)

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
        if notsat_num(j)==1
             fprintf('The %dst LMI is infeasible!\n',notsat_num(j));
        elseif notsat_num(j)==2
             fprintf('The %dnd LMI is infeasible!\n',notsat_num(j));
        elseif notsat_num(j)==3
             fprintf('The %drd LMI is infeasible!\n',notsat_num(j));
        else
             fprintf('The %dth LMI is infeasible!\n',notsat_num(j));
        end
    end
end