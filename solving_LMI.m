%%
clc;
clearvars 
LMIs=3; % LMI个数
A=[-2 0;0 -0.9];A1=0.899*[-1 0;-1 -1]; % 系统矩阵
dimn=length(A);
%%
setlmis([]) 
p = lmivar(1,[dimn 1]); % 定义变量P
q = lmivar(1,[dimn 1]); % 定义变量Q
%%
lmiterm([1 1 1 p],1,A,'s'); % LMI #1
lmiterm([1 1 1 q],1,1);
lmiterm([1 1 2 p],1,A1); 
lmiterm([1 2 2 q],-1,1);
%%
lmiterm([-2 1 1 p],1,1); % LMI #2
%%
lmiterm([-3 1 1 q],1,1); % LMI #3
%% LMI求解
lmis = getlmis;
[tmin,xfeas] = feasp(lmis,[0,0,10,0,0],0); % 使用feasp函数求解
% [tmin,xfeas] = feasp(lmis);
P = dec2mat(lmis,xfeas,p); %使用dec2mat从xfeas导出矩阵变量的可行值
Q = dec2mat(lmis,xfeas,q);
disp('P = ');disp(P);
disp('Q = ');disp(Q);
%% 检验求得的Xopt是否满足线性矩阵不等式约束  *************
evlmi=evallmi(lmis,xfeas); % 将求解的xfeas带入LMI中
Lhs=cell(1,LMIs);Rhs=cell(1,LMIs);
Lhs_Rhs=cell(1,LMIs);
%------------------------------------------检验每一个LMI的负定性
sum=0;label=zeros(1,LMIs);
store=cell(1,LMIs); % 记录各个LMI的左右特征值之差
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