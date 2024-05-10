%% 计算控制器增益 K
% Reference: Guang-Ren Duan, Hai-Hua Yu. LMIs in Control Systems: Analysis, Design and Applications[M]. Boca Raton: CRC Press, 2013, Page 201.
% Author          Date           Version     Modification
% Zhaoyong Liu    May-9-2024     1.0   

%%
clc; clear; clearvars;
LMIs = 1;  % the number of LMIs
n = 20;
e1 = 20*ones(n,1);
e2 = (20: -1: 1)';
Asp = spdiags([e1 e2],-1:0,n,n);
A = full(Asp);
B = [1 zeros(1,n-1)]';
C = [zeros(1,n-1) 1];
[dimn, dimm] = size(B);
[dimp, ~] = size(C);

save('A.mat','A')
save('B.mat','B')

%%
setlmis([]) 
p = lmivar(1,[dimn 1]); % define variable P
w = lmivar(2,[dimm dimn]);

%%
lmiterm([1 1 1 p],A,1,'s');  % LMI #1 
lmiterm([1 1 1 w],B,1,'s');

%%
lmiterm([-2 1 1 p],1,1);  % LMI #2

%% Solution of LMIs
lmis = getlmis;
% [tmin,xfeas] = feasp(lmis,[0,0,10,0,0],0); % solving with function "feasp"
[tmin,xfeas] = feasp(lmis);
P = dec2mat(lmis,xfeas,p); % export feasible value of matrix variable from "xfeas" using "dec2mat"
W = dec2mat(lmis,xfeas,w);

K = W*P^-1;
save('K.mat','K')
disp('K = ');disp(K);

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