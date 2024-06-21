%% Calculate filter parameters
% Reference: Guang-Ren Duan, Hai-Hua Yu. LMIs in Control Systems: Analysis, Design and Applications[M]. Boca Raton: CRC Press, 2013, Page 306.
% Author          Date           Version     Modification
% Zhaoyong Liu    Jun-18-2024     1.0   

%%
clc; clear; 
LMIs = 3; gamma = 1.0e-7;
A = [-4.4697 9.3324  8.7645;
     -9.6408 2.0236  2.2076;
     -6.6937 1.3595 -3.2241];
B = [0.2309;
    -0.4584;
    -0.8390];
C = [ 0.7813 -4.9301 -4.3341;
      4.9421 -3.4866 -4.1502];
D = [0.0074;
     -0.0088];
L =[-1.1604 1.7916 -1.4033;
    -2.0351 1.5369 -0.3777];

[dimn, dimp] = size(B);
[diml, ~] = size(C);
[dimm, ~] = size(L);

%%
setlmis([]) 
r = lmivar(1,[dimn 1]); % define variable R
x = lmivar(1,[dimn 1]);
m = lmivar(2,[dimn dimn]);
n = lmivar(2,[dimm dimn]);
z = lmivar(2,[dimn diml]);
df = lmivar(2,[dimm diml]);

%%
lmiterm([1 1 1 r],1,A,'s');  % LMI #1 
lmiterm([1 1 1 z],1,C,'s');
lmiterm([1 2 1 -m],1,1);
lmiterm([1 2 1 z],1,C);
lmiterm([1 2 1 x],1,A);
lmiterm([1 2 2 m],1,1,'s');
lmiterm([1 3 1 r],B',1);
lmiterm([1 3 1 -z],D',1);
lmiterm([1 3 2 x],B',1);
lmiterm([1 3 2 -z],D',1);
lmiterm([1 3 3 0],-gamma*eye(dimp));
lmiterm([1 4 1 0],L);
lmiterm([1 4 1 df],-1,C);
lmiterm([1 4 2 n],-1,1);
lmiterm([1 4 3 df],-1,D);
lmiterm([1 4 4 0],-gamma*eye(dimp));

%%
lmiterm([-2 1 1 x],1,1);  % LMI #2

%%
lmiterm([3 1 1 x],1,1);  % LMI #3
lmiterm([-3 1 1 r],1,1);

%% Solution of LMIs
lmis = getlmis;
[tmin, xopt]=feasp(lmis);

% export feasible value of matrix variable from "xfeas" using "dec2mat"
R = dec2mat(lmis,xopt,r);
X = dec2mat(lmis,xopt,x);
M = dec2mat(lmis,xopt,m);
N = dec2mat(lmis,xopt,n);
Z = dec2mat(lmis,xopt,z);
Df = dec2mat(lmis,xopt,df);

Af = X^-1*M;
Bf = X^-1*Z;
Cf = N;
save Af.mat Af
save Bf.mat Bf
save Cf.mat Cf
save Df.mat Df

%% Test whether the obtained solution satisfies the constraints of LMIs
evlmi = evallmi(lmis, xopt);  % bring the solved "xfeas" into LMIs
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