%% Joint state and switching signal observer
function [xhat_kplus1, sigmaHatk] = JointObsv(xhat_k,Y,matrices)
   
   A1 = matrices{1}{1}; A2 = matrices{1}{2};
   C1 = matrices{2}{1}; C2 = matrices{2}{2};
   L1 = matrices{3}{1}; L2 = matrices{3}{2};
   
   % mode estimation
   O11 = [C1; C1*A1]; O12 = [C1; C2*A1];
   O21 = [C2; C1*A2]; O22 = [C2; C2*A2];

   ImO11sym = colspace(sym(O11)); ImO12sym = colspace(sym(O12));  % Image space
   ImO21sym = colspace(sym(O21)); ImO22sym = colspace(sym(O22)); 

   ImO11 = double(ImO11sym); ImO12 = double(ImO12sym);
   ImO21 = double(ImO21sym); ImO22 = double(ImO22sym);

   if rank([ImO11, Y]) == rank(ImO11)
       sigmaHatk = 1;
   elseif rank([ImO12, Y]) == rank(ImO12)
       sigmaHatk = 2;
   elseif rank([ImO21, Y]) == rank(ImO21)
       sigmaHatk = 1;
   elseif rank([ImO22, Y]) == rank(ImO22)
       sigmaHatk = 2;
   end
 
   % state estimation
   if sigmaHatk == 1
       A = A1; C = C1; L = L1;
   elseif sigmaHatk == 2
       A = A2; C = C2; L = L2;
   end
   yk = Y(end);
   xhat_kplus1 = A*xhat_k + L*(yk - C*xhat_k);
    
end