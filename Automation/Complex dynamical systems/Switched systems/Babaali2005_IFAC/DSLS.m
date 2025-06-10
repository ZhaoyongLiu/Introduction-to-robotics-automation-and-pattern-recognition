%% Discrete-time switched linear system
function xkplus1 = DSLS(xk, swiSigk, matrices) 
   A1 = matrices{1}{1};
   A2 = matrices{1}{2};

   if swiSigk == 1
       A = A1;
   elseif swiSigk == 2
       A = A2;
   end

   xkplus1 = A*xk;

end