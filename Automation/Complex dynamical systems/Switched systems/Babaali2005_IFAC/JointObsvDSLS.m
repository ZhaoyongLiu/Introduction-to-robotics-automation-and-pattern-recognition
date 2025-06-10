%% Joint observer for DSLS
function [t, x, xhat, swiSigHat] = JointObsvDSLS(DSLS,JointObsv,tspan,swiSig,initialValues,matrices)
   
   k0 = tspan(1); kf = tspan(2);
   t = k0:1:kf;

   x0 = initialValues{1};
   swiSigHat0 = initialValues{2};
   xhat0 = initialValues{3};

   swiSigHat = zeros(2,kf-k0+1); swiSigHat(1,:) = k0:1:kf;
   swiSigHat(2,1) = swiSigHat0; 
   
   [dimn,~] = size(matrices{1}{1});
   x = zeros(dimn,kf-k0+1); x(:,1) = x0;
   xhat = zeros(dimn,kf-k0+1); xhat(:,1) = xhat0;

   [dimp,~] = size(matrices{2}{1});
   C1 = matrices{2}{1};
   C2 = matrices{2}{2};

   if swiSig(2,1) == 1
       y0 = C1*x0;
   elseif swiSig(2,1) == 2
       y0 = C2*x0;
   end
   y = zeros(dimp,kf-k0+1); y(:,1) = y0;

   if swiSigHat(2,1) == 1
       yhat0 = C1*xhat0;
   elseif swiSigHat(2,1) == 2
       yhat0 = C2*xhat0;
   end
   yhat = zeros(dimp,kf-k0+1); yhat(:,1) = yhat0;

   N = 2;  % The size of sliding window

   for i = k0+1:1:kf

       % state evolution
       xkplus1 = DSLS(x(:,i), swiSig(2,i), matrices);
       x(:,i+1) = xkplus1;
       if swiSig(2,i+1) == 1
             y(:,i+1) = C1*x(:,i+1);
       elseif swiSig(2,i+1) == 2
             y(:,i+1) = C2*x(:,i+1);
       end
       
       if i < N
           xhat(:,i+1) = xhat0;
           continue;
       end
       Y = [y(i-N+1); y(i)];   % [y(k-N+1),...,y(k)]'

       % joint observer
       [xhat_kplus1, sigmaHatk] = JointObsv(xhat(:,i),Y,matrices);
       
       swiSigHat(2,i) = sigmaHatk;
       xhat(:,i+1) = xhat_kplus1;

       if i == kf
           Y = [y(i+1-N+1); y(i+1)]; 
           [~, sigmaHatk] = JointObsv(xhat(:,i+1),Y,matrices);
           swiSigHat(2,i+1) = sigmaHatk;
           continue
       end

   end

end