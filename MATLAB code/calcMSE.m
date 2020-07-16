


function mse = calcMSE(CRdata, t, X)
   CRtest = X(1)*exp(X(2)*t)-X(3); 
   mse    = mean(abs(CRdata-CRtest).^2);
end