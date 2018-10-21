 function [rho] = ComputeRho(m,N,W_cor,i)

rho=0;
for j = 1:N
    rho=rho+m*W_cor(i,j); 
end

