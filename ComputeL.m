function [Ln1]=ComputeL(L,i,j,beta,m,v,rho,nabla_W)
L(1,beta,i)=L(1,beta,i)+(m/rho(1,j)*(v(1,j))*nabla_W(beta));  
L(2,beta,i)=L(2,beta,i)+(m/rho(1,j)*(v(2,j))*nabla_W(beta)); 

Ln1=L;