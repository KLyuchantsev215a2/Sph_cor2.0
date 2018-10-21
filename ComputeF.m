function [Fn1]=ComputeF(F,i,j,beta,m,v,rho,nabla_W)

F(1,beta,i)=(m/rho(1,j)*(v(1,j)-v(1,i))*nabla_W(beta));  
F(2,beta,i)=(m/rho(1,j)*(v(2,j)-v(2,i))*nabla_W(beta)); 

Fn1=F;