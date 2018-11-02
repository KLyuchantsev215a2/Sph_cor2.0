function [W_cor]=ComputeW_cor(N,x,m,h,rho,xi,i)
%kernel of SPH as a function
% of the coordinate particle i and j

% input:  %a number initial particle
          %b namber calculated particle 
          %h  blurring radius
          %x coordinate all particle 
% output: W = the force of the particle j effect on the initial i
sumW=0;
alpha=0;
betaij=zeros(2,2);
cormat=zeros(2,2);
            
           for j = 1:N
                   ri=xi-x(1:2,j);
                   sumW=sumW+m/rho(1,j)*ComputeW(i,j,x,h)*ri;
            end
            
            for beta=1:2
                for j = 1:N
                    ri=xi-x(1:2,j);
                    cormat=cormat+(ri*ri'*m/rho(1,j)*ComputeW(i,j,x,h));%51
                 end
            end
            
            betaij=cormat^(-1)*sumW;%rename
            
            for j = 1:N
                 ri=xi-x(1:2,j);
                 alpha=alpha+(m/rho(1,j)*(1+dot(betaij,ri))*ComputeW(i,j,x,h));%52
            end
            alpha=1/alpha;
            
            for j = 1:N
                 ri=xi-x(1:2,j);
                 W_cor(1,j)=ComputeW(i,j,x,h)*alpha*(1+dot(betaij,ri));%47
            end
       



