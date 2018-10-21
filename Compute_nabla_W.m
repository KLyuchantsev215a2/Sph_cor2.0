function [nabla_W]=Compute_nabla_W(i,j,x,h,beta)
%kernel of SPH as a function
% of the coordinate particle i and j

% input:  %a number initial particle
          %b namber calculated particle 
          %h  blurring radius
          %x coordinate all particle 
          %betta normal for the derivative 
% output: W = the force of the particle j effect on the initial i

nabla_W=0;
    
r=zeros(1,2);
rn=zeros(1,2);
r(1,1)=x(1,i)-x(1,j);
r(1,2)=x(2,i)-x(2,j);

q=norm(r,2)/h;
C=1/(pi*h*h);
rn=r/norm(r,2);
if (q>0) && (q<1)
  nabla_W=C*(10 / 7)*(-3*q+9/4*q*q)*(rn(1,beta))/h;%Review of Development of the Smooth Particle Hydrodynamics(SPH) MethodRade Vignjevic
elseif (q >= 1) && (q < 2)
   nabla_W = C*(10 / 7)*(-3/4)*(2 - q)^2*(rn(1,beta))/h;
end
