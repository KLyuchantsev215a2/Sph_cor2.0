%input:
    %rho = initial density
    %H = height of a plate
    %L = plate width
    %N = number of particles approximating the body for now 9+9=18
    %v_0 = initial velocity
%outpt:
    %dynamics of the body
    %mu,k the constants of the material for the generalized Hooke Law
    %E- elastic modulus for materials
    clear;
    
rho_0 =1.2;
v_0 = 0.001;
Time = 10;
mu = 10;             % модуль сдвига
k = 100;              % к-т объёмного сжатия
E=9*k*mu/(3*k+mu);   % модуль Юнга
sqn=12;
l=0.01;
N=sqn*sqn;
S=l*l;
%S=H*L;
m=rho_0*S/N;
h=1.2*(m/rho_0)^(1/2);
dt=0.0001;
fr=10;
  
    %coordinate(particle) initialization
    x=initialization_x(N,sqn,l);    
    v=initialization_v(N,sqn,v_0);
    x_old=x(1,1:N);
    
    rho=rho_0*ones(1,N);
    
    W_cor=zeros(N,N);
    W_old=zeros(N,N);
    nabla_W_old=zeros(2,N,N);
    nabla_W_cor=zeros(2,N,N);%2 dimension(x,y)
    
    F=zeros(2,2,N);
    T=zeros(2,2,N);
    SIG=zeros(2,2,N);
    
    nabla_W=[0,0];
    
         
    for i = 1:N
      F(1:2,1:2,i)=eye(2);
      SIG(1:2,1:2,i)=ComputeStress(F(1:2,1:2,i),mu,k);
    end
     
    
    
    for n = 1:fix(Time/dt)
                
        L=zeros(2,2,N);  
        W_cor=zeros(N,N);
        nabla_W_cor=zeros(2,N,N);
        
        
      
         
        for i = 1:N          
            sumW=0;
            
            for j = 1:N
                 sumW=sumW+m/rho(1,j)*ComputeW(i,j,x,h);
            end
            
            for j = 1:N
                 W_cor(i,j)=ComputeW(i,j,x,h)/sumW;
            end
        end
        
        for i = 1:N
            sumW=0;
            gamma=[0,0];
            Li=zeros(2,2);   % 
            diad=zeros(2,2);
            
            for j = 1:N
                   sumW=sumW+m/rho(1,j)*ComputeW(i,j,x,h);%bonet1999 57
            end
            
            for beta = 1:2
                 for j = 1:N
                     gamma(beta)=(gamma(beta)+m/rho(1,j)*Compute_nabla_W(i,j,x,h,beta));%bonet1999 57
                end
            end
            
            gamma(1:2)=gamma(1:2)/sumW;%bonet1999 57
     
           

       
            for beta = 1:2
                 for j = 1:N       
                    nabla_W_cor(beta,i,j)=(Compute_nabla_W(i,j,x,h,beta)-gamma(beta));%bonet1999 57
                end
            end
            
             for beta = 1:2
                 for j = 1:N       
                    nabla_W_cor(beta,i,j)=nabla_W_cor(beta,i,j)/sumW;%bonet1999 57
                end
            end
            
            
            
            for j = 1:N
                 diad(1,1)= nabla_W_cor(1,i,j)*x(1,j);
                 diad(1,2)= nabla_W_cor(1,i,j)*x(2,j);
                 diad(2,1)= nabla_W_cor(2,i,j)*x(1,j);
                 diad(2,2)= nabla_W_cor(2,i,j)*x(2,j);
                 Li=Li+m/rho(1,j)*diad;
            end
            
            Li=Li^(-1);
            
            
           for j = 1:N
            nabla_W_cor(1:2,i,j)=Li*nabla_W_cor(1:2,i,j);
           end

           % nabla_W_cor(1:2,i,j)=nabla_W_cor+W_cor_tmp;
            %nabla_W_cor=Compute_nabla_W_cor(i,N,x,m,h,rho,nabla_W_cor);
       end
        
     %   for i = 1:N
          %   for j = 1:N
          %     for beta = 1:2                 
           %          v=ComputeVelocity(i,j,beta,dt,m,v,rho,SIG,nabla_W_cor);
            %    end
           % end
         %end
         
          if n==4
              number=1;
             for i=1:N
                 W_cor=ComputeW_cor(i,N,x,m,h,rho,W_cor);
                  for j=1:N
                  W_old(i,j)=ComputeW(i,j,x,h);   
                  
                  for beta=1:2
                   nabla_W_old(beta,i,j)=Compute_nabla_W(i,j,x,h,beta);
                  end
                  
                  end
             end
                  x_coord(1:N) = x(1,1:N);
                  y_coord(1:N) = x(2,1:N);
                  subplot(2,2,1);
                  scatter(x_coord,y_coord);
                  tri=delaunay(x_coord,y_coord);
                  trisurf(tri,x_coord,y_coord,W_cor(1:N,number));
                  subplot(2,2,2);
                  trisurf(tri,x_coord,y_coord,W_old(1:N,number));
                  subplot(2,2,3);
                  trisurf(tri,x_coord,y_coord,nabla_W_cor(1,1:N,number));
                  subplot(2,2,4);
                  trisurf(tri,x_coord,y_coord,nabla_W_old(1,1:N,number));
          end 
          
         for i = 1:N
             for j = 1:N
                for beta = 1:2
                      nabla_W(1:2)=nabla_W_cor(1:2,i,j);%Compute_nabla_W(i,j,x,h,beta);% nabla_W_cor(beta,i,j);
                      L=ComputeL(L,i,j,beta,m,v,rho,nabla_W);  
                end
             end
         end
         
      %   for i = 1:N
          %   for j = 1:N
           %     for beta = 1:2
            %           nabla_W=Compute_nabla_W(i,j,x,h,beta);
           %            Fdot=ComputeF(F,i,j,beta,m,v,rho,nabla_W);  
           %    end
          %   end
         % end
      %   F=F+dt*Fdot;
          
         for i = 1:N
             %F(1:2,1:2,i)= F(1:2,1:2,i)+dt* L(1:2,1:2,i)*F(1:2,1:2,i);    
            dtLL = dt* L(1:2,1:2,i);
            F(1:2,1:2,i)= expm(dtLL)*F(1:2,1:2,i);
          end
        
        for i = 1:N
             SIG(1:2,1:2,i)=ComputeStress(F(1:2,1:2,i),mu,k);
        end
        
        for i = 1:N
             x(1,i)=x_old(1,i)-n*dt*fix((i-1)/sqn)*v_0;
             %x(2,i)=x(2,i)+dt*v(2,i);
        end
        
      
        
        for i = 1:N
             rho(1,i)=ComputeRho(m,N,W_cor,i);
        end
            
        
        x_coord(1:N) = x(1,1:N);
        y_coord(1:N) = x(2,1:N);
        subplot(2,2,1);
        scatter(x_coord,y_coord);
       
        
       detF=ones(1,N);
       for i = 1:N
              detF(1,i)=det(F(1:2,1:2,i));
       end
       
        tri=delaunay(x_coord,y_coord);
        subplot(2,2,2);
        trisurf(tri,x_coord,y_coord,detF(1,1:N));
        
        subplot(2,2,3);
        trisurf(tri,x_coord,y_coord,rho(1,1:N));
        
       errSIG=ones(1,N);
       for i = 1:N
               SIG15=SIG(1:2,1:2,2);
               errSIG(i)=norm((SIG(1:2,1:2,2)-SIG(1:2,1:2,i)));
       end
       
        subplot(2,2,4);
        trisurf(tri,x_coord,y_coord,errSIG);    
        pause(0.0000001);
        
        
    end
   % load mri 
  %  D = squeeze(D); 
  %  vtkwrite('mri.vtk', 'structured_points', 'mri', D)
  %   x = 1:100; 
   % y = sin(x); 
   % z = sqrt(x); 
   % vtkwrite('execute','polydata','lines',x,y,z); 
    main=F;    
