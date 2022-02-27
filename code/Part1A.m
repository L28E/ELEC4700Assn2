close all;
clear;

w=2;
l=3;
dx=0.05;
dy=0.05;
ny=w/dx;
nx=l/dy;

V0=1;

% G-matrix, relates the value of a node to all other nodes
G=sparse(nx*ny,nx*ny);
% F-vector, the boundary conditons
F = sparse(nx*ny,1);

% Mapping Equation, converts discretized 2D spacial domain to 1D vector in the "equation domain"
map=@(x,y) y+(x-1)*ny; 

% Populate G matrix
for x=1:nx
    for y=1:ny
        n=map(x,y); 
        if (x==1||y==1||x==nx||y==ny)
            % For boundaries, G matrix is set to 1             
            G(n,n)=1; 
        else               
            % For the rest of the G matrix, populate using finite difference
            G(n,n)=-2*(1/dx^2+1/dy^2);      
            G(n,n+1)=1/(dx^2);     % Relation to the node to the right (nxp)
            G(n,n-1)=1/(dx^2);     % Relation to the node to the left (nxm) 
            G(n,n+ny)=1/(dy^2);    % Relation to the node above (nyp)
            G(n,n-ny)=1/(dy^2);    % Relation to the node below (nym)        
        end
    end
end

% Populate F vector. All bounds are 0 except x=0
for x=1:nx
    for y=1:ny
      n=map(x,y);
      if (x==1)
        F(n)=V0;
      end 
    end      
end 
  
% Solve for voltage
V = G\F;
v_surf=zeros(nx,ny);
for x=1:nx
    for y=1:ny
      n=map(x,y);
      v_surf(x,y)=V(n);
    end      
end 

surf(v_surf','EdgeColor','none');
title('Q1A' );
ylabel('W');
xlabel('L');
colorbar;
