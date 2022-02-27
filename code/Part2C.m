close all;
clear;

w=2;
l=3;
V0=1;

width=[0.7 0.6 0.5 0.4 0.3 0.2];
curr=zeros(length(width),1);
size=zeros(length(width),1);

for m=1:length(width)
  dx=0.05;
  dy=dx;
  ny=round(w/dx);
  nx=round(l/dy);

  cMap=ones(nx,ny);
  
  middle=ny/2;
  upper=middle-ny*width(m)/2;
  lower=middle+ny*width(m)/2
  
  cMap(nx/3:2*nx/3  , 1:upper)=1e-2; 
  cMap(nx/3:2*nx/3 , lower:ny)=1e-2;

  % G-matrix, relates the value of a node to all other nodes
  G=sparse(nx*ny,nx*ny);
  % F-vector, the boundary conditons
  F = sparse(nx*ny,1);

  % Populate G matrix
  % Code from https://github.com/L28E/4700Code/blob/master/CondCode/GetCurrents.m
  for i = 1:nx
      for j = 1:ny
          n = j + (i - 1) * ny;

          if i == 1
              G(n, :) = 0;
              G(n, n) = 1;
              F(n)=V0;
          elseif i == nx
              G(n, :) = 0;
              G(n, n) = 1;
          elseif j == 1
              nxm = j + (i - 2) * ny;
              nxp = j + (i) * ny;
              nyp = j + 1 + (i - 1) * ny;

              rxm = (cMap(i, j) + cMap(i - 1, j)) / 2.0;
              rxp = (cMap(i, j) + cMap(i + 1, j)) / 2.0;
              ryp = (cMap(i, j) + cMap(i, j + 1)) / 2.0;

              G(n, n) = -(rxm+rxp+ryp);
              G(n, nxm) = rxm;
              G(n, nxp) = rxp;
              G(n, nyp) = ryp;

          elseif j ==  ny
              nxm = j + (i - 2) * ny;
              nxp = j + (i) * ny;
              nym = j - 1 + (i - 1) * ny;

              rxm = (cMap(i, j) + cMap(i - 1, j)) / 2.0;
              rxp = (cMap(i, j) + cMap(i + 1, j)) / 2.0;
              rym = (cMap(i, j) + cMap(i, j - 1)) / 2.0;

              G(n, n) = -(rxm + rxp + rym);
              G(n, nxm) = rxm;
              G(n, nxp) = rxp;
              G(n, nym) = rym;
          else
              nxm = j + (i-2)*ny;
              nxp = j + (i)*ny;
              nym = j-1 + (i-1)*ny;
              nyp = j+1 + (i-1)*ny;

              rxm = (cMap(i,j) + cMap(i-1,j))/2.0;
              rxp = (cMap(i,j) + cMap(i+1,j))/2.0;
              rym = (cMap(i,j) + cMap(i,j-1))/2.0;
              ryp = (cMap(i,j) + cMap(i,j+1))/2.0;

              G(n,n) = -(rxm+rxp+rym+ryp);
              G(n,nxm) = rxm;
              G(n,nxp) = rxp;
              G(n,nym) = rym;
              G(n,nyp) = ryp;
          end

      end
  end

  % Solve for voltage
  V = G\F;
  v_surf=zeros(nx,ny);
  for x=1:nx
      for y=1:ny
        n = y + (x - 1) * ny;
        v_surf(x,y)=V(n);
      end      
  end
  
  [Ex,Ey]=gradient(-v_surf');
  Jx=cMap'.*Ex;
  Jy=cMap'.*Ey;
  
%   % conductivity Plot 
%   figure()
%   subplot(2,2,1);
%   surf(cMap','EdgeColor','none');
%   title('σ(x,y)');
%   ylabel('W');
%   xlabel('L');
% 
%   % Voltage Plot 
%   subplot(2,2,2);
%   surf(v_surf','EdgeColor','none');
%   title('V(x,y)');
%   ylabel('W');
%   xlabel('L');
% 
%   % E field Plot   
%   subplot(2,2,3);
%   quiver(Ex,Ey);
%   title('E(x,y)' );
%   ylabel('W');
%   xlabel('L');
% 
%   % Current Density plot 
%   subplot(2,2,4);  
%   quiver(Jx,Jy);
%   title('J(x,y)' );
%   ylabel('W');
%   xlabel('L');

  % Calculate Current
  eFlowx = cMap' .* Ex;
  eFlowy = cMap' .* Ey;

  C0 = sum(eFlowx(:,1))
  Cnx = sum(eFlowx(:,nx))
  Curr(m) = (C0 + Cnx) * 0.5;
  size(m)=nx*ny;

end
figure()
plot(width,Curr)
title('Current vs. Channel Width' );
ylabel('Current (A)');
xlabel('Channel Width (proportion of width)');