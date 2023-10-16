(* ::Package:: *)

function COUPj14
% 1x1 coupled shcroedinger equation. Autor: Ruben Oncala
% Marletta's testexample for the SL12F code.
disp('**************************************************')
disp('Hybrid-Quarkonium meson states COUPj14:')
% disp('  Marletta''s testexample for the SL12F code')
disp('**************************************************')
format long;



  % the endpoints of the integration interval:
system.a=0.001;   
system.b=22; 
% parameters of the boundary conditions:
system.A1= eye (1);
system.A2= zeros(1,1);
system.B1= eye(1);
system.B2= zeros(1,1);
% function handle to the function returning the potential matrix
system.V=@potentialMatrix;


[m]=parameters
[j]=parameters1
disp(' ');
disp('Calculating eigenvalues with indices between 0 and 10...');
t=cputime;
[EigvData,meshData]=computeEigenvalues(system,0,10,5e-8);

disp(['Number of intervals in the mesh: ' num2str(length(meshData.h))]); % number of intervals in the mesh
% disp(['Time: ' num2str (cputime - t) ' sec']);
% disp(sprintf('k \t  E_k \t\t\t\t ErrEst \t status'))
for i=1:length(EigvData.eigenvalues)
    disp(sprintf('%-3.0f\t %-16.12f\t %-2.0e\t\t % d', EigvData.indices(i),EigvData.eigenvalues(i)/m,EigvData.errorEstimations(i),EigvData.status(i)));
end

end



function [m]=parameters
 m=1.4702;
 % m=4.8802;
end

function [j]=parameters1
j=2;
end 

function [lam0,lam1,lam3,k,pi]=parameters2
lam0=0.1;
lam1=0.1182;
lam3=0.2299;
k=0.187;
pi=3.1416;
end
%---------------------------------------------------



function r=potentialMatrix(x) % returns the potential matrix evaluated in x

[j]=parameters1;

  for i=1:4 
  r(1,1,i)=F(x(i),j-1);
  end
end

%--------------------------------------------------------------------------
function f1=Vg(x)
  [m]=parameters;
  if m==1.4702
  Eo=2.6984;
  else
  Eo=9.5325;    
  end
  f1=-0.489/x+0.187*x+Eo;  
  end

  function f2=Vs(x)
  [m]=parameters;
  if m==1.4702
  Eo=3.499539;
  else
  Eo=10.333696;    
  end 
  f2=0.0611/x+0.187*x+Eo;
  end

  function f3=Vp(x)
  [m]=parameters;
  if m==1.4702
  Eo=3.491169;
  else
  Eo=10.325269;    
  end 
  b1=0.0696430221609656;
  b2=-1.4593432845775876;
  a1=-0.06732994686962318;
  a2=0.014330609468130364;  
  f3=0.187*x+(0.0611/x)*(1.0+b1*x+b2*x.^2)/(1.0+a1*x+a2*x.^2)+Eo;
  end

  function f4=Vq(x)
  f4=Vp(x)-Vs(x);
  end
  
  function M1=F(x,j)
  [m]=parameters;
  M1=(j*(j+1))/(x.^2)+m*Vp(x);
  end

  function M2=R(x,j)
  [m]=parameters;
  M2=(j*(j+1))/(x.^2)+m*Vg(x);
  end

  function M3=A(x,j)
  [m]=parameters;
  M3=(j*(j-1))/(x.^2)+m*Vq(x)*(j+1)/(2*j+1)+m*Vs(x);
  end

  function M4=B(x,j)
  [m]=parameters;
  M4=m*Vq(x)*sqrt (j*(j+1))/(2*j+1);
  end

  function M5=CC(x,j)
  [m]=parameters;
  M5=((j+1)*(j+2))/(x.^2)+m*Vq(x)*(j)/(2*j+1)+m*Vs(x);
  end

  function M6=G(x,j)
  [m]=parameters;
  M6=2/(x.^2)+m*Vs(x);
  end

  function S1=Vps(x)
  [m]=parameters;
  [lam0,lam1,lam3,k,pi]=parameters2;
  k1=lam0*lam0/m;
  k2=(2*lam0*lam0/lam1)*sqrt(k/(pi*pi*pi));
  S1=-k1/(1+k2*x.^2);
  end
  
  function S2=Vss(x)
  [m]=parameters;
  [lam0,lam1,lam3,k,pi]=parameters2;
  kk1=lam0*lam0/m;
  kk2=(lam0*lam0/lam3)*k/(pi*pi);
  S2=-kk1/(1+kk2*x.^3);
  end

  function S3=Vqs(x)
  S3=-Vps(x)+Vss(x);
  end

  function L1=aap(j)
  [m]=parameters;
  L1=m*j/(2*j+1);
  end

  function L2=aa(j)
  [m]=parameters;
  L2=m*(j+1)/(2*j+1);
  end
  
  function L3=bb(j)
  [m]=parameters;
  L3=-m*sqrt (j*(j+1))/(2*j+1);
  end

  function G1=alpha(j)
  [m]=parameters;
  G1=m*sqrt((j*(j-1))/((2j-1)*(2*j+1)));
  end

  function G2=beta(j)
  [m]=parameters;
  G2=-m*j/sqrt((2*j-1)*(2*j+1));
  end

  function G3=gamma(j)
  [m]=parameters;
  G3=-m*(j+1)/sqrt((2*j+3)*(2*j+1));
  end

  function G4=delta(j)
  [m]=parameters;
  G4=m*sqrt(((j+2)*(j+1))/((2*j+3)*(2*j+1)));
  end

  function G5=omega(j)
  [m]=parameters;
  G5=-m*sqrt((2*j-1)/(2*j+1));
  end

  function G5=eta(j)
  [m]=parameters;
  G5=-m*sqrt((2*j+3)/(2*j+1));
  end
