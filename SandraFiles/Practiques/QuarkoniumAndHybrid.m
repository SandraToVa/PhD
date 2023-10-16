function Quarkonium
% 3x3 coupled shcroedinger equation. Creador subrutines i font: Ruben Oncala
% Autor: Sandra Tomàs

disp('**************************************************')
disp('Hybrid-Quarkonium energy espectrum:')
disp('**************************************************')
format long;


  % the endpoints of the integration interval:
system.a=0.001;   
system.b=22; 
% parameters of the boundary conditions:
system.A1= eye(3);
system.A2= zeros(3,3);
system.B1= eye(3);
system.B2= zeros(3,3);
% function handle to the function returning the potential matrix
system.V=@potentialMatrix;


[m]=parameters
[j]=parameters1
[s]=parameters3


disp(' ');
disp('Calculating eigenvalues with indices between 0 and 10:');
t=cputime;
[EigvData,meshData]=computeEigenvalues(system,0,10,5e-6);

% disp(['Number of intervals in the mesh: ' num2str(length(meshData.h))]); % number of intervals in the mesh
[lam0]=parameters2
% disp(['Time: ' num2str (cputime - t) ' sec']);
% disp(sprintf('k \t  E_k \t\t\t\t ErrEst \t status'))
for i=1:length(EigvData.eigenvalues)
    disp(sprintf('%-3.0f\t %-16.12f\t %-2.0e\t\t % d', EigvData.indices(i),EigvData.eigenvalues(i)/m,EigvData.errorEstimations(i),EigvData.status(i)));
end


end

% aqui escollirem els parametres j i m del sistema ---------------------------------------------

function [m]=parameters
 m=1.4702;      %Charm
 %m=4.8802;     %Bottom
end

% TOTAL ANGULAR MOMENTUM
function [j]=parameters1
j=1;
end 

% PERCENTATGE DE HYBRID (lam0=0 desactivem mixing)
function [lam0,lam1,lam3,k,pi]=parameters2
%lam0=0.6;
lam0=0.0;
lam1=0.1182;
lam3=0.2299;
k=0.187;
pi=3.1416;
end

% Escollim cas quarkonium o cas hybrid quarkonium
function [s]=parameters3
%s=0;     %Quarkonium
s=1;     %Hybrid
end
%-----------------------------------------------------------------------------------------------



function r=potentialMatrix(x) % returns the potential matrix evaluated in x

[j]=parameters1;
[s]=parameters3;

if s==0
  for i=1:4 
  % CAS QUARKONIUM S
    r(1,1,i)=R(x(i),j); 
    r(2,2,i)=R(x(i),j); 
    r(3,3,i)=R(x(i),j); 
    r(1,2,i)=0;  
    r(2,1,i)=0;
  
  end
end
  
if s==1
  % CAS HYBRIT QUARKONIUM P^(+,0,-)
  for i=1:4 
  if j==0 
      r(1,1,i)=G(x(i),j);
      r(2,2,i)=G(x(i),j);
      r(3,3,i)=G(x(i),j);
      r(1,2,i)=0;
      r(2,1,i)=0;
  else
      r(1,1,i)=A(x(i),j);
      r(2,2,i)=CC(x(i),j); 
      r(3,3,i)=F(x(i),j); 
      r(1,2,i)=B(x(i),j);
      r(2,1,i)=B(x(i),j);
  end
  
  end

end

 % NO CANVIA
  r(1,3,i)=0;
  r(3,1,i)=0;
  r(2,3,i)=0;
  r(3,2,i)=0;

end

%---------------------------------------------------------------------------------------------
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
  
  % MIXING POTENTIAL
  function S1=Vps(x) 
  [m]=parameters;
  [lam0,lam1,lam3,k,pi]=parameters2;
  k0=lam0^2/m;
  k1=((lam0^2*2*k^(0.5))/(lam1*pi^(3/2)))^(1/2);
  S1=k0*(+1-(k1*x)^2)/(1+(k1*x)^4);
  end
  
  function S2=Vss(x)
  [m]=parameters;
  [lam0,lam1,lam3,k,pi]=parameters2;
  k0=lam0^2/m;
  k2=((lam0^2*k)/(lam3*pi^2))^(1/3);
  S2=k0*(+1+(k2*x)^2)/(1+(k2*x)^5);
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
