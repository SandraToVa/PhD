function Spin1Jcal2_1
% Matriu 9x9 (capsa 1: 5x5) Acoblament spin per a s=1 j>1. Autor: Sandra Tomàs (Font: Ruben Oncala)
% INCORPORACIÓ DELS POTENCIALS QUE DEPENEN DE L'SPIN EN ELS HYBRIDS
% No incorpora mixing. Mirar QuarkoinumAndHybrid.m si es vol
%  
disp('**************************************************')
disp('Hybrid-Quarkonium meson states with spin coupling s=1(j>1):')
disp('**************************************************')
format long;


  % the endpoints of the integration interval:
system.a=0.001;   
system.b=22;

[m]=parameters
[j]=parameters1
[s]=parameters3

% parameters of the boundary conditions:
system.A1= eye(5);
system.A2= zeros(5,5);
system.B1= eye(5);
system.B2= zeros(5,5);
 

% function handle to the function returning the potential matrix
system.V=@potentialMatrix;


disp(' ');
disp('Calculating eigenvalues with indices between 0 and 10:');
t=cputime;
[EigvData,meshData]=computeEigenvalues(system,0,10,5e-6);

% disp(['Number of intervals in the mesh: ' num2str(length(meshData.h))]); % number of intervals in the mesh
%[lam0]=parameters2
% disp(['Time: ' num2str (cputime - t) ' sec']);
% disp(sprintf('k \t  E_k \t\t\t\t ErrEst \t status'))
for i=1:length(EigvData.eigenvalues)
    disp(sprintf('%-3.0f\t %-16.12f\t %-2.0e\t\t % d', EigvData.indices(i),EigvData.eigenvalues(i)/m,EigvData.errorEstimations(i),EigvData.status(i)));
end


end

% aqui escollirem els parametres L i m i S del sistema ---------------------------------------------

function [m]=parameters
 %m=1.4702;
 m=4.8802;
end

% TOTAL ANGULAR MOMENTUM (Jcal=J EN LO MIXING DESACTIVAT)
function [j]=parameters1
j=3;
end 

% ESTATS P^(+-0) HYBRIDS (s=1)
function [s]=parameters3
s=1;
end 

%function [lam0,lam1,lam3,k,pi]=parameters2
% PERCENTATGE DE HYBRID (lam0=0 desactivem mixing)
%lam0=0.6;
%lam1=0.1182;
%lam3=0.2299;
%k=0.187;
%pi=3.1416;
%end
%-----------------------------------------------------------------------------------------------



function r=potentialMatrix(x) % returns the potential matrix evaluated in x

[j]=parameters1;
[s]=parameters3;

  for i=1:4 
  if s==1 % CAS HYBRIT QUARKONIUM P^(+,0,-)
    if j>=2 
      r(1,1,i)=A(x(i),j-1)+Vhf(x(i))+2*I9(j);
      r(2,2,i)=CC(x(i),j-1)-FF(x(i),j)+2*I1(j);
      r(3,3,i)=F(x(i),j)-GG(x(i),j)+2*E5(j);
      r(4,4,i)=A(x(i),j+1)-HH(x(i),j)+2*A9(j);
      r(5,5,i)=CC(x(i),j+1)+Vhf(x(i))+2*A1(j);
      
      r(1,2,i)=B(x(i),j-1)+I3(j)+I7(j);
      r(2,1,i)=r(1,2,i);
      r(2,3,i)=FG(x(i),j)+F4(j)+G2(j);
      r(3,2,i)=r(2,3,i);
      r(3,4,i)=GH(x(i),j)+B8(j)+D6(j);
      r(4,3,i)=r(3,4,i);
      r(4,5,i)=B(x(i),j+1)+A3(j)+A7(j);
      r(5,4,i)=r(4,5,i);
      
      r(1,3,i)=F6(j);
      r(3,1,i)=r(1,3,i);
      r(1,4,i)=0;
      r(4,1,i)=r(1,4,i);
      r(1,5,i)=0;
      r(5,1,i)=r(1,5,i);
      
      r(2,4,i)=0;
      r(4,2,i)=r(2,4,i);
      r(2,5,i)=0;
      r(5,2,i)=r(2,5,i);
      
      r(3,5,i)=D4(j);
      r(5,3,i)=r(3,5,i);
      
    end
  end
  

  end
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

% EQUACIÓ 8
  function f2=Vs(x)
  [m]=parameters;
  if m==1.4702
  Eo=3.499539;
  else
  Eo=10.333696;    
  end 
  f2=0.0611/x+0.187*x+Eo;
  end
%EQUACIÓ 11
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

  
   % EFECTES HIPERFINS:
  
  % PRIMER POTENCIAL DE SPIN ----------------------------------------------
  
  function HF1=Vhf(x) 
   % Control de l'efecte del 1r potencial
  [m]=parameters;
  if m==1.4702
  HF1=0;
  else
  HF1=0;    
  end 
  
  end 
  
  function H1=BB(x,j)
  H1=Vhf(x)/(j+1);
  end
  function H2=BC(x,j)
  H2=Vhf(x)*((sqrt(j*(j+2)))/(j+1));
end
function H3=DD(x,j)
  H3=Vhf(x)/j;
  end
  function H4=DE(x,j)
  H4=Vhf(x)*(sqrt((j^2)-1)/j);
  end
  function H5=FF(x,j)
  H5=Vhf(x)*((j-1)/j);
  end
  function H6=FG(x,j)
  H6=Vhf(x)*(((j+1)*(sqrt(2*j-1)))/(j*(sqrt(2*j+1))));
  end
  function H7=GG(x,j)
  H7=Vhf(x)/(j*(j+1));
end
function H8=GH(x,j)
  H8=Vhf(x)*((j*sqrt(2*j+3))/((j+1)*sqrt(2*j+1)));
end  
function H9=HH(x,j)
  H9=Vhf(x)*((j+2)/(j+1));
end

% SEGON POTENCIAL DE SPIN -------------------------------------------------

function HF2=Vhf2(~) 
% Efecte del segon potencial
  [m]=parameters;
  if m==1.4702
  HF2=0.05;
  else
  HF2=0.009;    
  end 
  
  end 
 function AA1=A1(j)
  AA1=-Vhf2*(((3+j)/(9+6*j)));
 end
   function AA3=A3(j)
  AA3=Vhf2*((sqrt((1+j)*(2+j)))/(3+2*j));
   end
   function AA5=A5(j)
  AA5=-Vhf2/(3+3*j);
   end
   function AA7=A7(j)
  AA7=-Vhf2*(((2+j)^(3/2))/(sqrt(1+j)*(3+2*j)));
   end
   function AA9=A9(j)
  AA9=Vhf2*((j*(2+j))/(9+15*j+6*j^2));
   end
   function BB4=B4(j)
  BB4=-Vhf2*((sqrt(j)*((2+j)^(3/2)))/(3*(1+j)*(1+2*j)));
   end
   function BB6=B6(j)
  BB6=Vhf2*((j*sqrt(1+1/(1+j)))/(1+2*j));
   end
   function BB8=B8(j)
  BB8=Vhf2*(j/(3*(j+1)))*(((2*j+3)/(2*j+1))^(1/2));
   end
   function DD2=D2(j)
  DD2=Vhf2*((sqrt(j*(2+j)))/(3+3*j));
  end
   function DD4=D4(j)
  DD4=Vhf2*(j*(2+j))/(sqrt((1+j)*(2+j)*(1+2*j)*(3+2*j)));
  end
  function DD6=D6(j)
  DD6=-Vhf2*(j^2)/(3*(1+j)*sqrt(3+4*j*(2+j)));
   end
   function EE1=E1(j)
  EE1=-Vhf2*(2+j)/(3+9*j+6*(j^2));
   end
   function EE3=E3(j)
  EE3=Vhf2*(sqrt(j/(1+j)))/(1+2*j);
   end
   function EE5=E5(j)
  EE5=-Vhf2/(3*j*(1+j));
  end
   function EE7=E7(j)
  EE7=-Vhf2*(sqrt(1+1/j))/(1+2*j);
   end
  function EE9=E9(j)
  EE9=-Vhf2/(3*j)+Vhf2/(1+2*j);
   end
   function FF4=F4(j)
  FF4=-Vhf2*((1+j)^2)/(3*j*(sqrt(4*(j^2)-1)));
   end
   function FF6=F6(j)
  FF6=Vhf2*(j^2-1)/sqrt(j*(j-1)*(2*j-1)*(2*j+1));
  end
   function FF8=F8(j)
  FF8=Vhf2*(sqrt(j^2-1))/(3*j);
   end
   function GG2=G2(j)
  GG2=Vhf2*((1+j)*sqrt(1-2/(1+2*j)))/(3*j);
   end
   function GG4=G4(j)
  GG4=Vhf2*(sqrt((j-1)/j)*(1+j))/(1+2*j);
   end
   function GG6=G6(j)
  GG6=-Vhf2*(sqrt((j-1)/(j+1))*(j^2-1))/(3*j*(1+2*j));
  end
   function II1=I1(j)
  II1=Vhf2*(1-j^2)/(3*j-6*j^2);
   end
  function II3=I3(j)
  II3=-Vhf2*((j-1)^(3/2))/(sqrt(j)*(2*j-1));
   end
   function II5=I5(j)
  II5=Vhf2/(3*j);
   end
   function II7=I7(j)
  II7=Vhf2*sqrt((j-1)*j)/(2*j-1);
  end
   function II9=I9(j)
  II9=Vhf2*(2-j)/(-3+6*j);
  end