
function COUPj11
% 6x6 coupled shcroedinger equation. Autor: Ruben Oncala
% Marletta's testexample for the SL12F code.
disp('**************************************************')
disp('Hybrid-Quarkonium meson states COUPj11:')
% disp('  Marletta''s testexample for the SL12F code')
disp('**************************************************')
format long;


  % the endpoints of the integration interval:
system.a=0.001;   
system.b=22; 
% parameters of the boundary conditions:
system.A1= eye(6);
system.A2= zeros(6,6);
system.B1= eye(6);
system.B2= zeros(6,6);
% function handle to the function returning the potential matrix
system.V=@potentialMatrix;


[m]=parameters
[j]=parameters1
[lam0]=parameters2

disp(' ');
disp('Calculated eigenvalues with indices between 0 and 10:');
t=cputime;
[EigvData,meshData]=computeEigenvalues(system,0,15,5e-8);

% disp(['Number of intervals in the mesh: ' num2str(length(meshData.h))]); % number of intervals in the mesh
% disp(['Time: ' num2str (cputime - t) ' sec']);
% disp(sprintf('k \t  E_k \t\t\t\t ErrEst \t status'))
for i=1:length(EigvData.eigenvalues)
    disp(sprintf('%-3.0f\t %-16.12f\t %-2.0e\t\t % d', EigvData.indices(i),EigvData.eigenvalues(i)/m,EigvData.errorEstimations(i),EigvData.status(i)));
end


[x,Y0,YP0]=computeEigenfunction(system,meshData,EigvData.eigenvalues(1),1);
[x,Y1,YP1]=computeEigenfunction(system,meshData,EigvData.eigenvalues(2),1);
[x,Y2,YP2]=computeEigenfunction(system,meshData,EigvData.eigenvalues(3),1);
[x,Y3,YP3]=computeEigenfunction(system,meshData,EigvData.eigenvalues(4),1);
[x,Y4,YP4]=computeEigenfunction(system,meshData,EigvData.eigenvalues(5),1);
[x,Y5,YP5]=computeEigenfunction(system,meshData,EigvData.eigenvalues(6),1);
[x,Y6,YP6]=computeEigenfunction(system,meshData,EigvData.eigenvalues(7),1);
[x,Y7,YP7]=computeEigenfunction(system,meshData,EigvData.eigenvalues(8),1);
[x,Y8,YP8]=computeEigenfunction(system,meshData,EigvData.eigenvalues(9),1);
[x,Y9,YP9]=computeEigenfunction(system,meshData,EigvData.eigenvalues(10),1);
[x,Y10,YP10]=computeEigenfunction(system,meshData,EigvData.eigenvalues(11),1);
[x,Y11,YP11]=computeEigenfunction(system,meshData,EigvData.eigenvalues(12),1);
[x,Y12,YP12]=computeEigenfunction(system,meshData,EigvData.eigenvalues(13),1);
[x,Y13,YP13]=computeEigenfunction(system,meshData,EigvData.eigenvalues(14),1);
[x,Y14,YP14]=computeEigenfunction(system,meshData,EigvData.eigenvalues(15),1);
[x,Y15,YP15]=computeEigenfunction(system,meshData,EigvData.eigenvalues(16),1);


% calculem la proporcio de hybrid
100*(trapz (x,Y0(1,:,1))^2+trapz (x,Y0(2,:,1))^2+trapz (x,Y0(3,:,1))^2+trapz (x,Y0(4,:,1))^2+trapz (x,Y0(5,:,1))^2)/(trapz (x,Y0(1,:,1))^2+trapz (x,Y0(2,:,1))^2+trapz (x,Y0(3,:,1))^2+trapz (x,Y0(4,:,1))^2+trapz (x,Y0(5,:,1))^2+trapz (x,Y0(6,:,1))^2)
100*(trapz (x,Y1(1,:,1))^2+trapz (x,Y1(2,:,1))^2+trapz (x,Y1(3,:,1))^2+trapz (x,Y1(4,:,1))^2+trapz (x,Y1(5,:,1))^2)/(trapz (x,Y1(1,:,1))^2+trapz (x,Y1(2,:,1))^2+trapz (x,Y1(3,:,1))^2+trapz (x,Y1(4,:,1))^2+trapz (x,Y1(5,:,1))^2+trapz (x,Y1(6,:,1))^2)
100*(trapz (x,Y2(1,:,1))^2+trapz (x,Y2(2,:,1))^2+trapz (x,Y2(3,:,1))^2+trapz (x,Y2(4,:,1))^2+trapz (x,Y2(5,:,1))^2)/(trapz (x,Y2(1,:,1))^2+trapz (x,Y2(2,:,1))^2+trapz (x,Y2(3,:,1))^2+trapz (x,Y2(4,:,1))^2+trapz (x,Y2(5,:,1))^2+trapz (x,Y2(6,:,1))^2)
100*(trapz (x,Y3(1,:,1))^2+trapz (x,Y3(2,:,1))^2+trapz (x,Y3(3,:,1))^2+trapz (x,Y3(4,:,1))^2+trapz (x,Y3(5,:,1))^2)/(trapz (x,Y3(1,:,1))^2+trapz (x,Y3(2,:,1))^2+trapz (x,Y3(3,:,1))^2+trapz (x,Y3(4,:,1))^2+trapz (x,Y3(5,:,1))^2+trapz (x,Y3(6,:,1))^2)
100*(trapz (x,Y4(1,:,1))^2+trapz (x,Y4(2,:,1))^2+trapz (x,Y4(3,:,1))^2+trapz (x,Y4(4,:,1))^2+trapz (x,Y4(5,:,1))^2)/(trapz (x,Y4(1,:,1))^2+trapz (x,Y4(2,:,1))^2+trapz (x,Y4(3,:,1))^2+trapz (x,Y4(4,:,1))^2+trapz (x,Y4(5,:,1))^2+trapz (x,Y4(6,:,1))^2)
100*(trapz (x,Y5(1,:,1))^2+trapz (x,Y5(2,:,1))^2+trapz (x,Y5(3,:,1))^2+trapz (x,Y5(4,:,1))^2+trapz (x,Y5(5,:,1))^2)/(trapz (x,Y5(1,:,1))^2+trapz (x,Y5(2,:,1))^2+trapz (x,Y5(3,:,1))^2+trapz (x,Y5(4,:,1))^2+trapz (x,Y5(5,:,1))^2+trapz (x,Y5(6,:,1))^2)
100*(trapz (x,Y6(1,:,1))^2+trapz (x,Y6(2,:,1))^2+trapz (x,Y6(3,:,1))^2+trapz (x,Y6(4,:,1))^2+trapz (x,Y6(5,:,1))^2)/(trapz (x,Y6(1,:,1))^2+trapz (x,Y6(2,:,1))^2+trapz (x,Y6(3,:,1))^2+trapz (x,Y6(4,:,1))^2+trapz (x,Y6(5,:,1))^2+trapz (x,Y6(6,:,1))^2)
100*(trapz (x,Y7(1,:,1))^2+trapz (x,Y7(2,:,1))^2+trapz (x,Y7(3,:,1))^2+trapz (x,Y7(4,:,1))^2+trapz (x,Y7(5,:,1))^2)/(trapz (x,Y7(1,:,1))^2+trapz (x,Y7(2,:,1))^2+trapz (x,Y7(3,:,1))^2+trapz (x,Y7(4,:,1))^2+trapz (x,Y7(5,:,1))^2+trapz (x,Y7(6,:,1))^2)
100*(trapz (x,Y8(1,:,1))^2+trapz (x,Y8(2,:,1))^2+trapz (x,Y8(3,:,1))^2+trapz (x,Y8(4,:,1))^2+trapz (x,Y8(5,:,1))^2)/(trapz (x,Y8(1,:,1))^2+trapz (x,Y8(2,:,1))^2+trapz (x,Y8(3,:,1))^2+trapz (x,Y8(4,:,1))^2+trapz (x,Y8(5,:,1))^2+trapz (x,Y8(6,:,1))^2)
100*(trapz (x,Y9(1,:,1))^2+trapz (x,Y9(2,:,1))^2+trapz (x,Y9(3,:,1))^2+trapz (x,Y9(4,:,1))^2+trapz (x,Y9(5,:,1))^2)/(trapz (x,Y9(1,:,1))^2+trapz (x,Y9(2,:,1))^2+trapz (x,Y9(3,:,1))^2+trapz (x,Y9(4,:,1))^2+trapz (x,Y9(5,:,1))^2+trapz (x,Y9(6,:,1))^2)
100*(trapz (x,Y10(1,:,1))^2+trapz (x,Y10(2,:,1))^2+trapz (x,Y10(3,:,1))^2+trapz (x,Y10(4,:,1))^2+trapz (x,Y10(5,:,1))^2)/(trapz (x,Y10(1,:,1))^2+trapz (x,Y10(2,:,1))^2+trapz (x,Y10(3,:,1))^2+trapz (x,Y10(4,:,1))^2+trapz (x,Y10(5,:,1))^2+trapz (x,Y10(6,:,1))^2)
100*(trapz (x,Y11(1,:,1))^2+trapz (x,Y11(2,:,1))^2+trapz (x,Y11(3,:,1))^2+trapz (x,Y11(4,:,1))^2+trapz (x,Y11(5,:,1))^2)/(trapz (x,Y11(1,:,1))^2+trapz (x,Y11(2,:,1))^2+trapz (x,Y11(3,:,1))^2+trapz (x,Y11(4,:,1))^2+trapz (x,Y11(5,:,1))^2+trapz (x,Y11(6,:,1))^2)
100*(trapz (x,Y12(1,:,1))^2+trapz (x,Y12(2,:,1))^2+trapz (x,Y12(3,:,1))^2+trapz (x,Y12(4,:,1))^2+trapz (x,Y12(5,:,1))^2)/(trapz (x,Y12(1,:,1))^2+trapz (x,Y12(2,:,1))^2+trapz (x,Y12(3,:,1))^2+trapz (x,Y12(4,:,1))^2+trapz (x,Y12(5,:,1))^2+trapz (x,Y12(6,:,1))^2)
100*(trapz (x,Y13(1,:,1))^2+trapz (x,Y13(2,:,1))^2+trapz (x,Y13(3,:,1))^2+trapz (x,Y13(4,:,1))^2+trapz (x,Y13(5,:,1))^2)/(trapz (x,Y13(1,:,1))^2+trapz (x,Y13(2,:,1))^2+trapz (x,Y13(3,:,1))^2+trapz (x,Y13(4,:,1))^2+trapz (x,Y13(5,:,1))^2+trapz (x,Y13(6,:,1))^2)
100*(trapz (x,Y14(1,:,1))^2+trapz (x,Y14(2,:,1))^2+trapz (x,Y14(3,:,1))^2+trapz (x,Y14(4,:,1))^2+trapz (x,Y14(5,:,1))^2)/(trapz (x,Y14(1,:,1))^2+trapz (x,Y14(2,:,1))^2+trapz (x,Y14(3,:,1))^2+trapz (x,Y14(4,:,1))^2+trapz (x,Y14(5,:,1))^2+trapz (x,Y14(6,:,1))^2)
100*(trapz (x,Y15(1,:,1))^2+trapz (x,Y15(2,:,1))^2+trapz (x,Y15(3,:,1))^2+trapz (x,Y15(4,:,1))^2+trapz (x,Y15(5,:,1))^2)/(trapz (x,Y15(1,:,1))^2+trapz (x,Y15(2,:,1))^2+trapz (x,Y15(3,:,1))^2+trapz (x,Y15(4,:,1))^2+trapz (x,Y15(5,:,1))^2+trapz (x,Y15(6,:,1))^2)


k=2;
figure
hold on
plot(x,Y2(1,:,1),'b-*')
plot(x,Y2(2,:,1),'r-+')
plot(x,Y2(3,:,1),'m-o')
plot(x,Y2(4,:,1),'g.')
plot(x,Y2(5,:,1),'y-x')
plot(x,Y2(6,:,1),'k')
hold off
xlabel('r')
E=EigvData.eigenvalues(k+1)/m;
title(E)
legend('P_1^{--}(r)','P_1^{+-}(r)','P_1^{-+}(r)','P_1^{++}(r)','P_1^{00}(r)','S_0(r)')
%{
k=1;
figure
hold on
plot(x,Y1(1,:,1),'b')
plot(x,Y1(2,:,1),'r')
plot(x,Y1(3,:,1),'m')
plot(x,Y1(4,:,1),'g')
plot(x,Y1(5,:,1),'y')
plot(x,Y1(6,:,1),'k')
hold off
xlabel('x')
E=EigvData.eigenvalues(k+1)/m;
title(E)
legend('H-- (x)','H+-(x)','H-+(x)','H++(x)','H00 (x)','S00 (x)')

k=2;
figure
hold on
plot(x,Y2(1,:,1),'b')
plot(x,Y2(2,:,1),'r')
plot(x,Y2(3,:,1),'m')
plot(x,Y2(4,:,1),'g')
plot(x,Y2(5,:,1),'y')
plot(x,Y2(6,:,1),'k')
hold off
xlabel('x')
E=EigvData.eigenvalues(k+1)/m;
title(E)
legend('H-- (x)','H+-(x)','H-+(x)','H++(x)','H00 (x)','S00 (x)')

k=3;
figure
hold on
plot(x,Y3(1,:,1),'b')
plot(x,Y3(2,:,1),'r')
plot(x,Y3(3,:,1),'m')
plot(x,Y3(4,:,1),'g')
plot(x,Y3(5,:,1),'y')
plot(x,Y3(6,:,1),'k')
hold off
xlabel('x')
E=EigvData.eigenvalues(k+1)/m;
title(E)
legend('H-- (x)','H+-(x)','H-+(x)','H++(x)','H00 (x)','S00 (x)')

k=4;
figure
hold on
plot(x,Y4(1,:,1),'b')
plot(x,Y4(2,:,1),'r')
plot(x,Y4(3,:,1),'m')
plot(x,Y4(4,:,1),'g')
plot(x,Y4(5,:,1),'y')
plot(x,Y4(6,:,1),'k')
hold off
xlabel('x')
E=EigvData.eigenvalues(k+1)/m;
title(E)
legend('H-- (x)','H+-(x)','H-+(x)','H++(x)','H00 (x)','S00 (x)')

k=5;
figure
hold on
plot(x,Y5(1,:,1),'b')
plot(x,Y5(2,:,1),'r')
plot(x,Y5(3,:,1),'m')
plot(x,Y5(4,:,1),'g')
plot(x,Y5(5,:,1),'y')
plot(x,Y5(6,:,1),'k')
hold off
xlabel('x')
E=EigvData.eigenvalues(k+1)/m;
title(E)
legend('H-- (x)','H+-(x)','H-+(x)','H++(x)','H00 (x)','S00 (x)')

k=6;
figure
hold on
plot(x,Y6(1,:,1),'b')
plot(x,Y6(2,:,1),'r')
plot(x,Y6(3,:,1),'m')
plot(x,Y6(4,:,1),'g')
plot(x,Y6(5,:,1),'y')
plot(x,Y6(6,:,1),'k')
hold off
xlabel('x')
E=EigvData.eigenvalues(k+1)/m;
title(E)
legend('H-- (x)','H+-(x)','H-+(x)','H++(x)','H00 (x)','S00 (x)')

k=7;
figure
hold on
plot(x,Y7(1,:,1),'b')
plot(x,Y7(2,:,1),'r')
plot(x,Y7(3,:,1),'m')
plot(x,Y7(4,:,1),'g')
plot(x,Y7(5,:,1),'y')
plot(x,Y7(6,:,1),'k')
hold off
xlabel('x')
E=EigvData.eigenvalues(k+1)/m;
title(E)
legend('H-- (x)','H+-(x)','H-+(x)','H++(x)','H00 (x)','S00 (x)')

k=8;
figure
hold on
plot(x,Y8(1,:,1),'b')
plot(x,Y8(2,:,1),'r')
plot(x,Y8(3,:,1),'m')
plot(x,Y8(4,:,1),'g')
plot(x,Y8(5,:,1),'y')
plot(x,Y8(6,:,1),'k')
hold off
xlabel('x')
E=EigvData.eigenvalues(k+1)/m;
title(E)
legend('H-- (x)','H+-(x)','H-+(x)','H++(x)','H00 (x)','S00 (x)')

k=9;
figure
hold on
plot(x,Y9(1,:,1),'b')
plot(x,Y9(2,:,1),'r')
plot(x,Y9(3,:,1),'m')
plot(x,Y9(4,:,1),'g')
plot(x,Y9(5,:,1),'y')
plot(x,Y9(6,:,1),'k')
hold off
xlabel('x')
E=EigvData.eigenvalues(k+1)/m;
title(E)
legend('H-- (x)','H+-(x)','H-+(x)','H++(x)','H00 (x)','S00 (x)')

k=10;
figure
hold on
plot(x,Y10(1,:,1),'b')
plot(x,Y10(2,:,1),'r')
plot(x,Y10(3,:,1),'m')
plot(x,Y10(4,:,1),'g')
plot(x,Y10(5,:,1),'y')
plot(x,Y10(6,:,1),'k')
hold off
xlabel('x')
E=EigvData.eigenvalues(k+1)/m;
title(E)
legend('H-- (x)','H+-(x)','H-+(x)','H++(x)','H00 (x)','S00 (x)')

k=11;
figure
hold on
plot(x,Y11(1,:,1),'b')
plot(x,Y11(2,:,1),'r')
plot(x,Y11(3,:,1),'m')
plot(x,Y11(4,:,1),'g')
plot(x,Y11(5,:,1),'y')
plot(x,Y11(6,:,1),'k')
hold off
xlabel('x')
E=EigvData.eigenvalues(k+1)/m;
title(E)
legend('H-- (x)','H+-(x)','H-+(x)','H++(x)','H00 (x)','S00 (x)')

k=12;
figure
hold on
plot(x,Y12(1,:,1),'b')
plot(x,Y12(2,:,1),'r')
plot(x,Y12(3,:,1),'m')
plot(x,Y12(4,:,1),'g')
plot(x,Y12(5,:,1),'y')
plot(x,Y12(6,:,1),'k')
hold off
xlabel('x')
E=EigvData.eigenvalues(k+1)/m;
title(E)
legend('H-- (x)','H+-(x)','H-+(x)','H++(x)','H00 (x)','S00 (x)')


k=13;
figure
hold on
plot(x,Y13(1,:,1),'b')
plot(x,Y13(2,:,1),'r')
plot(x,Y13(3,:,1),'m')
plot(x,Y13(4,:,1),'g')
plot(x,Y13(5,:,1),'y')
plot(x,Y13(6,:,1),'k')
hold off
xlabel('x')
E=EigvData.eigenvalues(k+1)/m;
title(E)
legend('H-- (x)','H+-(x)','H-+(x)','H++(x)','H00 (x)','S00 (x)')


k=14;
figure
hold on
plot(x,Y14(1,:,1),'b')
plot(x,Y14(2,:,1),'r')
plot(x,Y14(3,:,1),'m')
plot(x,Y14(4,:,1),'g')
plot(x,Y14(5,:,1),'y')
plot(x,Y14(6,:,1),'k')
hold off
xlabel('x')
E=EigvData.eigenvalues(k+1)/m;
title(E)
legend('H-- (x)','H+-(x)','H-+(x)','H++(x)','H00 (x)','S00 (x)')


k=15;
figure
hold on
plot(x,Y15(1,:,1),'b')
plot(x,Y15(2,:,1),'r')
plot(x,Y15(3,:,1),'m')
plot(x,Y15(4,:,1),'g')
plot(x,Y15(5,:,1),'y')
plot(x,Y15(6,:,1),'k')
hold off
xlabel('x')
E=EigvData.eigenvalues(k+1)/m;
title(E)
legend('H-- (x)','H+-(x)','H-+(x)','H++(x)','H00 (x)','S00 (x)')


%}







end

% aqui escollirem els parametres L i m del sistema

function [m]=parameters
 m=1.4702;
  %m=4.8802;
end

function [j]=parameters1
j=0;
end 

function [lam0,lam1,lam3,k,pi]=parameters2
lam0=0.6;
lam1=0.5*0.1182;
lam3=0.2299;
k=0.187;
pi=3.1416;
end
%---------------------------------------------------



function r=potentialMatrix(x) % returns the potential matrix evaluated in x

[j]=parameters1;

  for i=1:4 
  r(1,1,i)=A(x(i),j-1);
  r(2,2,i)=CC(x(i),j-1);
  r(3,3,i)=A(x(i),j+1);
  r(4,4,i)=CC(x(i),j+1);
  r(5,5,i)=F(x(i),j);
  r(6,6,i)=R(x(i),j);

  r(1,2,i)=B(x(i),j-1);
  r(2,1,i)=r(1,2,i);
  r(1,3,i)=0;
  r(3,1,i)=r(1,3,i);
  r(1,4,i)=0;
  r(4,1,i)=r(1,4,i);

  r(2,3,i)=0;
  r(3,2,i)=r(2,3,i);
  r(2,4,i)=0;
  r(4,2,i)=r(2,4,i);
  r(3,4,i)=B(x(i),j+1);
  r(4,3,i)=r(3,4,i);

  r(1,5,i)=0;
  r(5,1,i)=r(1,5,i);
  r(2,5,i)=0;
  r(5,2,i)=r(2,5,i);
  r(3,5,i)=0;
  r(5,3,i)=r(3,5,i);
  r(4,5,i)=0;
  r(5,4,i)=r(4,5,i);

  r(1,6,i)=2*alpha(j)*Vqs(x(i));
  r(6,1,i)=r(1,6,i);
  r(2,6,i)=2*omega(j)*Vps(x(i))+2*beta(j)*Vqs(x(i));
  r(6,1,i)=r(2,6,i);
  r(3,6,i)=2*eta(j)*Vps(x(i))+2*gamma(j)*Vqs(x(i));
  r(6,3,i)=r(3,6,i);
  r(4,6,i)=2*delta(j)*Vqs(x(i));
  r(6,4,i)=r(4,6,i);
  r(5,6,i)=2*Vps(x(i));
  r(6,5,i)=r(5,6,i);



  end
end

%--------------------------------------------
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
