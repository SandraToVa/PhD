(* ::Package:: *)

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

k=0;
figure
hold on
plot(x,Y0(1,:,1),'b')
plot(x,Y0(2,:,1),'r')
plot(x,Y0(3,:,1),'m')
plot(x,Y0(4,:,1),'g')
hold off
xlabel('x')
E=EigvData.eigenvalues(k+1)/m;
title(E)
legend('S (x)','H+(x)','H-(x)','H0 (x)')

k=1;
figure
hold on
plot(x,Y1(1,:,1),'b')
plot(x,Y1(2,:,1),'r')
plot(x,Y1(3,:,1),'m')
plot(x,Y1(4,:,1),'g')
hold off
xlabel('x')
E=EigvData.eigenvalues(k+1)/m;
title(E)
legend('S (x)','H+(x)','H-(x)','H0 (x)')

k=2;
figure
hold on
plot(x,Y2(1,:,1),'b')
plot(x,Y2(2,:,1),'r')
plot(x,Y2(3,:,1),'m')
plot(x,Y2(4,:,1),'g')
hold off
xlabel('x')
E=EigvData.eigenvalues(k+1)/m;
title(E)
legend('S (x)','H+(x)','H-(x)','H0 (x)')

k=3;
figure
hold on
plot(x,Y3(1,:,1),'b')
plot(x,Y3(2,:,1),'r')
plot(x,Y3(3,:,1),'m')
plot(x,Y3(4,:,1),'g')
hold off
xlabel('x')
E=EigvData.eigenvalues(k+1)/m;
title(E)
legend('S (x)','H+(x)','H-(x)','H0 (x)')

k=4;
figure
hold on
plot(x,Y4(1,:,1),'b')
plot(x,Y4(2,:,1),'r')
plot(x,Y4(3,:,1),'m')
plot(x,Y4(4,:,1),'g')
hold off
xlabel('x')
E=EigvData.eigenvalues(k+1)/m;
title(E)
legend('S (x)','H+(x)','H-(x)','H0 (x)')

k=5;
figure
hold on
plot(x,Y5(1,:,1),'b')
plot(x,Y5(2,:,1),'r')
plot(x,Y5(3,:,1),'m')
plot(x,Y5(4,:,1),'g')
hold off
xlabel('x')
E=EigvData.eigenvalues(k+1)/m;
title(E)
legend('S (x)','H+(x)','H-(x)','H0 (x)')

k=6;
figure
hold on
plot(x,Y6(1,:,1),'b')
plot(x,Y6(2,:,1),'r')
plot(x,Y6(3,:,1),'m')
plot(x,Y6(4,:,1),'g')
hold off
xlabel('x')
E=EigvData.eigenvalues(k+1)/m;
title(E)
legend('S (x)','H+(x)','H-(x)','H0 (x)')

k=7;
figure
hold on
plot(x,Y7(1,:,1),'b')
plot(x,Y7(2,:,1),'r')
plot(x,Y7(3,:,1),'m')
plot(x,Y7(4,:,1),'g')
hold off
xlabel('x')
E=EigvData.eigenvalues(k+1)/m;
title(E)
legend('S (x)','H+(x)','H-(x)','H0 (x)')

k=8;
figure
hold on
plot(x,Y8(1,:,1),'b')
plot(x,Y8(2,:,1),'r')
plot(x,Y8(3,:,1),'m')
plot(x,Y8(4,:,1),'g')
hold off
xlabel('x')
E=EigvData.eigenvalues(k+1)/m;
title(E)
legend('S (x)','H+(x)','H-(x)','H0 (x)')

k=9;
figure
hold on
plot(x,Y9(1,:,1),'b')
plot(x,Y9(2,:,1),'r')
plot(x,Y9(3,:,1),'m')
plot(x,Y9(4,:,1),'g')
hold off
xlabel('x')
E=EigvData.eigenvalues(k+1)/m;
title(E)
legend('S (x)','H+(x)','H-(x)','H0 (x)')

k=10;
figure
hold on
plot(x,Y10(1,:,1),'b')
plot(x,Y10(2,:,1),'r')
plot(x,Y10(3,:,1),'m')
plot(x,Y10(4,:,1),'g')
hold off
xlabel('x')
E=EigvData.eigenvalues(k+1)/m;
title(E)
legend('S (x)','H+(x)','H-(x)','H0 (x)')

