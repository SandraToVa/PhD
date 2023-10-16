function [x,YM,YPM]=computeEigenfunction(system,md,E,multiplicity)
%computes the eigenfunctions and their first order derivatives corresponding
% to the eigenvalue E in the meshpoints of the mesh given by the input argument md.
%OUTPUT:    x : vector of x-values (=meshpoints) where the eigenfunctions were evaluated
%           YM: 3-dim matrix (n times nsteps times multiplicity)
%               contains information on the eigenfunctions                
%               YM(i,j,k) = i'th element (1<=i<=dimension n of the system) in
%               the k'th eigenfunction evaluated in the jth meshpoint .
%               Each YM(:,:,k) with ( 1<= k <= multiplicity ) represents an
%               eigenfunction.
%           YMP: similar for the first order derivatives

%NOTE: the eigenfunction and its derivative are only evaluated in the
%meshpoints (the values of these points are returned by the function 
%computeEigenfunction in the vector x). This may not be sufficient to produce a
%nice smooth plot of a higher eigenfunction. The evaluation of the
%eigenfunctions in user specified points is future work.

% see (V. Ledoux, M. Van Daele, Automatic computation of quantum-mechanical
% bound states and wavefunctions) for more details on the procedure.

addpath('./source')

%left hand part
PsiL=system.A2/system.A1;
Y(:,:,1) = system.A2;
neq = size(system.A1,1);
if multiplicity > neq
   error(SchrodSys:argChk, 'the multiplicity of an eigenvalue cannot be larger than the size of the system');   
end
PsiLe=zeros(neq,neq); %derivative w.r.t. E, used for normalization
YPL=system.A1;
YP(:,:,1)=YPL;
[Tm,Te]=getTransferMatrices(E,md,6); 
%Propagation of the left hand solution:
for k=1:md.matchingIndex 
    %entries of the transfer matrix (matrices U and W and their first
    %derivatives):
    u=Tm(1:neq,1:neq,k);   v=Tm(1:neq,neq+1:neq*2,k);
    up=Tm(neq+1:neq*2,1:neq,k); vp=Tm(neq+1:neq*2,neq+1:neq*2,k);
    %their derivatives w.r.t. E:
    ue=Te(1:neq,1:neq,k);   ve=Te(1:neq,neq+1:neq*2,k);
    upe=Te(neq+1:neq*2,1:neq,k); vpe=Te(neq+1:neq*2,neq+1:neq*2,k);
    %
    Psi0=PsiL;
    %CP propagation of Psi:
    PsiL=(u*PsiL+v)/(up*PsiL+vp);
    %propagation of the first derivative of Psi w.r.t. E:
    PsiLe=(ve+ue*Psi0+u*PsiLe)/(up*Psi0+vp)-PsiL*(vpe+upe*Psi0+up*PsiLe)/(up*Psi0+vp);
    %propagation of Y':
    YPL = (up*Psi0+vp)*YPL;
    Y(:,:,k+1)= PsiL*YPL; %wavefunction
    YP(:,:,k+1)=YPL;      %first order derivative  
end
% 
%right hand part
PsiR=system.B2/system.B1;
Y(:,:,length(md.h)+1) = system.B2;
PsiRe=zeros(neq,neq);
YPR=system.B1;
YP(:,:,length(md.h)+1)=YPL;
%Propagation of right hand solution
for k=length(md.h):-1:md.matchingIndex+1
    u=Tm(1:neq,1:neq,k);   v=Tm(1:neq,neq+1:neq*2,k);
    up=Tm(neq+1:neq*2,1:neq,k);   vp=Tm(neq+1:neq*2,neq+1:neq*2,k);
    ue=Te(1:neq,1:neq,k);   ve=Te(1:neq,neq+1:neq*2,k);
    upe=Te(neq+1:neq*2,1:neq,k); vpe=Te(neq+1:neq*2,neq+1:neq*2,k);
    Psi0=PsiR;
    PsiR=(vp'*PsiR-v')/(-up'*PsiR+u');
    PsiRe=(vpe'*Psi0+vp'*PsiRe-ve')/(-up'*Psi0+u')-PsiR*(-upe'*Psi0-up'*PsiRe+ue')/(-up'*Psi0+u'); 
    YPR = (-up'*Psi0+u')*YPR;
    Y(:,:,k)= PsiR*YPR;
    YP(:,:,k)=YPR;
end

%We have to make sure to match the left-hand half of an eigenfunction with
%the right-hand half of the same eigenfunction. A procedure is used which
%was discussed in [M.Marletta, Numer. Algorithms 4 (1993) 65].
[VL,D]=eig(YPR'*PsiL*YPL-(PsiR*YPR)'*YPL); 
[VR,D2]=eig((YPR'*PsiL*YPL-(PsiR*YPR)'*YPL)');
%determine the set of vectors v_L^j and v_R^j
d=abs(diag(D));
d2=abs(diag(D2));
for j=1:multiplicity
   [e,i ]=min(d);
   vl(:,j)=VL(:,i); 
   d(i)=max(d)+10;
   [e,i ]=min(d2);
   vr(:,j)=VR(:,i);
   d2(i)=max(d2)+10;
end
ZL=[PsiL*YPL*vl;YPL*vl];
ZR=[PsiR*YPR*vr;YPR*vr];
for j=1:multiplicity
  AT(:,j)=ZR\ZL(:,j); 
end
N1=YPL'*PsiLe*YPL;  %needed for normalization
N2=-YPR'*PsiRe*YPR;
YM=zeros(neq,length(md.h)+1,multiplicity);
YPM=zeros(neq,length(md.h)+1,multiplicity);
for j=1:multiplicity
  vrt=AT'*vr(:,j);
  vlt=vl(:,j);
  fnorm = sqrt(abs(vlt'*N1*vlt+vrt'*N2*vrt));
  for i=1:md.matchingIndex
    YM(:,i,j)=(Y(:,:,i)*vlt)./fnorm;
    YPM(:,i,j)=(YP(:,:,i)*vlt)./fnorm;
  end
  for i=md.matchingIndex+1:length(md.h)+1
    YM(:,i,j)=(Y(:,:,i)*vrt)./fnorm;
    YPM(:,i,j)=(YP(:,:,i)*vrt)./fnorm;
  end
end

x=[system.a system.a+cumsum(md.h)];

