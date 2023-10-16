function [E,Eref,multiplicity,status]=shootForEigenvalue(system,md,E1,E2,k)
%implements algorithm 1: shooting for the eigenvalue with index k

status=1;
matchingIndex=md.matchingIndex;

E=(E1+E2)/2; %starting trial value for E
it=0;
while true
  [M,Me] = computeMismatch(system,md,E,matchingIndex);
  newton = -M/Me;
  E0 = E;
  E = E + newton;
  if abs(newton) <= md.tol || abs(E-E0) < md.tol 
     break;  %eigenvalue approximation in E is sufficiently accurate
  end
  if abs(E1-E2)< md.tol
      E=E0;
      status=-1;
      break;
  end
  if E>E2 || E<E1     
      if getIndex(system,md.referenceMesh,E0)<=k
          E1=E0;
      else
          E2=E0 ;
      end
      E=(E1+E2)/2;
  end
  if it > 30 && mod(it,5)==0
      if matchingIndex<3*length(md.h)/4
        matchingIndex=matchingIndex+ceil((length(md.h)-matchingIndex)/10);
      else
        matchingIndex=matchingIndex-ceil((matchingIndex)/10);
      end
  end
  if it > 50
      %retry using the smallest eigenvalue of the matching matrix as
      %mismatch
      [E,Eref,multiplicity,status]=shootForEigenvalue2(system,md,E,E0,k); 
      status = -2;  %too many iterations
      return;
  end
  it=it+1;
end


%computation of reference ( = more accurate) eigenvalue (is used for error estimation)
Eref=E;
it=0;
tol=md.tol/10;
while true
  [M,Me,X] = computeMismatch(system,md.referenceMesh,Eref,matchingIndex*2);
  newton = -M/Me;
  E0 = Eref;
  Eref = Eref + newton;
  if abs(newton) <= tol || abs(Eref-E0) < tol 
     break;
  end
  if it> 50
      status = -2;
      break;
  end
  it=it+1;
end
% determine multiplicity of the eigenvalue:
tmp=abs(eig(X));
m=min(tmp);
multiplicity=max(length(tmp(abs((tmp-m)/m)<1e-14)),1);
end


function varargout = computeMismatch(system,md,E,matchingIndex)
PsiL=system.A2/system.A1;
PsiR=system.B2/system.B1;
neq = size(system.A1,1);
if nargout >= 2
   PsiLe=zeros(neq,neq);
   PsiRe=zeros(neq,neq);
   [Tm,Te]=getTransferMatrices(E,md,6); 
   %left propagation:
   for k=1:matchingIndex
     u=Tm(1:neq,1:neq,k);   v=Tm(1:neq,neq+1:neq*2,k);
     up=Tm(neq+1:neq*2,1:neq,k); vp=Tm(neq+1:neq*2,neq+1:neq*2,k);
     ue=Te(1:neq,1:neq,k);   ve=Te(1:neq,neq+1:neq*2,k);
     upe=Te(neq+1:neq*2,1:neq,k); vpe=Te(neq+1:neq*2,neq+1:neq*2,k);
     Psi0=PsiL;
     PsiL=(u*PsiL+v)/(up*PsiL+vp); %CA propagation of Psi
     PsiLe=(ve+ue*Psi0+u*PsiLe)/(up*Psi0+vp)-PsiL*(vpe+upe*Psi0+up*PsiLe)/(up*Psi0+vp);
   end
   %right propagation:
   for k=length(md.h):-1:matchingIndex+1
     u=Tm(1:neq,1:neq,k);   v=Tm(1:neq,neq+1:neq*2,k);
     up=Tm(neq+1:neq*2,1:neq,k);   vp=Tm(neq+1:neq*2,neq+1:neq*2,k);
     ue=Te(1:neq,1:neq,k);   ve=Te(1:neq,neq+1:neq*2,k);
     upe=Te(neq+1:neq*2,1:neq,k); vpe=Te(neq+1:neq*2,neq+1:neq*2,k);
     Psi0=PsiR;
     PsiR=(vp'*PsiR-v')/(-up'*PsiR+u');
     PsiRe=(vpe'*Psi0+vp'*PsiRe-ve')/(-up'*Psi0+u')-PsiR*(-upe'*Psi0-up'*PsiRe+ue')/(-up'*Psi0+u'); 
   end
else
   [Tm]=getTransferMatrices(E,md,6);  
   for k=1:matchingIndex
     u=Tm(1:neq,1:neq,k);   v=Tm(1:neq,neq+1:neq*2,k);
     up=Tm(neq+1:neq*2,1:neq,k); vp=Tm(neq+1:neq*2,neq+1:neq*2,k);
     PsiL=(u*PsiL+v)/(up*PsiL+vp);
   end
   for k=length(md.h):-1:matchingIndex+1
     u=Tm(1:neq,1:neq,k);   v=Tm(1:neq,neq+1:neq*2,k);
     up=Tm(neq+1:neq*2,1:neq,k);   vp=Tm(neq+1:neq*2,neq+1:neq*2,k);
     PsiR=(vp'*PsiR-v')/(-up'*PsiR+u');
   end
end
X=PsiL-PsiR;

%compute derivative wrt E (eigenvalue):
M=det(X);
varargout{1}=M;
if nargout>=2
    Xe=(PsiLe-PsiRe);
   % Me=dete(X,Xe)
    Me=0;
    for i=1:neq
      for j=1:neq
         Xt=X;
         Xt(i,:)=zeros(neq,1);
         Xt(:,j)=zeros(1,neq);
         Xt(i,j)=Xe(i,j);
         Me=Me+det(Xt);
      end
    end
    varargout{2}=Me;
end
if nargout ==3
   varargout{3}=X;
end
end

% 
% function dphi=dete(X,Xe)
% neq=size(X,1);
% if neq ==1
%     dphi = Xe;
%     return;
% end
% dphi=0;
% je=1;
% for ie=1:neq
%         Xt=[X(1:ie-1,2:end);X(ie+1:end,2:end)];
%         Xet=[Xe(1:ie-1,2:end);Xe(ie+1:end,2:end)];
%         s1=Xe(ie,je)*(-1)^(ie+je)*det(Xt);
%         s2=X(ie,je)*(-1)^(ie+je)*dete(Xt,Xet);              
%         dphi=dphi+s1+s2;
% end
% end


function [E,Eref,multiplicity,status]=shootForEigenvalue2(system,md,E,Eprev,k)
%Shooting for the eigenvalue with index k. now the smallest eigenvalue of
%the matching matrix is used as mismatch
%This method is called when the Newton-Raphson process applied on the
%det-matching function seems to fail or slowly converge (e.g. when the
%problem has an eigenvalue with an even multiplicity )

status=1;
matchingIndex=md.matchingIndex;

Mprev=computeMismatch2(system,md,Eprev,matchingIndex);
Es=E;
it=0;
while true
  M = computeMismatch2(system,md,E,matchingIndex);
  if abs(E-Eprev)<=md.tol
     break;  %eigenvalue approximation in E is sufficiently accurate
  end
  if it > 30 && mod(it,5)==0
      matchingIndex=matchingIndex+ceil((length(md.h)-matchingIndex)/10);
  end
  if it > 50
      status = -2;  %too many iterations
      break;
  end
  it=it+1;
  Eprev0=E;
  E= E-M*(E-Eprev)/(M-Mprev);%secant method
  Mprev = M;
  Eprev=Eprev0;
end

%computation of reference ( = more accurate) eigenvalue (is used for error estimation)
Eref=E;
Eprev=Es;
if abs(Eprev-Eref)<md.tol
    Eprev=Eref-md.tol;
end
Mprev=computeMismatch2(system,md.referenceMesh,Eprev,matchingIndex*2);
it=0;
tol=md.tol/10; %tol/10 to  make sure that reference eigenvalue is more accurate
while true
  [M,multiplicity] = computeMismatch2(system,md.referenceMesh,Eref,matchingIndex*2);
  Eprev0=Eref;
  Eref= Eref-M*(Eref-Eprev)/(M-Mprev);
  Mprev = M;
  Eprev=Eprev0;
  if abs(Eprev-Eref) <= tol 
      break;
  end
  if it> 50
      status = -2;
      break;
  end
  it=it+1;
end
end


function varargout = computeMismatch2(system,md,E,matchingIndex)
PsiL=system.A2/system.A1;
PsiR=system.B2/system.B1;
neq = size(system.A1,1);
if nargout >= 2
   PsiLe=zeros(neq,neq);
   PsiRe=zeros(neq,neq);
   [Tm,Te]=getTransferMatrices(E,md,6); 
   for k=1:matchingIndex
     u=Tm(1:neq,1:neq,k);   v=Tm(1:neq,neq+1:neq*2,k);
     up=Tm(neq+1:neq*2,1:neq,k); vp=Tm(neq+1:neq*2,neq+1:neq*2,k);
     ue=Te(1:neq,1:neq,k);   ve=Te(1:neq,neq+1:neq*2,k);
     upe=Te(neq+1:neq*2,1:neq,k); vpe=Te(neq+1:neq*2,neq+1:neq*2,k);
     Psi0=PsiL;
     PsiL=(u*PsiL+v)/(up*PsiL+vp);
     PsiLe=(ve+ue*Psi0+u*PsiLe)/(up*Psi0+vp)-PsiL*(vpe+upe*Psi0+up*PsiLe)/(up*Psi0+vp);
   end
   for k=length(md.h):-1:matchingIndex+1
     u=Tm(1:neq,1:neq,k);   v=Tm(1:neq,neq+1:neq*2,k);
     up=Tm(neq+1:neq*2,1:neq,k);   vp=Tm(neq+1:neq*2,neq+1:neq*2,k);
     ue=Te(1:neq,1:neq,k);   ve=Te(1:neq,neq+1:neq*2,k);
     upe=Te(neq+1:neq*2,1:neq,k); vpe=Te(neq+1:neq*2,neq+1:neq*2,k);
     Psi0=PsiR;
     PsiR=(vp'*PsiR-v')/(-up'*PsiR+u');
     PsiRe=(vpe'*Psi0+vp'*PsiRe-ve')/(-up'*Psi0+u')-PsiR*(-upe'*Psi0-up'*PsiRe+ue')/(-up'*Psi0+u'); 
   end
else
   [Tm]=getTransferMatrices(E,md,6);  
   for k=1:matchingIndex
     u=Tm(1:neq,1:neq,k);   v=Tm(1:neq,neq+1:neq*2,k);
     up=Tm(neq+1:neq*2,1:neq,k); vp=Tm(neq+1:neq*2,neq+1:neq*2,k);
     PsiL=(u*PsiL+v)/(up*PsiL+vp);
   end
   for k=length(md.h):-1:matchingIndex+1
     u=Tm(1:neq,1:neq,k);   v=Tm(1:neq,neq+1:neq*2,k);
     up=Tm(neq+1:neq*2,1:neq,k);   vp=Tm(neq+1:neq*2,neq+1:neq*2,k);
     PsiR=(vp'*PsiR-v')/(-up'*PsiR+u');
   end
end
X=PsiL-PsiR;
D=eig(X);
tmp=abs(D);
[M,i]=min(tmp);
M=D(i);

varargout{1}=M;
if nargout>=2
    % determine multiplicity of the eigenvalue:
    Ma=abs(M);
    varargout{2}=max(length(tmp(abs((tmp-Ma)/Ma)<1e-14)),1);
end
end
