function [index,status] = getIndex(system,md,E)
%computes the index of the eigenvalue E (matrix Prufer procedure)

% see (V. Ledoux, M. Van Daele, Automatic computation of quantum-mechanical
% bound states and wavefunctions) for details on the procedure.

tol=md.tol;
PsiL=system.A2/system.A1;  %initial values
PsiR=system.B2/system.B1;
neq = size(system.A1,1);   %size of the system
Tm=getTransferMatrices(E,md,2,true); %second order CA method is used
%Propagation from left to right
phi=atan(eig(PsiL)); phi(phi<0)=phi(phi<0)+pi;
SPsiL= sum(phi);  %S(\Psi_L)
status=1;
for k=1:length(md.h)
    u=Tm(1:neq,1:neq,k);   v=Tm(1:neq,neq+1:neq*2,k);
    theta=atan(diag(v/u)); 
    D=md.dimat(:,:,k); 
    up=Tm(neq+1:neq*2,1:neq,k); vp=Tm(neq+1:neq*2,neq+1:neq*2,k);
    Psi0=D'*PsiL*D;
    PsiL=(u*Psi0+v)/(up*Psi0+vp);
    rho=atan(eig(PsiL)); rho(rho<0)=rho(rho<0)+pi;
    Omega=eye(neq)/(u-v*Psi0)*(u*Psi0+v);
    tau=atan(eig(Omega)); tau(tau<0)=tau(tau<0)+pi;
    theta(theta<0)=theta(theta<0)+pi;
    if any(theta(theta<tol)) 
        status =-3;         %a zero eigenvalue may have been missed due to the accuracy of the computations being related to tol 
        warning('The Prufer indexing function may have produced an inaccurate result. Consider decreasing the input tolerance.');
    end
    vt=E-md.v0(:,k);
    d=sqrt(vt(vt>0))*md.h(k);   
    sum_sk= sum(floor(d/pi));
    sum_theta= sum(theta)+sum_sk*pi; 
    SPsiL= SPsiL + sum_theta+ sum(rho-tau);
    PsiL=D*PsiL*D';
end
phi=atan(eig(PsiR)); phi(phi<=0)=phi(phi<=0)+pi;
SPsiR= sum(phi);
w= atan(eig((PsiL-PsiR)/(eye(neq)+PsiL*PsiR))); 
w(w<0)=w(w<0)+pi;
index=round((SPsiL-SPsiR-sum(w))/pi+neq);

