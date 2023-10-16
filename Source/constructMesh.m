function md=constructMesh(system,tol)
% md=constructMesh(system,tol)
%  constructs a mesh corresponding to the user input tolerance tol
%  and computes some associated data (this data is energy-dependent.
%  Therefore it needs to be computed only once and is stored here 
%  (in the md-object) for repeated reuse)

tol = max(tol,1e-14);
md.tol=tol;
x = system.a;
h = (system.b-x)/100;
hmin= (system.b-system.a)/3; %3 is the minimum number of steps
hnew = 2*h;
i=1;
neq=size(system.A1,1);
while x < system.b
    it=0;
    while abs(hnew/h-1) > 0.01 && it< 7
      h=min([hnew,system.b-x,hmin]);  
      [v0,v,dimat] = computeVcoeffs(system.V,x,h);
      err = estimateError(v0,v,h);
      if err>eps
        hnew = h*(tol/err)^(1/5);
      else
        hnew= 5*h;  
      end
      if h == system.b-x && hnew>h
          break;
      end
      it = it+1;
    end  
    md.v1(:,:,i)=v(:,:,1); md.v2(:,:,i)=v(:,:,2);
    md.dimat(:,:,i) = dimat; md.v0(:,i) = v0;
    md.h(i) = h;
    x = x + h;
    hnew = h; h = h/2;
    i=i+1;
end

%it is better to compute md.cu,md.cv,... in a new loop, since the matrices can
%then be preallocated (this gives nice improvement in speed when neq is
%large)
md.cu=zeros(neq,neq,3,length(md.h));
md.cv=zeros(neq,neq,3,length(md.h));
md.cup=zeros(neq,neq,3,length(md.h));
md.cvp=zeros(neq,neq,3,length(md.h));
for i=1:length(md.h)
   [md.cu(:,:,:,i),md.cup(:,:,:,i),md.cv(:,:,:,i),md.cvp(:,:,:,i)]=computeCcoeffs(md.v0(:,i),md.v1(:,:,i),md.v2(:,:,i),md.h(i)); 
end



%a second mesh with meshpoints
%in between the meshpoints of the first mesh and twice so many meshpoints
hs = [0 md.h/2]+ [md.h/2 0];
x=system.a;
md.referenceMesh.tol=md.tol;
%preallocation is important when neq is large
md.referenceMesh.cu=zeros(neq,neq,3,length(hs)*2);
md.referenceMesh.cup=zeros(neq,neq,3,length(hs)*2);
md.referenceMesh.cv=zeros(neq,neq,3,length(hs)*2);
md.referenceMesh.cvp=zeros(neq,neq,3,length(hs)*2);
md.referenceMesh.v0=zeros(neq,length(hs)*2);
md.referenceMesh.v1=zeros(neq,neq,length(hs)*2);
md.referenceMesh.v2=zeros(neq,neq,length(hs)*2);
md.referenceMesh.dimat=zeros(neq,neq,length(hs)*2);
j=1;
for i=1:length(hs)
   md.referenceMesh.h(j)=hs(i)/2;
   [md.referenceMesh.v0(:,j),v,md.referenceMesh.dimat(:,:,j)] = computeVcoeffs(system.V,x,md.referenceMesh.h(j)); 
   md.referenceMesh.v1(:,:,j)=v(:,:,1); md.referenceMesh.v2(:,:,j)=v(:,:,2); 
   [md.referenceMesh.cu(:,:,:,j),md.referenceMesh.cup(:,:,:,j),md.referenceMesh.cv(:,:,:,j),md.referenceMesh.cvp(:,:,:,j)]=computeCcoeffs(md.referenceMesh.v0(:,j),v(:,:,1),v(:,:,2),hs(i)/2);
   x=x+md.referenceMesh.h(j);
   j=j+1;
   md.referenceMesh.h(j)=hs(i)/2;
   [md.referenceMesh.v0(:,j),v,md.referenceMesh.dimat(:,:,j)] = computeVcoeffs(system.V,x,md.referenceMesh.h(j)); 
   md.referenceMesh.v1(:,:,j)=v(:,:,1); md.referenceMesh.v2(:,:,j)=v(:,:,2); 
   [md.referenceMesh.cu(:,:,:,j),md.referenceMesh.cup(:,:,:,j),md.referenceMesh.cv(:,:,:,j),md.referenceMesh.cvp(:,:,:,j)]=computeCcoeffs(md.referenceMesh.v0(:,j),v(:,:,1),v(:,:,2),hs(i)/2);
   x=x+md.referenceMesh.h(j);
   j=j+1;
end


vs=md.v0;
[md.vmin,i]=min(vs(:));
md.matchingIndex=min(ceil(i/size(vs,1)),length(md.h));
md.referenceMesh.matchingIndex=md.matchingIndex*2;
end


function e_loc= estimateError(v0,v,h)
%returns error estimation

%eta-functions in Z=0:
eta0=1;
eta=[1/3,2/30,0.00952380952381,0.00105820105820]; %=get_eta(0,4)

v1=v(:,:,1);
v2=v(:,:,2);
v0 = diag(v0)*h^2;

v12=v1*v2-v2*v1;
v20 = v2*v0-v0*v2;
v10 = v1*v0-v0*v1;
v11=v1*v1;

err_u=v10*eta(1)/24+(-v11/24+4*v10/24-3*v20/24)*eta(2);
err_u=max(abs(err_u(:)));

err_up=(60*v2-v12)*eta0/120+(-60*v2-5*v20)*eta(1)/40+(+5*v12+120*v20)*eta(2)/40;
err_up=max(abs(err_up(:)))/h;

err_w=-(60*v2+v12)*eta(2)/120;
err_w=max(abs(err_w(:)))*h;

err_wp=(v10)*eta(1)/24+(-v11/24-v10/12+3*v20/24)*eta(2);
err_wp=max(abs(err_wp(:)));

e_loc = max([err_u err_up err_w err_wp]);
end

% 

function [cu,cup,cv,cvp]=computeCcoeffs(v0,v1,v2,h)
%computes the C^{U}_m, C^{W}_m matrices
v0=diag(v0)*h^2;
n=length(v0);
cu = zeros(n,n,3);
cup = zeros(n,n,3);
cv = zeros(n,n,3);
cvp =zeros(n,n,3);

v11=v1*v1;
v12=v1*v2-v2*v1;
v10=v1*v0-v0*v1;
v20=v2*v0-v0*v2;

%eta1:
cu(:,:,2)=-v1/2+1/24*v10;
%eta2:
cu(:,:,3)=-v11/24+4*v10/24-3*v20/24;

cup(:,:,1)=v2/2+v10/24-v12/120;
cup(:,:,2)=-3*v2/2-v11/24+v10/8-v20/8;
cup(:,:,3)=-7*v11/24+v12/8+3*v20;

cv(:,:,3)=-v2/2+v10/24;

cvp(:,:,2)=v1/2+v10/24;
cvp(:,:,3)=-v11/24-v10/12+3*v20/24;

end