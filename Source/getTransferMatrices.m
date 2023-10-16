function varargout = getTransferMatrices(varargin)
% [T,Te]=GETTRANSFERMATRICES(E,meshData,order,diag_form)
%computes and returns the transfer matrices of the CA method


E=varargin{1}; md = varargin{2};
diag_form = false; order = 2;
if nargin > 2
    order = varargin{3};  
    %order determines if the transfer matrix of the second order or the
    %sixth order method needs to be returned
    if nargin >3
        %diag_form determines if U^D or U needs to be returned
        diag_form= varargin{4};
    end
end
hm=md.h; v0m=md.v0;
len=size(v0m,1);

%preallocation:
T=zeros(len*2,len*2,length(hm));
if nargout==2
    Te=zeros(len*2,len*2,length(hm));
end

p = 2;
if order==6
    cum=md.cu; cupm=md.cup; cvm=md.cv; cvpm=md.cvp;
    p = 3;
end
if p==3
    wght1 = [ 0.001058201058201  0.000048100048100  0.000000925000925  0.000000010277788 ...
      0.000000000075572   0.000000000000398   0.000000000000002   0.000000000000000];
    wght2 = [ 0.009523809523810   0.000529100529101   0.000012025012025   0.000000154166821 ...
      0.000000001284724   0.000000000007557   0.000000000000033   0.000000000000000];
else
    wght1 = [0.066666666666667;   0.004761904761905;   0.000132275132275 ;  0.000002004168671; ...
        0.000000019270853;   0.000000000128472;   0.000000000000630;   0.000000000000002];
    wght2 = [0.333333333333333; 0.033333333333333;  0.001190476190476;  0.000022045855379; ...
             0.000000250521084; 0.000000001927085;  0.000000000010706; 0.000000000000045];     
end

for k=1:length(hm)
    v0=v0m(:,k);
    h=hm(k);
    if diag_form
       D=eye(len); 
    else
       D=md.dimat(:,:,k); 
    end
    DT=D';
    
    %prepare relevant one step quantities z's and eta's
    Zm = (v0-E)*h^2;   
%    [xsi,eta0,eta] = computeEtaMatrices(Zm,p);
    lZm=length(Zm); 
    for s=1:lZm
        Z=Zm(s);
        za = abs(Z);
        if za > 0.6    
            squ=sqrt(Z)*j;
            eta0(s)=sin(squ)/squ;
            xsi(s)=cos(squ);
            eta(s,1)=(xsi(s)-eta0(s))/Z;
            eta(s,2)=(eta0(s)-3*eta(s,1))/Z;
            eta(s,3)=(eta(s,1)-5*eta(s,2))/Z;
        elseif p ==3
            hlp=(Z.^(0:length(wght1)-1))';   
            eta(s,4) = wght1*hlp;
            eta(s,3) = wght2*hlp;
            for i=2:-1:1
                eta(s,i) = Z * eta(s,i+2) + (2*i+3) * eta(s,i+1);
            end
            eta0(s) = Z * eta(s,2) + 3*eta(s,1);  
            xsi(s) = Z * eta(s,1) + eta0(s);
       else %if p ==2
            hlp=(Z.^(0:length(wght1)-1));   
            eta(s,2) = hlp*wght1;
            eta(s,1) = hlp*wght2;
            eta0(s) = Z * eta(s,2) + 3*eta(s,1);  
            xsi(s) = Z * eta(s,1) + eta0(s);      
        end
    end
   
    %entries of the second order transfer matrix:
    u= diag(xsi);           %U
    up = diag(Zm'.*eta0);   %U'
    v = diag(eta0);         %W
    vp= diag(xsi);          %W'
    
    if order == 6
      %extra correction terms are computed  
      cu=cum(:,:,:,k); cup=cupm(:,:,:,k);
      cv=cvm(:,:,:,k); cvp=cvpm(:,:,:,k);
      tmpu=cu(:,:,2)*diag(eta(:,1))+cu(:,:,3)*diag(eta(:,2));
      tmpup=cup(:,:,1)*diag(eta0(:))+cup(:,:,2)*diag(eta(:,1))+cup(:,:,3)*diag(eta(:,2));
      tmpv = cv(:,:,3)*diag(eta(:,2));  
      tmpvp=cvp(:,:,2)*diag(eta(:,1))+cvp(:,:,3)*diag(eta(:,2)); 
      
      
      u  = u + tmpu;
      up = up + tmpup;
      v = v + tmpv;
      vp = vp + tmpvp;
    end

    u  = D*u*DT;
    up = D*(up)/h*DT;
    v = D*(v)*h*DT;  
    vp = D*vp*DT;        

    T(:,:,k)=[u, v ; up, vp];
   
    if nargout ==2
       %derivatives with respect to E need to be computed 
       ue = diag(eta0);  
       uep = diag(xsi + eta0);
       ve = diag(eta(:,1));
       vep = diag(eta0);
       if order == 6
         tmpue=cu(:,:,1)*diag(eta(:,1))+cu(:,:,2)*diag(eta(:,2))+cu(:,:,3)*diag(eta(:,3)); 
         tmpupe=cup(:,:,1)*diag(eta(:,1))+cup(:,:,2)*diag(eta(:,2))+cup(:,:,3)*diag(eta(:,3));
         tmpve = cv(:,:,2)*diag(eta(:,2))+cv(:,:,3)*diag(eta(:,3));  
         tmpvpe=cvp(:,:,1)*diag(eta(:,1))+cvp(:,:,2)*diag(eta(:,2))+cvp(:,:,3)*diag(eta(:,3)); 
         
         ue = ue + tmpue;  
         uep = uep + tmpupe;   
         ve = ve + tmpve;
         vep = vep + tmpvpe;
       end

       ue = D*(-ue*h^2/2)*DT; 
       uep = D*(-uep*h/2)*DT; 
       ve = D*(-ve*(h^3)/2)*DT; 
       vep = D*(-vep*h^2/2)*DT; 
       Te(:,:,k)=[ue, ve ; uep, vep];
    end
end
varargout{1}=T;
if nargout==2
    varargout{2}=Te;
end
end

% function [xsi,eta0,eta] = computeEtaMatrices(Zm,m)
% if m==3
%     wght1 = [ 0.001058201058201  0.000048100048100  0.000000925000925  0.000000010277788 ...
%       0.000000000075572   0.000000000000398   0.000000000000002   0.000000000000000];
%     wght2 = [ 0.009523809523810   0.000529100529101   0.000012025012025   0.000000154166821 ...
%       0.000000001284724   0.000000000007557   0.000000000000033   0.000000000000000];
% else
%     wght1 = [0.066666666666667;   0.004761904761905;   0.000132275132275 ;  0.000002004168671; ...
%         0.000000019270853;   0.000000000128472;   0.000000000000630;   0.000000000000002];
%     wght2 = [0.333333333333333; 0.033333333333333;  0.001190476190476;  0.000022045855379; ...
%              0.000000250521084; 0.000000001927085;  0.000000000010706; 0.000000000000045];     
% end
% for k=1:length(Zm)
% Z=Zm(k);
% za = abs(Z);
% 
% if za > 1    
%     squ=sqrt(Z)*j;
%     eta0(k)=sin(squ)/squ;
%     xsi(k)=cos(squ);
%     eta(k,1)=(xsi(k)-eta0(k))/Z;
%     eta(k,2)=(eta0(k)-3*eta(k,1))/Z;
%     eta(k,3)=(eta(k,1)-5*eta(k,2))/Z;
% elseif m ==3
%      hlp=(Z.^(0:length(wght1)-1))';   
%     eta(k,4) = wght1*hlp;
%     eta(k,3) = wght2*hlp;
%     for i=2:-1:1
%         eta(k,i) = Z * eta(k,i+2) + (2*i+3) * eta(k,i+1);
%     end
%     eta0(k) = Z * eta(k,2) + 3*eta(k,1);  
%     xsi(k) = Z * eta(k,1) + eta0(k);
% else %if m ==2
%     hlp=(Z.^(0:length(wght1)-1));   
%     eta(k,2) = hlp*wght1;
%     eta(k,1) = hlp*wght2;
%     eta0(k) = Z * eta(k,2) + 3*eta(k,1);  
%     xsi(k) = Z * eta(k,1) + eta0(k);      
% end
% end
% end
