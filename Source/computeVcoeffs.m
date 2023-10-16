function varargout=computeVcoeffs(V,x,h)
%[v0,v,D]=computeVcoeffs(V,x,h)
%   computes the V-coefficients needed by the CP method for the mesh interval [x, x+h]
%      INPUT: V: function handle to potential
%             x, h: determine mesh interval
%      OUTPUT: v0: the diagonal elements of the diagonal matrix V_0^D 
%              v:  v(:,:,1)= V_1^D, v(:,:,2)= V_2^D
%
%[v0,D]=computeVcoeffs(V,x,h)    % can be used when only second order
%                                % algorithm will be applied and V_1 and 
%                                % V_2 are not needed
if nargout == 2
    V0 = V(x+h/2);  %for the second order algorithm, replacing V by the midpoint value suffices
    [eigenvectors,eigenvalues]=eig(V0);
    varargout{1} = diag(eigenvalues);
    varargout{2}=eigenvectors;
else
    xx=[0.339981043584856,0.861136311594053]; %Gauss quadrature points
    w=[0.652145154862546,0.347854845137454]; %Gauss weights
    a=x; b=x+h;
    m = (b+a)/2; d = h/2;
    l1= legendre((1-xx)/2); 
    l2= legendre((1+xx)/2);
    dx = d*xx;  
    ffm=V([m+dx m-dx]); n = size(ffm,1);
    ress=zeros(n,n,3);
    for i=1:n
       for j=1:i
            fp = ffm(i,j,1:2);
            fm = ffm(i,j,3:end);
            res = ((w.*fm(:)')*l1+(w.*fp(:)')*l2);
            ress(i,j,:)=res;
            ress(j,i,:)=res;
        end
    end
    qpot=ress*h/2;
    V0 = qpot(:,:,1)/h;
    [eigenvectors,eigenvalues]=eig(V0); %v0*eigenvectors=eigenvectors*eigenvalues
    v0 = diag(eigenvalues);
    %Next compute v1,v2 in the workbasis
    D = eigenvectors; 
    varargout{1}=v0;
    varargout{3}=D;
    v(:,:,1)=D'*(h*qpot(:,:,2)*3)*D;
    v(:,:,2)=D'*(h*qpot(:,:,3)*5)*D;
    varargout{2}=v;
end
end

function res = legendre(gammas)
res(:,1)=ones(1,length(gammas));
res(:,2) = -1 + 2*gammas;
x = 2*gammas-1;
res(:,3) = (3*x'.*res(:,2)-res(:,1))/2;
end