function [eigvData,meshData]=computeEigenvalues(system,kmin,kmax,tol)
%[EigvData,meshData]=computeEigenvalues(system,kmin,kmax,tol)
% computes the eigenvalues of a coupled-channel Schrodinger system
% INPUT:
%       - system: a structure containing information on the Schrodinger
%                 problem
%                 system.a, system.b: endpoints of the integration interval
%                 system.A1,system.A2:
%                 system.B1,system.B2: parameters of the boundary conditions
%                 system.V: function handle to the potential matrix
%                           function
%       - kmin , kmax: lowest and highest index of the range of eigenvalues
%                      to be computed
%       - tol: user input tolerance (~ relative error in the eigenvalue
%              approximations)
% OUTPUT:
%       - eigvData: a Matlab structure 
%                   eigvData.eigenvalues = vector containing the eigenvalue
%                                          approximations
%                   eigvData.status = vector of same length as eigvData.eigenvalues
%                            eigvData.status(i)= 1 : no difficulties detected in the approximation of the i-th eigenvalue  
%                            eigvData.status(i)= -1: the newton-raphson process didn't converge to an eigenvalue, try a lower input tolerance   
%                            eigvData.status(i)= -2: too many iterations in the newton-raphson process. The requested
%                                                    accuracy may not have been reached.
%                   eigvData.indices = eigenvalue indices
%       
%       - meshData: structure collecting all data associated to the mesh
%                   constructed
%                   meshData.tol : expected accuracy for the results
%                                  computed over this mesh
%                   meshData.v0: vector of V_0 matrices, length of the vector
%                                equals the number of mesh intervals
%                   meshData.v1: vector of {\bar V}_1 matrices, length of the vector
%                                equals the number of mesh intervals
%                   meshData.v2: vector of {\bar V}_2 matrices, length of the vector
%                                equals the number of mesh intervals
%                   meshData.h: vector containing the length of the mesh intervals
%                   meshData.dimat: vector of D matrices (orthogonal matrix
%                                  used by the diagonalization process), the length of the vector
%                                  equals the number of mesh intervals
%                   meshData.cu,cup,cv,cvp: C-coefficients used by the
%                                           sixth order CP algorithm

%
addpath('./source')


%check validity of the input parameters
if kmin > kmax || kmin < 0
   error(SchrodSys:argChk, 'No eigenvalue is found in the given interval.\n Check the following conditions: kmin < kmax and kmin>0');   
end
if tol<eps
   error(SchrodSys:argChk, 'The user input tolerance cannot be smaller than the floating point accuracy');    
end

%construct a mesh corresponding to the user input tolerance tol
meshData = constructMesh(system,tol);

%compute eigenvalues with index between kmin and kmax
initE=meshData.vmin;
step = (mean(meshData.v0(:))-initE)/100;
k=kmin;
while k<=kmax
  %getBracket uses the indexing function to determine good starting values
  %for the shooting process
  [Elow,Eup,status1] = getBracket(system,meshData.referenceMesh,k,initE,step);
  %shoot for the eigenvalue with index k
  [E,Eref,multiplicity,status2] = shootForEigenvalue(system,meshData,Elow,Eup,k);
  eigvData.eigenvalues(k-kmin+1:k-kmin+multiplicity)=E;
  eigvData.errorEstimations(k-kmin+1:k-kmin+multiplicity)=abs(Eref-E);
  eigvData.status(k-kmin+1:k-kmin+multiplicity)=min(status1,status2);
  eigvData.indices(k-kmin+1:k-kmin+multiplicity)=k:k+multiplicity-1;
  step=max((E-initE)/10,step/10); initE=E-max(tol*10,1e-9); 
  k=k+multiplicity;
end
end



