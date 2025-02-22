%Valor r0
setr0(3.964)
%massa
load("dades.mat","m_c","m_b")
setm_q(m_c)


%values (autentic)
Min=0.27;
Mfin=0.57;

%Now we compute I_theta for 4s->3s for diferent M
%|Ef-Ei|>M

%number of spaces:
k=50;
%M row:
nM=(Mfin-Min)/k;
M=Min:nM:Mfin;
Itheta=zeros(1,k+1);

Itheta0=zeros(1,k+1);
Itheta1=zeros(1,k+1);
Itheta2=zeros(1,k+1);

n=1;
while n<=k+1
    %For transitions with only one option like DtoS or StoS
    %Itheta(n)=ComputeExpValM(I,F,M(n),@DtoStrans);

    %For transitions with diferent m, DtoD we must do the sum and average
    %over all states (D0toD0 + etc)
    %example for l=2->l=2
    %Itheta(n)=ComputeExpValM(I,F,M(n),@D0toD0trans);
    %Itheta(n)=Itheta(n)+2*ComputeExpValM(I,F,M(n),@D1toD1trans);
    %Itheta(n)=Itheta(n)+2*ComputeExpValM(I,F,M(n),@D2toD2trans);
    %Itheta(n)=Itheta(n)/5; 
    %Here it is 5 because 5=(2*l+1),l=2
    %2* because we have m=+-1 -> m=+-1 and m=+-2 -> m=+-2

    %Itheta0(n)=ComputeExpValM(I,F,M(n),@D0toD0trans)^2;
    %Itheta1(n)=ComputeExpValM(I,F,M(n),@D1toD1trans)^2;
    %Itheta2(n)=ComputeExpValM(I,F,M(n),@D2toD2trans)^2;
   
    transition=I_thetaFunctions('HQp0toD2');
    ExpValFunc=ExpValFunctions('HQ');

    Itheta0(n)=ExpValFunc(1,1,M(n),transition,0,2,false);
    Itheta(n)=Itheta0(n)^2;

    %transition0=I_thetaFunctions('QQD0toD0');
    %transition1=I_thetaFunctions('QQD1toD1');
    %transition2=I_thetaFunctions('QQD2toD2');

    %Itheta0(n)=ExpValFunc(3,1,M(n),transition0,2,2);
    %Itheta1(n)=ExpValFunc(3,1,M(n),transition1,2,2);
    %Itheta2(n)=ExpValFunc(3,1,M(n),transition2,2,2);
    %Itheta(n)=(Itheta0(n)^2+2*Itheta1(n)^2+2*Itheta2(n)^2)/5;

    n=n+1;

end

%graphic
scatter(M,Itheta)
Mtrans=M.';
Itrans=Itheta.';

 



%VALOR DE LA MASSA
%m=1.4702; charm
%m=4.8802; bottom
function x0=m_q
global v0
x0=v0;
end

function setm_q(val0)
global v0
v0=val0;
end

% VALOR DEL SPIN
function x7=spin
global v7
x7=v7;
end 

function setspin(val7)
global v7
v7 = val7;
end

% VALOR DEL r0
function x3=r0
global v3
x3=v3;
end 

function setr0(val3)
global v3
v3 = val3;
end

