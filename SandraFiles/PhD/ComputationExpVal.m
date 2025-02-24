%Valor r0
setr0(3.964)
%massa
load("dades.mat","m_c","m_b")
setm_q(m_b)

%values (autentic)
Min=0.27;
Mfin=0.62; %Has to be changed by hand depnending on the trasnition

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
IthetaN=zeros(1,k+1); %N is negative =-1
IthetaM=zeros(1,k+1); %M is dobule negative =-2

n=1;
while n<=k+1
%Inser here the transition following the guide of TransitionsAdded.m
%
transition0=I_thetaFunctions('HQsdtoP0');
transition1=I_thetaFunctions('HQsdtoP1');
transitionN=I_thetaFunctions('HQsdtoPn'); %n of negative

ExpValFunc=ExpValFunctions('HQ');

Itheta0(n)=ExpValFunc(3,2,M(n),transition0,1,1,true);
Itheta1(n)=ExpValFunc(3,2,M(n),transition1,1,1,true);
IthetaN(n)=ExpValFunc(3,2,M(n),transitionN,1,1,true);

%The N is not n in n(s/d)_J: remember the ordering of states for (s/d)1 and p1

Itheta(n)=( 2*(Itheta1(n)-IthetaN(n))*conj(Itheta1(n)-IthetaN(n)) )/3;


%%%%%
fprintf('Number of iteration: %d / %d\n', n, k);
    n=n+1;

end

% Crear el scatter plot
figure; % Crear una nueva figura
scatter(M, Itheta, 60, 'filled'); % Scatter plot con puntos llenos y tamaño 100

% Cambiar los nombres de los ejes
xlabel('Dipion mass (GeV)');
ylabel('I_{\theta}^2 (1/GeV)');

Mtrans=M';
Itrans=Itheta';

 



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

