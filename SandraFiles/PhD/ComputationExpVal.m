%Valor r0
setr0(3.964)
%massa
load("dades.mat","m_c","m_b")
setm_q(m_c)
 

%hybrid
setspin(1);
sin=spin;
[Einicial,Yi,xi]=QuarkoniumS0J0(m_q,spin);

%quarkonium
setspin(0)
sfin=spin;
[Efinal,Yf,xf]=QuarkoniumS0J0(m_q,spin);

%Operator that we want to find the expectated value
%WHICH x??? xf or xi because they are different!!!
if length(xi)~=length(xf)
    disp('Change tolerance of hybrid:');
    length(xi)-length(xf)
else
    x=xi; %The same for xi and xf same x
    Op=diag(x);

    %PROBLEM WITH THE FUNCTION, x for hybrid and x for quarkonium have diferent
    %dimension. Due to ConstrucMesh?????
    V=ExpectedValue(1,2,Op,@QuarkoniumS0J0,@QuarkoniumS0J0,m_q,sin,sfin);
end

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

