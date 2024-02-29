%This is because some states of 1+- (p0) where wrongly identified.
%I am using the hole potentials and the cosntants found
%Programa auxiliar que uso per a calcular mes estats del espectre i saber
%en quina posició estan

%I need to define k1, k2 and r0
%Per al bottomonium: 
setk1(0.03411)
setk2(-0.00085)
setr0(3.964)
setL1(0.07965)
setL3(0.3105)
load("dades.mat","m_b","m_c")
setm_q(m_c)
setspin(1)

[E,W,x]=Spin1Jcal0_2(m_q,spin);

%I am identifing 1+- estates with the 1,2,3 estates found:

%p1=aux(1)
%(p/f)2=aux(2)
%p0=aux(3)   El que possiblement està mal identificat!

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

% VALOR DEL Vhf
function x1=k1
global v1
x1=v1;
end 

function setk1(val1)
global v1
v1 = val1;
end

% VALOR DEL Vhf2
function x2=k2
global v2
x2=v2;
end 

function setk2(val2)
global v2
v2 = val2;
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

% VALOR DEL -g\Lambda' (el signe negatiu ja esta en la Vsp)
function x4=L1
global v4
x4=v4;
end 

function setL1(val4)
global v4
v4 = val4;
end

% VALOR DEL -g\Lambda''' (per a la Vpp)
function x5=L3
global v5
x5=v5;
end 

function setL3(val5)
global v5
v5 = val5;
end