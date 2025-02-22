%Valor r0
setr0(3.964)
setk1(0.01926)
setk2(0.00453)
setL1(0.07965)
setL3(0.3105)
load("dades.mat","m_c","m_b")
setm_q(m_b)

%Programa que calcula los valors del espectre per als nivells del c
% 2(s/d)2 and 1d2
c=zeros(8);

%Vector en los valors de la energia que necesito
aux=QuarkoniumS0J1(m_q);
c(1)=aux(3);
aux=Spin1Jcal0_1(m_q);
c(2)=aux(2);
aux=Spin1Jcal1_2(m_q);
c(3)=aux(3);
aux=Spin1Jcal2_1(m_q);
c(4)=aux(3);
aux=QuarkoniumS0J2(m_q);
c(5)=aux(2);
aux=Spin1Jcal1_2(m_q);
c(6)=aux(2);
aux=Spin1Jcal2_1(m_q);
c(7)=aux(2);
aux=Spin1Jcal3_2(m_q);
c(8)=aux(1);

c

%Per al bottom son los estats c + los estats b

if m_q==4.8802
    b=zeros(8);
    aux=QuarkoniumS0J1(m_q);
    b(1)=aux(4);
    aux=Spin1Jcal0_2(m_q);
    b(2)=aux(2);
    aux=Spin1Jcal1_1(m_q);
    b(3)=aux(3);
    aux=Spin1Jcal2_2(m_q);
    b(4)=aux(3);
    aux=QuarkoniumS0J2(m_q);
    b(5)=aux(3);
    aux=Spin1Jcal1_1(m_q);
    b(6)=aux(4);
    aux=Spin1Jcal2_2(m_q);
    b(7)=aux(5);
    aux=Spin1Jcal3_1(m_q);
    b(8)=aux(3);

    b
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
