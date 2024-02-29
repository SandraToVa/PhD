%Valor r0
setr0(3.964)
setk1(0.11455)
setk2(0.00385)
setL1(0.07965)
setL3(0.3105)
load("dades.mat","m_c","m_b")
setm_q(m_c)
setspin(1)

%Totes les E en GeV
%Vector de valors de E calculades
a=zeros(10);
%Vector de valors de E teòriques
t(1)=4.1106;
t(2)=4.1756;
t(3)=4.2386;
t(4)=4.4396;
t(5)=4.1116;
t(6)=4.1786;
t(7)=4.5136;
t(8)=4.1436;
t(9)=4.2306;
t(10)=4.2516;
%Vector error de energies teòriques
e(1)=0.0276;
e(2)=0.0186;
e(3)=0.0266;
e(4)=0.0466;
e(5)=0.0236;
e(6)=0.0276;
e(7)=0.0536;
e(8)=0.0256;
e(9)=0.0326;
e(10)=0.0346;

%Vector en los valors de la energia que necesito
[aux,W,x]=Spin1Jcal0_2(m_q,spin);
a(1)=aux(1);
[aux,W,x]=Spin1Jcal2_2(m_q,spin);
a(2)=aux(1);
[aux,W,x]=Spin1Jcal2_2(m_q,spin);
a(3)=aux(2);
[aux,W,x]=QuarkoniumS0J0(m_q,spin);
a(4)=aux(1);
[aux,W,x]=Spin1Jcal1_1(m_q,spin);
a(5)=aux(1);
[aux,W,x]=Spin1Jcal1_1(m_q,spin);
a(6)=aux(2);
[aux,W,x]=Spin1Jcal1_1(m_q,spin);
if m_q==1.4702
    a(11)=aux(3); %If charmonium
end
if m_q==4.8802
    a(11)=aux(5); %If bottomium
end
[aux,W,x]=QuarkoniumS0J1(m_q,spin);
a(8)=aux(2);
[aux,W,x]=QuarkoniumS0J2(m_q,spin);
a(9)=aux(1);
[aux,W,x]=Spin1Jcal3_1(m_q,spin);
a(10)=aux(1);
    
chi=0;
for i=1:10
    chi=chi+((a(i)-t(i))^2)/((e(i))^2);
end


chi
a

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

% VALOR DEL SPIN
function x4=spin
global v4
x4=v4;
end 

function setspin(val4)
global v4
v4 = val4;
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




