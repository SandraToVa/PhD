%This is because some states of 1+- (p0) where wrongly identified.
%I am using the hole potentials and the cosntants found
%Programa auxiliar que uso per a calcular mes estats del espectre i saber
%en quina posició estan

%I need to define k1, k2 and r0
%Per al bottomonium: 
%setk1(-0.075)
%setk2(0.001)
setr0(3.964)
setL1(0.059)
setL3(-0.230)
load("dades.mat","m_b","m_c")
setm_q(m_c)
setspin(1)
% l= interpolació (0), llagures distncies (1), bad long distances (2), Vhf
% partial A (3), Vhf2 partial B (4)
setl(0)
%Espectre en la A=0 i B=0 - z(i) -
setk1(0)
setk2(0)

z=zeros(1,14);
[aux,~,~]=QuarkoniumS0J1(m_q,spin);
z(1)=aux(1);
z(5)=aux(2);
[aux,~,~]=Spin1Jcal0_1(m_q,spin);
z(2)=aux(1);
[aux,~,~]=Spin1Jcal1_2(m_q,spin);
z(3)=aux(1);
[aux,~,~]=Spin1Jcal2_1(m_q,spin);
z(4)=aux(1);
[aux,~,~]=Spin1Jcal0_2(m_q,spin);
z(6)=aux(1);
[aux,~,~]=Spin1Jcal2_2(m_q,spin);
z(8)=aux(1);
%[aux,W,x]=Spin1Jcal2_2(m_q,spin);
z(11)=aux(2);
[aux,~,~]=QuarkoniumS0J0(m_q,spin);
z(13)=aux(1);
[aux,~,~]=Spin1Jcal1_1(m_q,spin);
z(7)=aux(1);
%[aux,W,x]=Spin1Jcal1_1(m_q,spin);
z(10)=aux(2);
%[aux,W,x]=Spin1Jcal1_1(m_q,spin);
if m_q==1.4702
    z(14)=aux(3); %If charmonium
end
if m_q==4.8802
    z(14)=aux(5); %If bottomium
end
%[aux,W,x]=QuarkoniumS0J1(m_q,spin);
%c(12)=aux(2);
[aux,~,~]=QuarkoniumS0J2(m_q,spin);
z(9)=aux(1);
[aux,~,~]=Spin1Jcal3_1(m_q,spin);
z(12)=aux(1);

%%
[ES1J1p,~,~]=Spin1Jcal1_1(m_q,spin);

[ES0J1,~,~]=QuarkoniumS0J1(m_q,spin);
[ES1J0s,~,~]=Spin1Jcal0_1(m_q,spin);
[ES1J1s,~,~]=Spin1Jcal1_2(m_q,spin);
[ES1J2s,~,~]=Spin1Jcal2_1(m_q,spin);
[ES1J0p,~,~]=Spin1Jcal0_2(m_q,spin);

[ES1J2p,~,~]=Spin1Jcal2_2(m_q,spin);
[ES0J2,~,~]=QuarkoniumS0J2(m_q,spin);
[ES1J3p,~,~]=Spin1Jcal3_1(m_q,spin);
[ES0J0,~,~]=QuarkoniumS0J0(m_q,spin);
[ES1J3s,~,~]=Spin1Jcal3_1(m_q,spin);


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

%Usant la interpolació o només les llargues distàncies
%l=1; llargues distàncies
%l=0; interpolació
function x6=l
global v6
x6=v6;
end

function setl(val6)
global v6
v6=val6;
end