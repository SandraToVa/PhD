
%Same as IntervalsConfian but for the extra states for charm and bottom
%Instead of computing the spectrum with the avergae values of the intervals
%this way we compute it with the linearization

%Valor r0
setr0(3.964)
setL1(0.059)
setL3(-0.230)
load("dades.mat","m_c","m_b")
setm_q(m_c)
setspin(1)
% l= interpolació (0), llagures distncies (1), bad long distances (else)
setl(0)

%Primer calculem l'espectre en la A i B òptimes - f(i) -
setk1(-0.0848)
setk2(0.0021)

[aux,~,~]=QuarkoniumS0J1(m_q,spin);
f(1)=aux(3);
if m_q==4.8802
    f(9)=aux(4);
end
[aux,~,~]=Spin1Jcal0_1(m_q,spin);
f(2)=aux(2);
[aux,~,~]=Spin1Jcal1_2(m_q,spin);
f(3)=aux(3);
f(6)=aux(2);
[aux,~,~]=Spin1Jcal2_1(m_q,spin);
f(4)=aux(3);
f(7)=aux(2);
[aux,~,~]=QuarkoniumS0J2(m_q,spin);
f(5)=aux(2);
if m_q==4.8802
    f(13)=aux(3);
end
%[aux,~,~]=Spin1Jcal1_2(m_q,spin);
%f(6)=aux(2);
%[aux,~,~]=Spin1Jcal2_1(m_q,spin);
%f(7)=aux(2);
[aux,~,~]=Spin1Jcal3_2(m_q,spin);
f(8)=aux(1);

if m_q==4.8802
    %[aux,~,~]=QuarkoniumS0J1(m_q,spin);
    %f(1)=aux(4);
    [aux,~,~]=Spin1Jcal0_2(m_q,spin);
    f(10)=aux(2);
    [aux,~,~]=Spin1Jcal1_1(m_q,spin);
    f(11)=aux(3);
    f(14)=aux(4);
    [aux,~,~]=Spin1Jcal2_2(m_q,spin);
    f(12)=aux(3);
    f(15)=aux(5);
    %[aux,~,~]=QuarkoniumS0J2(m_q,spin);
    %f(13)=aux(3);
    %[aux,~,~]=Spin1Jcal1_1(m_q,spin);
    %f(14)=aux(4);
    %[aux,~,~]=Spin1Jcal2_2(m_q,spin);
    %f(15)=aux(5);
    [aux,~,~]=Spin1Jcal3_1(m_q,spin);
    f(16)=aux(3);   
end


fprintf('f,');
 
%Espectre en la A=0 i B=0 - z(i) -
setk1(0)
setk2(0)

[aux,~,~]=QuarkoniumS0J1(m_q,spin);
c(1)=aux(3);
if m_q==4.8802
    c(9)=aux(4);
end
[aux,~,~]=Spin1Jcal0_1(m_q,spin);
c(2)=aux(2);
[aux,~,~]=Spin1Jcal1_2(m_q,spin);
c(3)=aux(3);
c(6)=aux(2);
[aux,~,~]=Spin1Jcal2_1(m_q,spin);
c(4)=aux(3);
c(7)=aux(2);
[aux,~,~]=QuarkoniumS0J2(m_q,spin);
c(5)=aux(2);
if m_q==4.8802
    c(13)=aux(3);
end
%[aux,~,~]=Spin1Jcal1_2(m_q,spin);
%c(6)=aux(2);
%[aux,~,~]=Spin1Jcal2_1(m_q,spin);
%c(7)=aux(2);
[aux,~,~]=Spin1Jcal3_2(m_q,spin);
c(8)=aux(1);

if m_q==4.8802
    %[aux,~,~]=QuarkoniumS0J1(m_q,spin);
    %c(1)=aux(4);
    [aux,~,~]=Spin1Jcal0_2(m_q,spin);
    c(10)=aux(2);
    [aux,~,~]=Spin1Jcal1_1(m_q,spin);
    c(11)=aux(3);
    c(14)=aux(4);
    [aux,~,~]=Spin1Jcal2_2(m_q,spin);
    c(12)=aux(3);
    c(15)=aux(5);
    %[aux,~,~]=QuarkoniumS0J2(m_q,spin);
    %c(13)=aux(3);
    %[aux,~,~]=Spin1Jcal1_1(m_q,spin);
    %c(14)=aux(4);
    %[aux,~,~]=Spin1Jcal2_2(m_q,spin);
    %c(15)=aux(5);
    [aux,~,~]=Spin1Jcal3_1(m_q,spin);
    c(16)=aux(3);   
end
    
fprintf('c,');

%Espectre A+DeltaA,B - fa(i) - 

setk1(-0.08479)
setk2(0.0021)

[aux,~,~]=QuarkoniumS0J1(m_q,spin);
fa(1)=aux(3);
if m_q==4.8802
    fa(9)=aux(4);
end
[aux,~,~]=Spin1Jcal0_1(m_q,spin);
fa(2)=aux(2);
[aux,~,~]=Spin1Jcal1_2(m_q,spin);
fa(3)=aux(3);
fa(6)=aux(2);
[aux,~,~]=Spin1Jcal2_1(m_q,spin);
fa(4)=aux(3);
fa(7)=aux(2);
[aux,~,~]=QuarkoniumS0J2(m_q,spin);
fa(5)=aux(2);
if m_q==4.8802
    fa(13)=aux(3);
end
%[aux,~,~]=Spin1Jcal1_2(m_q,spin);
%fa(6)=aux(2);
%[aux,~,~]=Spin1Jcal2_1(m_q,spin);
%fa(7)=aux(2);
[aux,~,~]=Spin1Jcal3_2(m_q,spin);
fa(8)=aux(1);

if m_q==4.8802
    %[aux,~,~]=QuarkoniumS0J1(m_q,spin);
    %fa(1)=aux(4);
    [aux,~,~]=Spin1Jcal0_2(m_q,spin);
    fa(10)=aux(2);
    [aux,~,~]=Spin1Jcal1_1(m_q,spin);
    fa(11)=aux(3);
    fa(14)=aux(4);
    [aux,~,~]=Spin1Jcal2_2(m_q,spin);
    fa(12)=aux(3);
    fa(15)=aux(5);
    %[aux,~,~]=QuarkoniumS0J2(m_q,spin);
    %fa(13)=aux(3);
    %[aux,~,~]=Spin1Jcal1_1(m_q,spin);
    %fa(14)=aux(4);
    %[aux,~,~]=Spin1Jcal2_2(m_q,spin);
    %fa(15)=aux(5);
    [aux,~,~]=Spin1Jcal3_1(m_q,spin);
    fa(16)=aux(3);   
end

fprintf('fa,');

%Espectre A,B+DeltaB - fb(i) -

setk1(-0.0848)
setk2(0.00211)

[aux,~,~]=QuarkoniumS0J1(m_q,spin);
fb(1)=aux(3);
if m_q==4.8802
    fb(9)=aux(4);
end
[aux,~,~]=Spin1Jcal0_1(m_q,spin);
fb(2)=aux(2);
[aux,~,~]=Spin1Jcal1_2(m_q,spin);
fb(3)=aux(3);
fb(6)=aux(2);
[aux,~,~]=Spin1Jcal2_1(m_q,spin);
fb(4)=aux(3);
fb(7)=aux(2);
[aux,~,~]=QuarkoniumS0J2(m_q,spin);
fb(5)=aux(2);
if m_q==4.8802
    fb(13)=aux(3);
end
%[aux,~,~]=Spin1Jcal1_2(m_q,spin);
%fb(6)=aux(2);
%[aux,~,~]=Spin1Jcal2_1(m_q,spin);
%fb(7)=aux(2);
[aux,~,~]=Spin1Jcal3_2(m_q,spin);
fb(8)=aux(1);

if m_q==4.8802
    %[aux,~,~]=QuarkoniumS0J1(m_q,spin);
    %fb(1)=aux(4);
    [aux,~,~]=Spin1Jcal0_2(m_q,spin);
    fb(10)=aux(2);
    [aux,~,~]=Spin1Jcal1_1(m_q,spin);
    fb(11)=aux(3);
    fb(14)=aux(4);
    [aux,~,~]=Spin1Jcal2_2(m_q,spin);
    fb(12)=aux(3);
    fb(15)=aux(5);
    %[aux,~,~]=QuarkoniumS0J2(m_q,spin);
    %fb(13)=aux(3);
    %[aux,~,~]=Spin1Jcal1_1(m_q,spin);
    %fb(14)=aux(4);
    %[aux,~,~]=Spin1Jcal2_2(m_q,spin);
    %fb(15)=aux(5);
    [aux,~,~]=Spin1Jcal3_1(m_q,spin);
    fb(16)=aux(3);   
end
    
fprintf('fb,');  

% Vector de a_i i b_i - a(i) i b(i) - 
if m_q==1.4702
    a=zeros(1,8);
    b=zeros(1,8);
    for i=1:8
        a(i)=(fa(i)-f(i))/0.00001;
        b(i)=(fb(i)-f(i))/0.00001;
    end
end
if m_q==4.8802
    a=zeros(1,16);
    b=zeros(1,16);
    for i=1:16
        a(i)=(fa(i)-f(i))/0.00001;
        b(i)=(fb(i)-f(i))/0.00001;
    end
end

fprintf('a,b'); 

%%-----------------------------------------------------------------------

%Les dades de l'interval les hem trobat executant el codi per als valors
%normal de l'espectre. Usem els valors de ci que obtenim alli i els posem
%aqui com a dades

%Intervals de confiança en level=0.68 ja que és el corresponent a 1 sigma.
%p00, p10 i p01
%Son els extrem dels valors de l'interval ja trobats
ci(1,2)=-0.093330978878651;
ci(2,2)=-0.0562619988852521;
ci(1,3)=-0.00748000664716397;
ci(2,3)=0.0100221041916237;


%k1=p10 i k2=p01
%Calcular primer en k1 i k2 i despues en los intervals resultats trobar los
%millors valors per a p10 i p01
%columna del mig interval A p10
%tercer columna intervla B p01
p10=-0.075;
p01=0.001;


%Ara caclulem l'espectre amb E_i=c_i+b_iB+a_iA
%Tambe calculem lo interval de confiança

%Interval de confiança de p10 i p01 són ci
DeltaA=ci(2,2)-p10;
deltaA=p10-ci(1,2);
DeltaB=ci(2,3)-p01;
deltaB=p01-ci(1,3);

fprintf(',deltes A i B');

if m_q==1.4702
    e=zeros(1,8);
    Deltae=zeros(1,8);
    deltae=zeros(1,8);
    for i=1:8
        e(i)=c(i)+ b(i)*p01 + a(i)*p10;
        Deltae(i)=abs(a(i))*DeltaA + abs(b(i))*DeltaB;
        deltae(i)=abs(a(i))*deltaA + abs(b(i))*deltaB;
    end
end
if m_q==4.8802
    e=zeros(1,16);
    Deltae=zeros(1,16);
    deltae=zeros(1,16);
    for i=1:16
        e(i)=c(i)+ b(i)*p01 + a(i)*p10;
        Deltae(i)=abs(a(i))*DeltaA + abs(b(i))*DeltaB;
        deltae(i)=abs(a(i))*deltaA + abs(b(i))*deltaB;
    end
end

E=transpose(e);
DeltaE=transpose(Deltae);
deltaE=transpose(deltae);

fprintf(',E,deltes E');



%-----------------------------------------------------------------------

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