%Correspondence of the states of charmonium and bottomonium with
%the states computed with the matlab program

%Valor r0
setr0(3.964)
%Constants ajustades
setk1(0.11455)
setk2(0.00385)
%Valors Lambda g
setL1(0.07965)
setL3(0.3105)
% l= interpolació (0), llagures distncies (1), bad long distances (else)
setl(0)
%massa
load("dades.mat","m_c","m_b")
setm_q(m_c)
setspin(1) %s=1 for hybrids


%Vector en los valors de la energia que necesito
%HYPERFINE STATES
a=zeros(14,1);

%HYBRID NON HYPERFINE SPLITTING STATES
h=zeros(5,1);

%(s/d)_1%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   1-- (quarkonium)
[aux,~,~]=QuarkoniumS0J1(m_q,spin);
a(1)=aux(1);
h(1)=aux(1);
t(1)=4.0296;    %valor teòric
%   0-+ (hybrids)
[aux,~,~]=Spin1Jcal0_1(m_q,spin);
a(2)=aux(1);
t(2)=3.8976;    %valor teòric
%   1-+ (hybrids)
[aux,~,~]=Spin1Jcal1_2(m_q,spin);
a(3)=aux(1);
t(3)=3.9286;    %valor teòric
%   2-+ (hybrids)
[aux,~,~]=Spin1Jcal2_1(m_q,spin);
a(4)=aux(1);
t(4)=4.0746;    %valor teòric


%p_1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   1++ (quarkonium)
[aux,~,~]=QuarkoniumS0J1(m_q,spin);
a(12)=aux(2);
h(2)=aux(2);
t(12)=4.1436;   %valor teòric
%   0+- (hybrids)
[aux,~,~]=Spin1Jcal0_2(m_q,spin);
a(5)=aux(1);
t(5)=4.1106;    %valor teòric
%   1+- (hybrids)
[aux,~,~]=Spin1Jcal1_1(m_q,spin);
a(9)=aux(1);
t(9)=4.1116;    %valor teòric
%   2+- (hybrids)
[aux,~,~]=Spin1Jcal2_2(m_q,spin);
a(6)=aux(1);
t(6)=4.1756;    %valor teòric


%(p/f)_2%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   2++ (quarkonium)
[aux,~,~]=QuarkoniumS0J2(m_q,spin);
a(13)=aux(1);
h(3)=aux(1);
t(13)=4.2306;   %valor teòric
%   1+- (hybrids)
[aux,~,~]=Spin1Jcal1_1(m_q,spin);
a(10)=aux(2);
t(10)=4.1786;   %valor teòric
%   2+- (hybrids)
[aux,~,~]=Spin1Jcal2_2(m_q,spin);
a(7)=aux(2);
t(7)=4.2386;    %valor teòric
%   3+- (hybrids)
[aux,~,~]=Spin1Jcal3_1(m_q,spin);
a(14)=aux(1);
t(14)=4.2516;   %valor teòric


%p_0%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   0++ (quarkonium)
[aux,~,~]=QuarkoniumS0J0(m_q,spin);
a(8)=aux(1);
h(4)=aux(1);
t(8)=4.4396;    %valor teòric
%   1+- (hybrid)
[aux,~,~]=Spin1Jcal1_1(m_q,spin);
if m_q==1.4702
    a(11)=aux(3); %If charmonium
end
if m_q==4.8802
    a(11)=aux(5); %If bottomium
end
t(11)=4.5136;   %valor teòric charmonium

%d_2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[aux,~,~]=QuarkoniumS0J2(m_q,spin);
h(5)=aux(2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Now for each set of hybrid states, we want diferent energy levels (n)

%   p0 is the only coming from QuraknoiumS0J0
[aux,Y0,~]=QuarkoniumS0J0(m_q,spin);
p0=aux;

%   (s/d)1 shares with p1 the script QuraknoiumS0J1
[aux,Y1,~]=QuarkoniumS0J1(m_q,spin);
%It is not straigth forwad which states are p1 os (s/d)1 so we must chek
%We know that are (s/d)1:
positionsSD1= [1, 3, 5, 6, 8];
sd1=aux(positionsSD1);
%We know that are p1:
positionsP1= [2, 4, 7];
p1=aux(positionsP1);

%   (p/f)2 shares with d2 the script QuraknoiumS0J2
[aux,Y2,~]=QuarkoniumS0J2(m_q,spin);
%We know that are (s/d)1:
positionsPF2= [1, 3, 5, 6];
pf2=aux(positionsPF2);
%We know that are p1:
positionsD2= [2, 4, 7];
d2=aux(positionsD2);

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
