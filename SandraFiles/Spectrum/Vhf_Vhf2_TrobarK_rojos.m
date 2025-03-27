%Valor r0
setr0(3.964)
setL1(0.07965)
setL3(0.3105)
load("dades.mat","m_c","m_b")
setm_q(m_c)
setspin(1)

% Ajust amb estats rojos només

%Totes les E en GeV
%Vector de valors de E calculades
a=zeros(14);
%Vector de valors de E teòriques
t(1)=4.0296;
t(2)=3.8976;
t(3)=3.9286;
t(4)=4.0746;

%Vector error de energies teòriques
e(1)=0.0176;
e(2)=0.0186;
e(3)=0.0236;
e(4)=0.0216;

I1=1;
I2=1;
results(I1,I2)=0;
%Programa que busca la k òptima per a la chi^2
for ka1=0.088:0.0001:0.090
    setk1(ka1);
    I2=1;
    for ka2=0.0040:0.0001:0.0045
        setk2(ka2);
        %Vector en los valors de la energia que necesito
        [aux,W,x]=QuarkoniumS0J1(m_q,spin);
        a(1)=aux(1);
        [aux,W,x]=Spin1Jcal0_1(m_q,spin);
        a(2)=aux(1);
        [aux,W,x]=Spin1Jcal1_2(m_q,spin);
        a(3)=aux(1);
        [aux,W,x]=Spin1Jcal2_1(m_q,spin);
        a(4)=aux(1);
    
        chi=0;
        for i=1:4
            chi=chi+((a(i)-t(i))^2)/((e(i))^2);
        end
        I2
        results(I1,I2)=chi;
        I2=I2+1;
        
    end 
    
    I1
    I1=I1+1;
    
end

[m,p]=min(results);
[~,c]=min(m);
f=p(c);

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
function x4=spin
global v4
x4=v4;
end 

function setspin(val4)
global v4
v4 = val4;
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

