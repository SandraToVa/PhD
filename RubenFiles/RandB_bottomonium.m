%CALUL DEL ESPECTRE BOTTOMONIUM
%Valor r0
setr0(3.964)
setk1(0.035398)
setk2(0.001717)


%Totes les E en GeV
%Vector de valors de E calculades
a=zeros(14);

%Vector en los valors de la energia que necesito
aux=QuarkoniumS0J1(4.8802,1,0);
a(1)=aux(1);
aux=Spin1Jcal0_1(4.8802,0,1);
a(2)=aux(1);
aux=Spin1Jcal1_2(4.8802,1,1);
a(3)=aux(1);
aux=Spin1Jcal2_1(4.8802,2,1);
a(4)=aux(1);
aux=Spin1Jcal0_2(4.8802,0,1);
a(5)=aux(1);
aux=Spin1Jcal2_2(4.8802,2,1);
a(6)=aux(1);
aux=Spin1Jcal2_2(4.8802,2,1);
a(7)=aux(2);
aux=QuarkoniumS0J0(4.8802,0,1);
a(8)=aux(5);
aux=Spin1Jcal1_1(4.8802,1,1);
a(9)=aux(1);
aux=Spin1Jcal1_1(4.8802,1,1);
a(10)=aux(2);
aux=Spin1Jcal1_1(4.8802,1,1);
a(11)=aux(5);
aux=QuarkoniumS0J1(4.8802,1,1);
a(12)=aux(2);
aux=QuarkoniumS0J2(4.8802,2,1);
a(13)=aux(1);
aux=Spin1Jcal3_1(4.8802,3,1);
a(14)=aux(1);
    
a

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