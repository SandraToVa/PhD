%I want to compute the hyperfine spectrum of Quarkonium (no hybrid) at long
%distances. First I compute it as in the paper of Rubern Oncala

%En los codis de QuarkoniumS0J#.m calculo los estats de energia.
%En el paper on posa L es la J/j dels codis. Per tant els estats #s se
%calculen en QuarkoniumS0J0.m; els #p en QuarkoniumS0J1.m; etc. S'ha de
%reproduir la Taula V quan se usa lo potencial sencer.

%El potencial que usem es el Eq. (A2) que es el Vg en els codis.

%A llargues distàncies la part de (A2) que es k/r s'ha de llevar i ens
%quedem en lo altre.

load("dades.mat","m_b","m_c")
setm_q(m_c)
setspin(0) %for quarkonium s=0

[aux_s,mesh,syst]=QuarkoniumS0J0(m_q,spin);
[aux_p,mesh,syst]=QuarkoniumS0J1(m_q,spin);
[aux_d,mesh,syst]=QuarkoniumS0J2(m_q,spin);


%Espectre
a_s=aux_s;
a_p=aux_p;
a_d=aux_d;

%Diferencia
d_s=zeros(6);
d_p=zeros(6);
d_d=zeros(6);

for i=1:1:5
    d_s(i)=a_s(i+1)-a_s(i);
    d_p(i)=a_p(i+1)-a_p(i);
    d_d(i)=a_d(i+1)-a_d(i);
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
function x4=spin
global v4
x4=v4;
end 

function setspin(val4)
global v4
v4 = val4;
end

