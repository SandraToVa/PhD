m_c=1.4702;
m_b=4.8802;
L1=0.07965; %no ho uso pq va massa lento
%mínima=0.059
%màxima=0.07965  %OPTIMA

L3=0.3105; %ho definixo al codi que executo
%mínima=0,23
%màxima=0.3105  %OPTIMA

CF_c=1.12155;
CF_b=0.87897;
kappa=0.187;

r0=3.964;

%Faltairia posar los valors dels papers de l'espectre de c per tindre'ls


save("dades.mat","m_c","m_b","kappa","CF_b","CF_c","L3","L1")