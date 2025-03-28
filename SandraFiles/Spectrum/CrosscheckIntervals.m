% Codi per comprovar si la manera de trobar intervals de confiança es
% correcta

%In IntervalsConfian.m we computed the spectrum with E_i=c_i+b_iB+a_iA
%with A p10 & B p01 

%Now we want to compute the a_i and b_i modifing the hamiltonian (aka the
%Vhf and Vhf2) using a Mecanical Quantica theorem: c_i + a_iA + b_iB =
%<funcions d'ona amb A i B optimes | operator of H partial A or B | funcions
% d'ona amb A i B optimes>; aquest nou amiltonia es la Vhf derivada respecta
%A (per a a_i) i la Vhf2 derivada respecte B (per a b_i)

%For each state I compute the wave functions of the state with the normal
%hiperfine hamiltonian
%For the corresponding state I create a matrix with the partial of the
%hyperfine potential multipied by the matrix (8) and (9) in the paper
% For p1 and (p/f)2 (2+-) the result is very different we see that they are
% very mixed so instead of the partial matrix i used the full matrix
% including both states in both of them

% CONCLUSION: the reasults are not similar enough. Linearizing the problem
% was an error. We need to compute the errors in another way. Finding A and
% B for lattice data +- lattice errors.


setr0(3.964)
setL1(0.059)
setL3(-0.230)
load("dades.mat","m_c","m_b")
setm_q(m_c)
setspin(1)
% l= interpolació (0), llagures distncies (1), bad long distances (2)
setl(0)

%Primer calculem l'espectre en la A i B òptimes
setk1(-0.0848)
setk2(0.0021)

% Initialize cell array to store variable-sized matrices
wf = cell(1, 14);         %for each cell there will be a matrix(wf made by diferent states at diferent x positions)
Hmatrix_a = cell(1, 14);    %for each cell an array of matrices of the H marix at each position x formed by Hvector
Hmatrix_b = cell(1, 14);
xvector = cell(1, 14);
E = zeros(1,14);
final_a = zeros(1,14);
final_b = zeros(1,14);

% All of s=0 states do not depend on Vhf nor Vhf2 therefore their
% derivative with respect to A or B of the potencials is 0 and the hmatrix
% is just zero

[aux,auxwf,x]=QuarkoniumS0J1(m_q,spin);
% 1--
E(1)=aux(1);
wf{1} = auxwf(1:2, :, 1);
xvector{1}=x;
Hvector_a = zeros(2,2,length(x));
Hvector_b = zeros(2,2,length(x));
Hmatrix_a{1} = Hvector_a;
Hmatrix_b{1} = Hvector_a;
[final_a(1),final_b(1)]=computeABstate(Hvector_a,Hvector_b,wf{1});

% 1++
E(5)=aux(2);  
wf{5}=auxwf(3, :, 2);
xvector{5}=x;
Hvector_a = zeros(1,1,length(x));
Hvector_b = zeros(1,1,length(x));
Hmatrix_a{5} = Hvector_a;
Hmatrix_b{5} = Hvector_a;
[final_a(5),final_b(5)]=computeABstate(Hvector_a,Hvector_b,wf{5});

% 2++
[aux,auxwf,x]=QuarkoniumS0J2(m_q,spin);
E(9)=aux(1);
wf{9}=auxwf(1:2, :, 1);
xvector{9}=x;
Hvector_a = zeros(2,2,length(x));
Hvector_b = zeros(2,2,length(x));
Hmatrix_a{9} = Hvector_a;
Hmatrix_b{9} = Hvector_a;
[final_a(9),final_b(9)]=computeABstate(Hvector_a,Hvector_b,wf{9});

% 0 ++
[aux,auxwf,x]=QuarkoniumS0J0(m_q,spin);
E(13)=aux(1);
wf{13}=auxwf(1, :, 1);
xvector{13}=x;
Hvector_a = zeros(1,1,length(x));
Hvector_b = zeros(1,1,length(x));
Hmatrix_a{13} = Hvector_a;
Hmatrix_b{13} = Hvector_a;
[final_a(13),final_b(13)]=computeABstate(Hvector_a,Hvector_b,wf{13});

% State 0-+ (s/d)1
[aux,auxwf,x]=Spin1Jcal0_1(m_q,spin);
E(2)=aux(1);
wf{2}=auxwf(:,:,1);
xvector{2}=x;
aux_matrix = [4, 0; 0, -2];     %matrix that multiplies Vhf
aux_matrix2 = -2*[0, -sqrt(2)/3; -sqrt(2)/3, -2/3];   %matrix that multiplies Vhf2
[Hmatrix_a{2}, Hmatrix_b{2}] = computeHmatrix(2,x,aux_matrix,aux_matrix2);
[final_a(2),final_b(2)]=computeABstate(Hmatrix_a{2},Hmatrix_b{2},wf{2});

% State 1-+ (s/d)1
[aux,auxwf,x]=Spin1Jcal1_2(m_q,spin);
E(3)=aux(1);
%wf{3}=auxwf(:, :, 1);
wf{3}=auxwf(1:2, :, 1);
xvector{3}=x;  
aux_matrix = [2, 0; 0, -1];    %matrix that multiplies Vhf
aux_matrix2 = -2*[0, (sqrt(1/2)-sqrt(2))/3;(sqrt(1/2)-sqrt(2))/3, 2/3];   %matrix that multiplies Vhf2
%aux_matrix = [2, 0, 0; 0, -1, -sqrt(3); 0, -sqrt(3), 1];
%aux_matrix2 = -2*[0, (sqrt(1/2)-sqrt(2))/3, sqrt(1/6);(sqrt(1/2)-sqrt(2))/3, 2/3, 0; sqrt(1/6), 0, -1/3];
[Hmatrix_a{3}, Hmatrix_b{3}] = computeHmatrix(2,x,aux_matrix,aux_matrix2);
[final_a(3),final_b(3)]=computeABstate(Hmatrix_a{3},Hmatrix_b{3},wf{3});

% State 2-+ (s/d)1
[aux,auxwf,x]=Spin1Jcal2_1(m_q,spin);
E(4)=aux(1);
wf{4}=auxwf(1:2, :, 1);
%wf{4}=auxwf(:, :, 1);
xvector{4}=x;
aux_matrix = [-2, 0; 0, 1];     %matrix that multiplies Vhf
aux_matrix2 = -(2/3)*[0, 1;1, 1];  %matrix that multiplies Vhf2
%aux_matrix = [-2, 0, 0, 0, 0; 0, 1, -3*sqrt(3/5), 0, 0; 0, -3*sqrt(3/5), 1/3, -(4/3)*sqrt(7/5), 0; 0, 0,  -(4/3)*sqrt(7/5), 8/3, 0; 0, 0, 0, 0, 1];
%aux_matrix2 = -(2/3)*[0, 1, 0, 0, 0;1, 1, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0];
[Hmatrix_a{4}, Hmatrix_b{4}] = computeHmatrix(2,x,aux_matrix,aux_matrix2);
[final_a(4),final_b(4)]=computeABstate(Hmatrix_a{4},Hmatrix_b{4},wf{4});

disp('States (s/d)1 computed');

% State 0+- p1
[aux,auxwf,x]=Spin1Jcal0_2(m_q,spin);
E(6)=aux(1);
wf{6}=auxwf(:, :, 1);
xvector{6}=x;
aux_matrix = 2;    %matrix that multiplies Vhf
aux_matrix2 = 4/3;   %matrix that multiplies Vhf2
[Hmatrix_a{6}, Hmatrix_b{6}] = computeHmatrix(1,x,aux_matrix,aux_matrix2);
[final_a(6),final_b(6)]=computeABstate(Hmatrix_a{6},Hmatrix_b{6},wf{6});

% State 1+- p1
% Per que doni millor s'hauria d'incloure el mixing en (p/f)2 i p0 (1+-)
[aux,auxwf,x]=Spin1Jcal1_1(m_q,spin);
E(7)=aux(1);
wf{7}=auxwf(2, :, 1);
xvector{7}=x;
aux_matrix = 1;     %matrix that multiplies Vhf
aux_matrix2 = 2/3;   %matrix that multiplies Vhf2
[Hmatrix_a{7}, Hmatrix_b{7}] = computeHmatrix(1,x,aux_matrix,aux_matrix2);
[final_a(7),final_b(7)]=computeABstate(Hmatrix_a{7},Hmatrix_b{7},wf{7});

% State 2+- p1
% Incloen lo mixing en (p/f)2 (2+-)
[aux,auxwf,x]=Spin1Jcal2_2(m_q,spin);
E(8)=aux(1);
wf{8}=auxwf(1:3, :, 1);
xvector{8}=x;
aux_matrix = [-1, -sqrt(3), 0; -sqrt(3), 1, 0; 0, 0, -2/3];      %matrix that multiplies Vhf
aux_matrix2 = [-2/3, 0, 0; 0, 0, 0; 0, 0, 0];   %matrix that multiplies Vhf2
[Hmatrix_a{8}, Hmatrix_b{8}] = computeHmatrix(3,x,aux_matrix,aux_matrix2);
[final_a(8),final_b(8)]=computeABstate(Hmatrix_a{8},Hmatrix_b{8},wf{8});

disp('States p1 computed');

% State 1+- (p/f)2
% Per que doni millor s'hauria d'incloure el mixing en p1 i p0 (1+-)
[aux,auxwf,x]=Spin1Jcal1_1(m_q,spin);
E(10)=aux(2);
wf{10}=auxwf(3:4, :, 2);
xvector{10}=x;
aux_matrix = [3, 0; 0, -2];    %matrix that multiplies Vhf
aux_matrix2 = -(2/5)*[1, -sqrt(3/2);-sqrt(3/2), -8/3];   %matrix that multiplies Vhf2
[Hmatrix_a{10}, Hmatrix_b{10}] = computeHmatrix(2,x,aux_matrix,aux_matrix2);
[final_a(10),final_b(10)]=computeABstate(Hmatrix_a{10},Hmatrix_b{10},wf{10});

% State 2+- (p/f)2
% Incloen lo mixing en p1 (2+-)
[aux,auxwf,x]=Spin1Jcal2_2(m_q,spin);
E(11)=aux(2);
wf{11}=auxwf(1:3, :, 2);
xvector{11}=x;
aux_matrix = [-1, -sqrt(3), 0; -sqrt(3), 1, 0; 0, 0, -2/3];    %matrix that multiplies Vhf
aux_matrix2 = -(2/5)*[0, 0, 0; 0, 1/3, sqrt(2/3)-sqrt(3/2);0, sqrt(2/3)-sqrt(3/2), -8/9];  %matrix that multiplies Vhf2
[Hmatrix_a{11}, Hmatrix_b{11}] = computeHmatrix(3,x,aux_matrix,aux_matrix2);
[final_a(11),final_b(11)]=computeABstate(Hmatrix_a{11},Hmatrix_b{11},wf{11});

% State 3+- (p/f)2
[aux,auxwf,x]=Spin1Jcal3_1(m_q,spin);
E(12)=aux(1);
wf{12}=auxwf(1:2, :, 1);
xvector{12}=x;
aux_matrix = -2*[1, 0; 0, -2/3];    %matrix that multiplies Vhf
aux_matrix2 = -(2/5)*[-2/3, sqrt(2/3);sqrt(2/3), 16/9];  %matrix that multiplies Vhf2
[Hmatrix_a{12}, Hmatrix_b{12}] = computeHmatrix(2,x,aux_matrix,aux_matrix2);
[final_a(12),final_b(12)]=computeABstate(Hmatrix_a{12},Hmatrix_b{12},wf{12});

disp('States (p/f)2 computed');

% State 1+- p1
% Per tal que no doni zero s'hauria d'incloure el mixing en p1 i (p/f)2
% (1+-)
[aux,auxwf,x]=Spin1Jcal1_1(m_q,spin);
E(14)=aux(3);
wf{14}=auxwf(1, :, 3);
xvector{14}=x;
aux_matrix = 0;     %matrix that multiplies Vhf
aux_matrix2 = 0;   %matrix that multiplies Vhf2
[Hmatrix_a{14}, Hmatrix_b{14}] = computeHmatrix(1,x,aux_matrix,aux_matrix2);
[final_a(14),final_b(14)]=computeABstate(Hmatrix_a{14},Hmatrix_b{14},wf{14});

disp('States p0 computed');


function [Hmatrix_a,Hmatrix_b]=computeHmatrix(dim,x,aux_matrix,aux_matrix2)
%Input:     dim = dimension of the rellevant! wf =1 for p_n and d_n states
%                 and 2 for (s/d) and (p/f)
%           x = array of meshpoints in where the matrix are computed
%           aux_matrices = matrices corresponding to eq (8) and (9) of the
%                          paper for each state have to be defined
%
%Output:    Hmatrices = matrices corresponding to the partial of A and B of
%                       the hamiltonian that will be multiplied only for the
%                       rellevant wf

    Hvector_a = zeros(dim,dim,length(x));
    Hvector_b = zeros(dim,dim,length(x));
    Vhf = zeros(1,length(x));
    Vhf2 = zeros(1,length(x));
    for x_el=1:length(x)
        Vhf(x_el) = 1 / (1 + (x(x_el) / r0)^5);
        Vhf2(x_el) = x(x_el)^2 / (1 + (x(x_el) / r0)^7);
        Hvector_a(:,:,x_el) = aux_matrix .* Vhf(x_el);
        Hvector_b(:,:,x_el) = aux_matrix2 .* Vhf2(x_el);
    end 
    Hmatrix_a = Hvector_a;
    Hmatrix_b = Hvector_b;
end

function [a,b]=computeABstate(Hmatrix_a,Hmatrix_b,wf)
%Input:     Hmatrices = matrices that will be sanwihed between the
%                       wavefunctions
%           wf = only the rellevant row of the wave function that will be
%                multiplied in front and behind the Hmatrix
%
%Output:    a,b = values of the a and b of each state

    [rows,cols] = size(wf);

    % Initialize accumulators
    sum_Ha_weighted = 0;
    sum_Hb_weighted = 0;
    sum_wf_norm = 0;

    % Loop through each column
    for i = 1:cols  % Iterate over 'cols' dimension
        wf_t = wf(:, i)';  % Transpose wf(:, i) to be a row vector (1xN)
        H_a = Hmatrix_a(:, :, i);  % Extract the corresponding H matrix (NxN)
        H_b = Hmatrix_b(:, :, i);  % Extract the corresponding H matrix (NxN)
    
        % Compute wf_t * Hmatrix * wf
        term1_a = wf_t * H_a * wf(:, i);
        term1_b = wf_t * H_b * wf(:, i);
        sum_Ha_weighted = sum_Ha_weighted + term1_a;  % Accumulate sum
        sum_Hb_weighted = sum_Hb_weighted + term1_b;
    
        % Compute wf_t * wf (without H)
        term2 = wf_t * wf(:, i);
        sum_wf_norm = sum_wf_norm + term2;  % Accumulate sum
    end

    % Compute the final ratio
    a = sum_Ha_weighted / sum_wf_norm;
    b = sum_Hb_weighted / sum_wf_norm;

end


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