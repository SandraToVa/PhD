% Codi guardar en un fitxer txt les wf

%%Valor r0
setr0(3.964)
setL1(0.059)
setL3(-0.230)
load("dades.mat","m_c","m_b")
setm_q(m_c)
setspin(1)
% l= interpolació (0), llagures distncies (1), bad long distances (2)
setl(0)

%Primer calculem l'espectre en la A i B òptimes - f(i) -
setk1(-0.0848)
setk2(0.0021)

% Initialize cell array to store variable-sized matrices
wf = {};
x_arrey = {};
%(s/d)1
[aux,auxwf,x]=QuarkoniumS0J1(m_q,spin);
f(1)=aux(1);
wf{1}=auxwf(:,:,1);
x_arrey{1}=x;
f(5)=aux(2);  %p1
wf{5}=auxwf(:,:,2); %p1
x_arrey{5}=x;
[aux,auxwf,x]=Spin1Jcal0_1(m_q,spin);
f(2)=aux(1);
wf{2}=auxwf(:,:,1);
x_arrey{2}=x;
[aux,auxwf,x]=Spin1Jcal1_2(m_q,spin);
f(3)=aux(1);
wf{3}=auxwf(:,:,1);
x_arrey{3}=x;

[aux,auxwf,x]=Spin1Jcal2_1(m_q,spin);
f(4)=aux(1);
wf{4}=auxwf(:,:,1);
x_arrey{4}=x;
%p1
[aux,auxwf,x]=Spin1Jcal0_2(m_q,spin);
f(6)=aux(1);
wf{6}=auxwf(:,:,1);
x_arrey{6}=x;
[aux,auxwf,x]=Spin1Jcal1_1(m_q,spin);
f(7)=aux(1);
wf{7}=auxwf(:,:,1);
x_arrey{7}=x;
f(10)=aux(2); %(p/f)2
wf{10}=auxwf(:,:,2); %(p/f)2
x_arrey{10}=x;
% p0
if m_q==1.4702
    f(14)=aux(3); %If charmonium
    wf{14}=auxwf(:,:,3);
end
if m_q==4.8802
    f(14)=aux(5); %If bottomium
    wf{14}=auxwf(:,:,5);
end
x_arrey{14}=x;
[aux,auxwf,x]=Spin1Jcal2_2(m_q,spin);
f(8)=aux(1);
wf{8}=auxwf(:,:,1);
x_arrey{8}=x;
f(11)=aux(2); %(p/f)2
wf{11}=auxwf(:,:,2); %(p/f)2
x_arrey{11}=x;
%(p/f)2
[aux,auxwf,x]=QuarkoniumS0J2(m_q,spin);
f(9)=aux(1);
wf{9}=auxwf(:,:,1);
x_arrey{9}=x;
[aux,auxwf,x]=Spin1Jcal3_1(m_q,spin);
f(12)=aux(1);
wf{12}=auxwf(:,:,1);
x_arrey{12}=x;
%p0
[aux,auxwf,x]=QuarkoniumS0J0(m_q,spin);
f(13)=aux(1);
wf{13}=auxwf(:,:,1);
x_arrey{13}=x;

disp('Wave funtions computed');

%%

% Assume wf is a cell array where each wf{n} is a (3 or 4, variable length) matrix
save_path = '/Users/sandra/Downloads'; % Change this to your desired folder
filename = fullfile(save_path, 'wf_a_data.txt'); % Combine path and filename
fid = fopen(filename, 'w'); % Open file for writing

for n = 1:14  % Loop through each state (third dimension)
    data = wf{n}; % Extract the n-th matrix
    [rows, cols] = size(data); % Get its size
    fprintf(fid, 'State %d: Size [%d x %d]\n', n, rows, cols); % Write header
    
    for i = 1:rows
        fprintf(fid, '%f ', data(i, :)); % Print row values space-separated
        fprintf(fid, '\n'); % New line after each row
    end
    
    fprintf(fid, '\n'); % Extra newline to separate states
end

fclose(fid); % Close the file
disp('Data saved successfully!');


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