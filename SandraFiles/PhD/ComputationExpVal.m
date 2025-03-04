%Valor r0
setr0(3.964)
%massa
load("dades.mat","m_c","m_b")
setm_q(m_c)

%values (autentic)
Min=0.27;
Mfin=0.62; %Has to be changed by hand depnending on the trasnition

%Now we compute I_theta for 4s->3s for diferent M
%|Ef-Ei|>M

%number of spaces:
k=50;
%M row:
nM=(Mfin-Min)/k;
M=Min:nM:Mfin;

%spin average de els termes involucrats enla integral
I_if_square=zeros(1,k+1);
I_if_c_square=zeros(1,k+1);
I_if_0c_square=zeros(1,k+1);
I_if_s_square=zeros(1,k+1);

%directament les I_if sense square
I_if0=zeros(1,k+1);
I_if_c0=zeros(1,k+1);
I_if1=zeros(1,k+1);
I_if_c1=zeros(1,k+1);
I_if_s1=zeros(1,k+1);
I_if2=zeros(1,k+1);
I_if_c2=zeros(1,k+1);
I_if_s2=zeros(1,k+1);

n=1;
while n<=k+1
%Inser here the transition following the guide of TransitionsAdded.m
%%%%%%%%%%%%%%%%%%%%%%%%%

%només 1 contribució Ji=0,m=0->Jf=0,m'=0
%al spin average contribueixen els termes I^2, I_c^2 i I*I_c

%canviar
%estats inciials i finals
Ni=n;
Nf=n';

transition=FormFactor_ItoF('QQS0toS0_F/');
transitionC=FormFactor_ItoF('QQS0toS0_Fc');

ExpValFunc=ExpValFunctions('QQ');

I_if(n)=ExpValFunc(Ni,Nf,M(n),transition,0,0);
I_if_c(n)=ExpValFunc(Ni,Nf,M(n),transition_c,0,0);

%el spin average = com nomes hi ha 1 m possible per a la transicio
%directament per a cada

I_if_square(n)=I_if(n)^2;
I_if_c_square(n)=I_if_c(n)^2;
I_if_0c_square(n)=I_if(n)*I_if_c(n);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Number of iteration: %d / %d\n', n, k);
    n=n+1;

end


%Aquetses variables es passen a fitxers notebook
Mtrans=M';
Itrans=I_if_square';
Itrans_c=I_if_c_square';
Itrans_0c=I_if_0c_square';
Itrans_s=I_if_s_square';

% The following code is to create .txt files where the variables are stored
baseName = 'QQcharm_3s-2s';
folderPath = '/Users/sandra/Documents/Doctorat/Projectes\ PhD/Transicions\ a\ 2\ pions/Ordre\ més\ baix\ del\ Lagrangià/TransitionData';

fileName = fullfile(folderPath, sprintf('%s.txt', baseName));
fileName_c = fullfile(folderPath, sprintf('%s_c.txt', baseName));
fileName_0c = fullfile(folderPath, sprintf('%s_0c.txt', baseName));
fileName_s = fullfile(folderPath, sprintf('%s_s.txt', baseName));

%Make sure the path exists
if ~exist(folderPath, 'dir')
    fprintf('Error: The folder "%s" does not exist.\n', folderPath);
    return; 
end

% Open to write
file = fopen(fileName, 'w');
file_c = fopen(fileName_c, 'w');
file_0c = fopen(fileName_0c, 'w');

% Writhe the heading
fprintf(file,    '%6s %12s\n', 'M', 'I_if^2');
fprintf(file_c,  '%6s %12s\n', 'M', 'I_if_c^2');
fprintf(file_0c, '%6s %12s\n', 'M', 'I_if*I_if_c');

% Wirthe the two columns
for i = 1:length(M)
    fprintf(file,    '%.4f %12.8f\n', M(i), I_if_square(i));
    fprintf(file_c,  '%.4f %12.8f\n', M(i), I_if_c_square(i));
    fprintf(file_0c, '%.4f %12.8f\n', M(i), I_if_0c_square(i));
end


fclose(file);
fclose(file_c);
fclose(file_0c);

%Not all transitions need a _s file so only created for thoes which need it
%(p->p and d->s)
if any(I_if_s_square)
    file_s = fopen(fileName_s, 'w');
    fprintf(file_s, '%6s %12s\n', 'M', 'I_if_s^2');
    for i = 1:length(M)
        fprintf(file_s, '%.4f %12.8f\n', M(i), I_if_s_square(i));
    end
    fclose(file_s);
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
function x7=spin
global v7
x7=v7;
end 

function setspin(val7)
global v7
v7 = val7;
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

