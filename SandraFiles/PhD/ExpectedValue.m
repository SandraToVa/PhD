
function V=ExpectedValue(i,f,Op,initial,final,m,sin,sfin)
%Function to compute espected values using wavefunctions and any observable
%Output:        - Returns a vector of positions = columns of wavefunctions initial
%               and final (if the hybrid computed is a 3x3 matrix, 3 colums, if it is a
%               7x7 matrix, 7 columns)
%               This vector is the expectation value of
%               final wavefunction * operator * initial wavefunction
%Input:         - i=position in the spectum of eneries of the hybrid or
%               quarkonium that we want the inital wavefunction to have
%               - f= same as i but for final wavefucntion
%               - Op= operator in matrix form (it has to be a diagonal
%               matrix) the rowas and columns the same as the
%               wavefunctions. This is the operator that we want to find
%               the expectation value of
%               - initial= initial wavefunction form the scripts of "Spin..."
%               for hybrid and "Quarkonium..." for quarkonium
%               - final= same as initial for final state
%               - m= mass of the quark. Charm or Bottom

%I have to modify the tolerance in "initial" and "final" for the mesh to
%coincide. And for the x and lenght of the 2nd dimension of W to coincide!
setspin(sin)
[Ei,Wi,x]=initial(m,spin);

setspin(sfin)
[Ef,Wf,x]=final(m,spin);

%The vale of our initial state
disp('Initial state i')
i
disp('With energy:')
Ei(i)
disp('Wavefunction:')
Yi=Wi(:,:,i)

%The length of the wavefunction i 
[numRi,numCi] = size(Yi);


%The vale of our final state
disp('Initial final f')
f
disp('With energy:')
Ef(f)
disp('Wavefunction:')
Yf=Wf(:,:,f)

%The length of the wavefunction f 
[numRf,numCf] = size(Yf);


%Now we normalize the wavefunctions
Yi_n= Yi/norm(Yi);
Yf_n = Yf/norm(Yf);
 %He de normalitzar la x, crec q o


%numRi=numRf and numCi=numCf
V=zeros(numRi,numRf);
for N=1:numRi
    %Muliply each row for any row of the other wavefunction normalized
    for M=1:numRf
        Yi_f=Yi_n(N,:); %row vector
        Yf_f=Yf_n(M,:); %row vector

        %With out wavefunctions computed now we multiply for our operator
        %(matrix) Yf*Op*Yi
        %Auxiliar matrix(vector) for the first multiplication. It is a row
        %vector
        V1=Yf_f*Op;
        %Now our second multiplication is for each element in V1 and each
        %element in Yi
        V2=zeros(1,length(V1)); %row of zeros
        for el=1:length(V1)
           %Second auxiliar vector a row vector. Final multiplication of the
           %Yf*Op*Yi
           V2(el)=V1(el)*Yi_f(el);
        end
        %The integrand is exactly V2
        I=V2;
        V(N,M)=trapz(x,I);  
    end
end

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

