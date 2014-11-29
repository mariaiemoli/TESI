function u0=rieman(ul,ur,nx,x,x_r)
% U0=INIT(UX,N) calcola la condizione iniziale per
%       i metodi per le equazioni iperboliche, applicando
%       le condizioni al contorno della variabile globale BC
%       Usa N intervalli su [-1,1], con U0(I)=UX(X(I)).
% U0=INIT(UX,N,AB) calcola la condizione iniziale 
%       sull'intervallo AB = [A,B].

global bc

a=x(1);
b=x(nx);

h=(b-a)/(nx);

if nargin < 5
    x_r=0;
end

for i=1:nx
    if x(i) <= x_r
        u0(i)=ul;
    else
        u0(i)=ur;
    end
end

