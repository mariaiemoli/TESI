function u0=init(ux,n,ab)
% U0=INIT(UX,N) calcola la condizione iniziale per
%       i metodi per le equazioni iperboliche, applicando
%       le condizioni al contorno della variabile globale BC
%       Usa N intervalli su [-1,1], con U0(I)=UX(X(I)).
% U0=INIT(UX,N,AB) calcola la condizione iniziale 
%       sull'intervallo AB = [A,B].

global bc
if nargin < 3
    ab=[-1,1];
end
a=ab(1); b=ab(2); h=(b-a)/(n);
x=a-h/2:h:b+h/2; 
u0=feval(ux,x);
% Aggiunge le condizioni al contorno
if bc=='f'
   u0(1)=u0(2); u0(n+2)=u0(n+1);
elseif bc=='p'
   u0(1)=u0(n+1); u0(n+2)=u0(2);
else
   display('bc non e'' stato impostato')
   return
end