function [u,x]=god_e(u0,flux,flux1,t,cfl,ab)
%   ATTENZIONE: questa e' una versione entropica del metodo 
%               di Godunov
%   [u,x]=GOD_E(U0,FLUX,T,CFL)  Calcola la soluzione del 
% problema
%   u_t+(flux(u))_x=0 su [-1,1] con il metodo di Godunov.
%   FLUX e' la funzione di flusso.
%   U0 e' il vettore che contiene le condizioni iniziali,T e'
%   l'istante finale, CFL e' una stima di max|f'(u)|
%  [u,x]=GOD_NED(U0,FLUX,T,CFL,AB)  Calcola la soluzione
%  del problema
%   u_t+(flux(u))_x=0 sull'intervallo AB=[A,B]
%   Le condizioni al contorno sono contenute nella variabile globale
%   BC. I valori possibili sono: 
%   BC='p'  Condizioni al contorno periodiche
%   BC='f'  Condizioni al contorno free-flow
%   LANDA0 (global) e' lo scalare t.c. dt=LANDA0*h/CFL

global bc landa0

a=ab(1)
b=ab(2)
n=length(u0)-2
h=(b-a)/(n)
dt=landa0*h/abs(cfl)
nt=fix(t/dt)+1 % arrotonda (T/DT) all'intero immediatamente superiore
dt = t/nt 
landa=dt/h 


for kt = 1:nt
    % Calcola il flusso: memorizza f(i+1/2) in f(i)
    fl=feval(flux,u0);
    for i=1:n+1
        s = (fl(i+1) - fl(i))*(u0(i+1)-u0(i));
        if s >= 0
            f(i) = fl(i);
        else
            f(i) = fl(i+1);
        end
       % entropy fix: corregge il flusso se c'e' una rarefazione
        % transonica
        fl1=feval(flux1,u0);
        if fl1(i) <0 & fl1(i+1)> 0
            % trova il valore di u per il quale f'(u)=0
            if u0(i) >= u0(i+1)
                us=fzero(flux1,[u0(i+1),u0(i)])
                display('us >0')
            else
                us=fzero(flux1,[u0(i),u0(i+1)]);
                display('us <0')
            end 
            f(i) = feval(flux,us); 
        end
    end    
    % Calcola la soluzione
    for i=2:n+1
        u(i) = u0(i) - landa*(f(i)-f(i-1));
    end
% Calcola le condizioni al contorno
 if bc=='f'
        u(1) = u0(1); u(n+2) = u0(n+2);
    elseif bc=='p'
        u(1) = u(n+1); u(n+2) = u(2);
    else
        display('bc non e'' stato impostato')
        return
    end
    u0=u;   % aggiorna u0 per il prossimo passo
end

x=linspace(a-h/2,b+h/2,n+2);