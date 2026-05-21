clc;
clear;
addpath("helper");
%%%% Progetto 1 Calcolo Numerico %%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Analisi con metodi diretti %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
disp("%%% Analisi con metodi diretti %%%")
disp("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
%% Punto 1
% Definire la matrice delle ammettenze A ed il termine noto b
Gin = 1 / 60;
Vin = 15;
% Rk = Rin / 2k
% Gk = 2k * Gin
G1 = 2 * 1 * Gin;
G2 = 2 * 2 * Gin;
G3 = 2 * 3 * Gin;
G4 = 2 * 4 * Gin;
G5 = 2 * 5 * Gin;
G6 = 2 * 6 * Gin;
disp("Punto 1")
format short e
disp("Matrice delle Ammettenze:")

A = [G1+G2+Gin, -G1, 0 , -G2, 0, 0, 0, 0, 0; %Riga 1
     -G1, 2*G1 + G2, -G1, 0, -G2, 0 , 0 , 0 , 0;%Riga 2
     0, -G1, G2+G1, 0, 0, -G2, 0, 0 , 0;%Riga 3
     -G2, 0, 0, G4+G3+G2, -G3, 0, -G4, 0, 0;%Riga 4
     0, -G2, 0, -G3, G4 + 2*G3 + G2, -G3, 0, -G4, 0; %Riga 5
     0, 0, -G2, 0, -G3, G4+G3+G2, 0, 0, -G4;
     0, 0, 0, -G4, 0, 0, G6+G5+G4, -G5, 0;
     0, 0, 0, 0, -G4, 0, -G5, G6+ 2*G5+ G4, -G5;
     0, 0, 0, 0, 0, -G4, 0, -G5, G4+ G5 + G6;]
disp("Termine Noto:")
b = [Vin * Gin; 0; 0; 0; 0; 0; 0; 0; 0]

%% Punto 2
disp("Punto 2")
n = size(A,1);
format short e
% Eseguire il comando format short e. Verificare l'esistenza ed unicit`a
% della fattorizzazione LU di A con Lii = 1, i= 1,...,9.

% Dalla teoria so che la fattorizzazione LU senza pivoting esiste ed è
% unica sse tutti i determinanti delle sottomatrici princiapli sono diversi
% da zero.
% Ho anche delle condizioni sufficienti che mi garantiscono applicabilità
% del meg per ottenere fattorizzazione LU, ovvero matrice a dominanza
% diagonale per righe o per colonne, oppure matrice simmetrica e definita
% positiva

% Verifico se la matrice è simmetrica e definita positiva
if ( A == A')
    if( all(eig(A) > 0))
        disp("La matrice è simmetrica e definita positiva, fattorizzazione" + ...
            " LU esiste ed è unica!")
    end
end
% In questo caso avrei già trovato una condizione sufficiente per
% esistenza ed unicità della fattorizzazione LU, per completezza verifico
% anche le altre condizioni.

% verifico se vale la dominanza diagonale stretta per righe o colonne
rowDominance = 1;
for i = 1 : n
    rowSum = sum(abs(A(i, :))) - abs(A(i, i));
    if abs(A(i, i)) <= rowSum
        rowDominance = 0;
        disp("Matrice non a dominanza diagonale di riga stretta")
        break;
    end
end

% Non vale la condizione di dominanza di riga stretta, provo per colonne
columnDominance = 1;
for j = 1 : n
    columnSum = sum(abs(A(:, j))) - abs(A(j, j));
    if abs(A(j, j)) <= columnSum
        columnDominance = 0;
        disp("Matrice non a dominanza diagonale di colonna stretta")
        break;
    end
end

% Una delle condizioni sufficienti regge, se cosi non fosse dovrei passare
% per la condizione necessaria e sufficiente, ovvero
% Calcolare i determinanti delle sottomatrici principali
% Questa verifica è molto onerosa in termini di calcoli, qui è riportata
% solo per completezza
determinants = zeros(n, 1);
for k = 1:n
    subMatrix = A(1:k, 1:k);
    determinants(k) = det(subMatrix);
    if(determinants(k) == 0)
        disp("Determinante della sottomatrice uguale a zero")
    end
end

%Le ipotesi sono vere pertanto posso concludere che esiste ed è unica la
%fattorizzazione LU, stampo i determinanti delle sottomatrici principali
disp("Determinanti delle sottomatrici principali:");
disp(determinants);

disp("Tutti diversi da zero, fattorizzazione LU esiste ed è unica!")
%% Punto 3
disp("Punto 3")
% Calcolare la fattorizzazione LU utilizzando il comando lu di Matlab.
%  Riportare il determinante della matrice U e verificare se ´e stato
%  eseguito il pivoting.

[L, U, P] = lu(A); % Calcolare la fattorizzazione LU
detU = det(U); % Determinante della matrice U
% Se la matrice di permutazine P è l'identità non ho eseguito nessuno swap
% di riga, pertanto nessun pivoting
if isequal(P, eye(n))
    disp("Nessun pivoting eseguito.");
else
    disp("Pivoting eseguito.");
end
disp("Determinante della matrice U:");
disp(detU);
% Nessun pivoting eseguito.
% Determinante della matrice U: 1.5599e-06
%% Punto 4
disp("Punto 4")
% Eseguire il comando format long e. Si consideri la soluzione del sistema
% Ax= b. (1)
% Utilizzare la fattorizzazione LU di A calcolata
% precedentemente per risolvere il sistema (1) con
% le funzioni fwsub.m e bksub.m utilizzate nel Lab 2.
%  Memorizzare nel vettore xc la soluzione
% del sistema triangolare superiore.
% Riportare la norma euclidea del xc cos`ı ottenuto

format long e

% Ricordo che Ax = b sse (Ly = b) and (Ux = y)
%Procedo pertanto a risolvere il sistema (1) in due step, il primo è quello
% di risolvere il sistema triangolare inferiore con sostituzione in avanti
% Ly = b, per farlo usa la funzione fwsub.m

y = fwsub(L, b);

% Posso procedere a risolvere il sistema triangolare superiore Ux = y
% con la sostituzione all'indietro tramite la funzione bckwsub.m
xc = bckwsub(U,y);

% Calcolo la norma euclidea di xc e la stampo
normXc = norm(xc); % Calcolo della norma euclidea di xc
disp("Norma euclidea di xc:");
disp(normXc);
% Norma euclidea di xc: 3.918021398826541e+00

%% Punto 5
disp("Punto 5")
% Risolvere il sistema (1) utilizzando il comando \di Matlab
% e memorizzare nel vettore xm il
% risultato ottenuto. Riportare la norma euclidee del vettore xm.

%Risolvo con comando built in di matlab
xm = A \ b;
normXm = norm(xm);
disp("Norma euclidea di xm:");
disp(normXm);

%% Punto 6
format short e
disp("Punto 6")
% Eseguire il comando format short e. Calcolare la norma
% infinito della diﬀerenza tra xm e
% xc e riportare il risultato ottenuto.

%Ricordo che la norma infinito ritorna il sup delle compoenti del vettore

InfNorm = norm(xm-xc, inf);
disp("Norma infinito della differenza di xm-xc")
disp(InfNorm);

% Noto che la norma infinto della differenza è dell'ordine dell'epsilon
% macchina, lo posso quindi considerare un risultato più che accettabile
% quello ottenuto tramite fattorizzaazione LU

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Analisi con metodi iterativi %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
disp("%%% Analisi con metodi iterativi %%%")
disp("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")

%% Punto 1
disp("Punto 1")
% Verificare che la matrice A `e simmetrica e definita positiva utilizzando
% gli opportuni comandi
% Matlab. Riportare qui sotto l'esito della verifica e i comandi utilizzati.

% Per verifica che A è simmetrica e definita positiva uso il confronto
% element wise con la trasposta e successivamente calcolo gli autovalori
% verificando che siano tutti strettamente maggiori di zero, come da
% definzione

if( A == A')
    if(all(eig(A)>0))
        disp("La matrice è simmetrica e definita positiva")
    else
        disp("La matrice non è SDP")
    end

end
% La matrice risulta effettivamente simmetrica e definita positiva

%% Punto 2
disp("Punto 2")
% Eseguire il comando format long e. Si ponga x0=zeros(9,1), toll=1e-11 e
% nmax=1000. Si utilizzi la function jacobi.m utilizzata nel Lab 3 per
% risolvere il sistema (1) con il metodo di Jacobi e si memorizzi la
% soluzione calcolata nel vettore xJ. Riportare il numero di iterazioni
% kJ eﬀettuate per raggiungere la precisione prescritta.
format long e

x0 = zeros(9,1);
toll = 1e-11;
nmax=1000;

%Risolvo con jacobi
[xJ,kJ] = jacobi(A,b,x0,toll,nmax);

%Riporto il numero di iterazioni
disp("Numero di iterazioni del metodo di jacobi:")
disp(kJ);

%% Punto3
disp("Punto 3")
% Eseguire il comando format short e. Si calcoli e si riporti
% l'errore relativo commesso
format short e

%Riporto errore relativo
disp("True Relative Error")
true_rel_err_j = norm(xm-xJ) / norm(xm)

% true_rel_err = 9.7672e-12
% Noto che avevamo richiesto una precisione dell'ordine di 1e-11 rispetto
% al residuo relativo, so che vale la relazione di maggiorazione:
% true_rel_err <=  K(A) * res_rel e mi torna
% con il risultato numerico

% True Relative Error
%
% true_rel_err_j =
%
%    9.7672e-12

%% Punto 4
disp("Punto 4")
% Si calcoli il raggio spettrale ρJ della matrice di iterazione BJ
% utilizzata al punto precedente
% e si confronti l'errore relativo precedentemente calcolato con ρ kJ
% J , dove kJ `e il numero di
% iterazioni eﬀettuate dal metodo di Jacobi.
% Riportare ρJ e motivare i risultati ottenuti alla
% luce della teoria.


% calcolo il raggio spettrale di BJ definito come il massimo autovalore della
% matrice di iterazione di jacobi
% Calcolo il raggio spettrale della matrice di iterazione di jacobi, parto
% dalla scomposizione di A = D-E(lower) - F(upper)
D = diag(diag(A));
E = tril(A, -1);
F = triu(A, 1);

% Bj = D \ (E+F) = D \ (D - A) = (I - D\A)

BJ = (eye(n) - D\A);

% Calcolo il raggio spettrale della matrice di iterazione definito come il
% sup degli autovalori
rhoJ = max(abs(eig(BJ)));
disp("Raggio spettrale ρJ della matrice di iterazione di Jacobi:");
disp(rhoJ);

% Confronto errore relativo ottenuto con il raggio spettrale di Jacobi
% Parto dalla derivazione del bound sull'errore all'iterazione k, dove
% e_k = xm - xk.
% Per il metodo di Jacobi vale e_k = BJ^k * e_0.
% In disuguaglianza vale norm(e_k) <= norm(BJ^k) * norm(e_0).
% Inoltre, per k grande, il comportamento di norm(BJ^k) è governato dal
% raggio spettrale rhoJ(So che vale per k grande che norm(Bj^k)^ 1/k ~ rhoj)
% Pertanto, asintoticamente per k grande, posso scrivere
% norm(e_k) <= rhoJ^k * norm(e_0).
% Noto però che e_0 = xm - x0 = xm.
% Inoltre norm(xm - xk) / norm(xm) = true_rel_err_j.
% Quindi per k grande true_rel_err_j <= rhoj^k
% Per questo motivo confronto true_rel_err_j con rhoJ^kJ.


% Mostro il fatto che regge la discussione teorica a livello numerico
disp("Confronto tra errore relativo e raggio spettrale elevato a kJ:");
disp(true_rel_err_j)
disp(rhoJ^kJ)
%Infatti numericamnte ottengo che true_rel_err_j <= rhoj^k

%% Punto 5
disp("Punto 5")
% Ripetere i tre punti precedenti utilizzando il metodo del gradiente
% e x0, toll e nmax come sopra. Per approssimare la soluzione del sistema
% (1) si utilizzi la funzione graddyn.m utilizzata
% nel Lab 4 e si memorizzi la soluzione calcolata nel vettore xG.
% Si riporti il numero di iterazioni
% kG eﬀettuate per raggiungere la precisione
% prescritta e l'errore relativo commesso

%Procedo con i tre punti di prima
format long e
x0 = zeros(9,1);
toll = 1e-11;
nmax=1000;
%Risolvo con metodo del gradiente (richardson dinamico)
[xG, kG, res_rel_gdy] = graddyn(A,b,x0,nmax,toll);
%Riporto il numero di iterazioni
disp("Numero di iterazioni del metodo di richardson dinamico:")
disp(kG);
% Riporto errore relativo con questo metodo
format short e
disp("Errore relativo:")
true_relative_error_g = norm(xm -xG) / norm(xm)

% Procedo a calcolare la matrice d'iterazione del metodo di richardson
% stazionario rispetto al parametro alpha ottimale.
% Considero la matrice di
% itereazione di richardson stazionario perchè non dipende dall'iterazione
% come invece accade con il metodo di richardson dinamico. So che una volta
% scelto il parametro ottimale le prestazioni di tale metodo hanno un bound
% teorico che è dettato dalle caratterstiche stesse della matrice, ovver
% dipende dal suo condizionamento. % Quello che mi aspetto di ottenere con
% richardson dinamico è un comportamento confrontabile con quello del
% richardson stazionario con parametro ottimale alpha. In particolare preso
% il parametro ottimale alpha definito come alpha = 2 / (lambda_max(A) + lambda_min(A)),
% la sua matrice di iterazione ha raggio spettrale pari a
% rho_r = (K(A) - 1) / (K(A) + 1).
% % Usando la stessa idea del punto precedente, confronto true_rel_err_g con
% rho_r^kG
% ottengo che true_rel_err_g <= rho_r^kG, e dal momento che tale raggio
% spettrale fa riferiment a richardson stazionario, tale bound sarà
% pessimistico, mi aspetto di fare meglio.
kA = cond(A);
rhoR = (kA - 1 )/ (kA + 1)
% Ora mi aspetto di vedere a livello numerico verificata la relazione

disp("Confronto tra errore relativo e raggio spettrale elevato a kG:");
disp(true_relative_error_g );
disp(rhoR^kG);
%% Punto 6
disp("Punto 6")
%  Si consideri P la matrice diagonale avente i medesimi coeﬃcienti diagonali di A (ossia P
% corrisponde alla parte diagonale di A). Si risolva il sistema lineare (1) applicando il metodo
% del gradiente precondizionato con la matrice P come precondizionatore e medesimi valori
% di x0, toll e nmax. Si utilizzi la funzione gradprec.m utilizzata nel Lab 4 e si memorizzi
% la soluzione calcolata nel vettore xGP. Si riporti il numero di iterazioni kGP eﬀettuate per
% raggiungere la precisione prescritta e l'errore relativo commesso

% Calcolo il precondizionatore
P = diag(diag(A));

% Calcolo la soluzione con il metodo del gradiente precondizionato
[xGP, kGP, res_rel_gpd ] = gradprec(A,b,P,x0,nmax,toll);

% Riporto il numero di iterazioni del metodo di gradiente precondizionato
disp("Numero di iterazioni del metodo di gradiente precondizionato:")
disp(kGP);

% Riporto errore relativo con il metodo di gradiente precondizionato
true_relative_error_gp = norm(xm - xGP) / norm(xm);
disp("Errore relativo con il metodo di gradiente precondizionato:")
disp(true_relative_error_gp);

%% Punto 7
disp("Punto 7")
% Si risolva il sistema lineare (1) applicando il metodo del gradiente
% coniugato (non precon-
% dizionato) con medesimi valori di x0, toll e nmax.
% Si utilizzi la funzione gc.m utilizzata nel
% Lab 4 e si memorizzi la soluzione calcolata nel vettore xGC.
% Si riporti il numero di iterazioni
% kGC eﬀettuate per raggiungere la precisione prescritta e
% l'errore relativo commesso

% Calcolo la soluzione con il metodo del gradiente coniugato, mi aspetto
% che finisca in un numero di iterazioni inferiore
[xGC, kGC, res_rel_gc] = gc(A, b, x0, nmax, toll);

% Riporto il numero di iterazioni del metodo di gradiente coniugato
disp("Numero di iterazioni del metodo di gradiente coniugato:")
disp(kGC);

% Riporto errore relativo con il metodo del gradiente coniugato
true_relative_error_gc = norm(xm - xGC) / norm(xm);
disp("Errore relativo con il metodo di gradiente coniugato:")
disp(true_relative_error_gc);


%% Punto 8
% Disegnare su un grafico l'andamento del residuo normalizzato ∥r(k)∥
% ∥b∥ in funzione delle iter-
% azioni knei tre casi (gradiente, gradiente precondizionato e gradiente coniugato),
% e confrontare
% le curve ottenute. Caricare l'immagine ottenuta in formato .png.

% Ho salvato tutti e tre i vettori dell'andamento del residuo relativo
res_rel_gc; %Gradiente coniugato
res_rel_gpd;%Gradiente precondizionatp
res_rel_gdy; %Gradiente normale

% Plotto l'andamento del residuo relativo in scala semi-log sulle y
% e le iterazioni sulle x


figure
semilogy(0:kG, res_rel_gdy, 'b-', 0:kGP, res_rel_gpd, 'g-', 0:kGC, res_rel_gc, 'r-', 'LineWidth', 1.2)
grid on
xlabel('Iterazioni')
ylabel('Residuo relativo')
legend('Gradiente dinamico', 'Gradiente precondizionato', 'Gradiente coniugato')

%% Punto 9
% Calcolare i numeri di condizionamento delle matrici A e P−1A e commentare i precedenti
% risultati ottenuti alla luce della teoria.

%Calcolo numero di precondizionamento di A
kA = cond(A);
%Calcolo condizionamento della matrice precondizionata p^-1 * A
kA_prec = cond(P \ A);
format short
%Stampo i valori
disp("Numero di condizionamento della matrice A:")
disp(kA);
disp("Numero di condizionamento della matrice P^-1 * A:")
disp(kA_prec);

% Il precondizionatore P = diag(A) riduce il condizionamento da K(A) ~ 28.4
% a K(P^-1 A) ~ 14.8. Questo spiega perché il gradiente precondizionato
% converge in 160 iterazioni vs 348 del gradiente dinamico: il bound sulla
% velocità di convergenza dipende da (K-1)/(K+1), che decresce al diminuire
% di K. Più è piccolo più il bound stringe in fretta al crescere delle
% iterazioni.