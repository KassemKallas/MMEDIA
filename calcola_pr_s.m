function [ pri_si ] = calcola_pr_s( rj, si, eps, Pmal, alpha )
if length(eps) == 1
    eps = [eps eps];
end;
eps0 = eps(1);
eps1 = eps(2);
if si == 0
    eps = eps0;
else
    eps = eps1;
end;
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
% rj = vettore dei report 
% si = stato al tempo i
% eps = errore di misura
% Pmal = Probabilita' di flipping
% alpha = percentuale di bizantini

m = length(rj);
eta = eps*(1-Pmal)+(1-eps)*Pmal;
deltas = zeros(m,1);
deltas(rj == si) = 1;
pri_si = ((1-eps)*deltas+eps*(~deltas))*(1-alpha)+((1-eta)*deltas+eta*(~deltas))*alpha;
end

