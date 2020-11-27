function [ R, s, s_hat ] = crea_reports_Markov_synch_biz( n, m, alpha, Pmal, eps, rho, rho_hat, s1, s1_hat, flag_byz_behavior)

%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
% n = numero di osservazioni temporali
% m = numero di nodi
% alpha = percentuale di bizantini
% Pmal = Probabilita' di flipping
% eps = errore di misura
% rho = probabilita' del modello
% s1 = stato uniziale
rho = 1-rho;
rho_hat = 1-rho_hat;
if length(rho) == 1
    rho = [rho rho];
end;
if length(rho_hat) == 1
    rho_hat = [rho_hat rho_hat];
end;
if length(eps) == 1
    eps = [eps eps];
end;
eps0 = eps(1);
eps1 = eps(2);
R = zeros(m,n);
s = zeros(1,n);
s_hat = zeros(1,n);

if s1 == -1
    s1 = randi(2)-1;
end
if s1_hat == -1
    s1_hat = randi(2)-1;
end
s(1) = s1;
s_hat(1) = s1_hat;

for k = 2:n
    cs = rand;
    if ~s(k-1)
        rho_e = rho(1);
    else
        rho_e = rho(2);
    end;
    if cs < rho_e
        s(k) = ~s(k-1);
    else
        s(k) = s(k-1);
    end;  
    cs = rand;
    if ~s_hat(k-1)
        rho_e = rho_hat(1);
    else
        rho_e = rho_hat(2);
    end;
    if cs < rho_e
        s_hat(k) = ~s_hat(k-1);
    else
        s_hat(k) = s_hat(k-1);
    end;      
end
Noise_dec = rand(m,n);
indx0 = find(s == 0);
indx1 = find(s == 1);
U = repmat(s,m,1);
Noise_dec_d0 = Noise_dec;
Noise_dec_d0(:,indx1) = 1;
Noise_dec_d1 = Noise_dec;
Noise_dec_d1(:,indx0) = 1;
if length(indx0) > 0
    indx_flip0 = find( Noise_dec_d0 < eps0 );    
    U(indx_flip0) = ~U(indx_flip0);
end;
if length(indx1) > 0
    indx_flip1 = find( Noise_dec_d1 < eps1 );
    U(indx_flip1) = ~U(indx_flip1);
end;
R = U;
%Per comodita' i bizantini sono i primi

Num_B = round(alpha*m);
for biz_indx = 1:Num_B
    for state = 1:n
        if flag_byz_behavior == 1
            Noise_flip = rand;
            if s_hat(state) == 0
                if Noise_flip < Pmal(1)
                    R(biz_indx,state) = ~R(biz_indx,state);
                end;
            else
                if Noise_flip < Pmal(2)
                    R(biz_indx,state) = ~R(biz_indx,state);
                end;
            end;
        else
            Noise_dec = rand;
            if s_hat(state) == 0
                if Noise_dec < eps0
                    R(biz_indx,state) = ~s_hat(state);
                else
                    R(biz_indx,state) = s_hat(state);
                end;
            else
                if Noise_dec < eps1
                    R(biz_indx,state) = ~s_hat(state);
                else
                    R(biz_indx,state) = s_hat(state);
                end;
            end;
        end;
    end;
end;
