function [ pri_si ] = calcola_pr_s_synch( rj, si, si_hat, hj, eps, Pmal, flag_byz_behavior )
m = length(rj);
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
if hj == 1
    %Honest node
    deltas = zeros(m,1);
    deltas(rj == si) = 1;
    pri_si = (1-eps)*deltas+eps*(~deltas);
else
    if flag_byz_behavior == 1
        if si_hat == 0
            eta = eps*(1-Pmal(1))+(1-eps)*Pmal(1);
        else
            eta = eps*(1-Pmal(2))+(1-eps)*Pmal(2);
        end;
        deltas = zeros(m,1);
        deltas(rj == si) = 1;
        pri_si = (1-eta)*deltas+eta*(~deltas);
    else
        if si_hat == 0
            eps = eps0;
            deltas = zeros(m,1);
            deltas(rj == si_hat) = 1;
            pri_si = (1-eps)*deltas+eps*(~deltas);
        else
            eps = eps1;
            deltas = zeros(m,1);
            deltas(rj == si_hat) = 1;
            pri_si = (1-eps)*deltas+eps*(~deltas);
        end;
    end;
end;


