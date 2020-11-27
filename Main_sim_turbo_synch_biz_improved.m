rng('default');
rng(12);
Niteraz = 1000;

flag_dec_model = 1;%1 = new approach (synch) - 0 = old approach (no synch)
rho = 0.5;
rho_hat = 0.5;
rho_hat_dec = rho_hat;
flag_byz_behavior = 0;%1 = classic approach but syncronized (use 2 different Pmal) - 0 = zero mutual information (reports depend on s_hat only)
Pmal = [0.5 0.5];
eps = 0.1;
alpha = 0.45;
s1 = -1;%-1 = unknown
s1_hat = -1;%-1 = unknown

m = 20; %numero osservazioni
n = 20; %numero nodi
Num_iteraz_turbo = 5;
ER = zeros(Num_iteraz_turbo,Niteraz);
ER_OTT_bound = zeros(1,Niteraz);
P_byz_byz = zeros(Num_iteraz_turbo,Niteraz);
P_byz_honest = zeros(Num_iteraz_turbo,Niteraz);
Nbyz = round(alpha*n);

for IT = 1:Niteraz
    tau_l = 0.5*ones(m,1);
    tau_r = 0.5*ones(m,1);
    tau_l_hat = 0.5*ones(m,1);
    tau_r_hat = 0.5*ones(m,1);
    
    phi_l = 0.5*ones(m,1);
    phi_r = 0.5*ones(m,1);
    phi_l_hat = 0.5*ones(m,1);
    phi_r_hat = 0.5*ones(m,1);
    
    nu_u = 0.5*ones(m,n);
    nu_d = 0.5*ones(m,n);
    nu_u_hat = 0.5*ones(m,n);
    nu_d_hat = 0.5*ones(m,n);
    
    lambda_d = 0.5*ones(n,m);
    lambda_u = 0.5*ones(n,m);
    omega_d = 0.5*ones(n,1);
    omega_u = 0.5*ones(n,1);
    
    % Previous one
    
    
    MESS = zeros(2,m);
    MESS_hat = zeros(2,m);
    if rem(IT,50) == 0
        fprintf('Iterazione %d su %d\n',IT,Niteraz);
    end;
    
    [ R, s, s_hat ] = crea_reports_Markov_synch_biz( m, n, alpha, Pmal, eps, rho, rho_hat, s1, s1_hat, flag_byz_behavior);
    
    %R = ALL_GEN(IT).R;
    %s = ALL_GEN(IT).s;
    
    %[ R, s, Nbyz ] = crea_reports_Markov_stat( m, n, alpha, Pmal, eps, rho, s1);
    % Iterazioni turbo
    % Inizializzazione omega
    omega_u = alpha;
    EST = zeros(1,m);
    EST_hat = zeros(1,m);
    for tit = 1:Num_iteraz_turbo
        if flag_dec_model == 1
            [ lambda_u ] = calc_lambda_u( omega_u, lambda_d, -1 );
            
            [ nu_u ] = calc_nu_u_synch( R, lambda_u, nu_d_hat, eps, Pmal, flag_byz_behavior, 1 );
            [ tau_r ] = calc_tau_r( phi_r, nu_u );
            [ tau_l ] = calc_tau_l( phi_l, nu_u );
            [ phi_r ] = calc_phi_r( rho, tau_r, 0.5 );
            [ phi_l ] = calc_phi_l( rho, tau_l );
            [ nu_d ] = calc_nu_d( phi_r, phi_l, nu_u );
            
            [ nu_u_hat ] = calc_nu_u_synch( R, lambda_u, nu_d, eps, Pmal, flag_byz_behavior, 0 );
            [ tau_r_hat ] = calc_tau_r( phi_r_hat, nu_u_hat );
            [ tau_l_hat ] = calc_tau_l( phi_l_hat, nu_u_hat );
            [ phi_r_hat ] = calc_phi_r( rho_hat_dec, tau_r_hat, 0.5 );
            [ phi_l_hat ] = calc_phi_l( rho_hat_dec, tau_l_hat );
            [ nu_d_hat ] = calc_nu_d( phi_r_hat, phi_l_hat, nu_u_hat );
            
            [ lambda_d ] = calc_lambda_d_synch( R, nu_d, nu_d_hat, eps, Pmal, flag_byz_behavior );
            [ omega_d ] = calc_omega_d( lambda_d );
            
            MESS(1,:) = exp(sum(log(nu_u+1e-20),2)+log(phi_r+1e-20)+log(phi_l+1e-20));
            MESS(2,:) = exp(sum(log(1-nu_u+1e-20),2)+log(1-phi_r+1e-20)+log(1-phi_l+1e-20));
            MESS_hat(1,:) = exp(sum(log(nu_u_hat+1e-20),2)+log(phi_r_hat+1e-20)+log(phi_l_hat+1e-20));
            MESS_hat(2,:) = exp(sum(log(1-nu_u_hat+1e-20),2)+log(1-phi_r_hat+1e-20)+log(1-phi_l_hat+1e-20));
            PS = MESS./repmat(sum(MESS),2,1);
            PS_hat = MESS_hat./repmat(sum(MESS_hat),2,1);
            MESS1_byz = exp(sum(log(lambda_d+1e-20),2)+log(omega_u));
            MESS2_byz = exp(sum(log(1-lambda_d+1e-20),2)+log(1-omega_u));
            P_byz = MESS1_byz./(MESS1_byz+MESS2_byz);
            
            PS_ALL(tit,IT).vals = PS;
            PS_ALL_hat(tit,IT).vals = PS_hat;
            P_byz_ALL(tit,IT).vals = P_byz;
            %PS_ALL2(tit,IT).val = PS;
            [maxv indxm] = max(MESS);
            EST = indxm-1;
            indx = find(abs(PS(1,:)-0.5) < 1e-8);
            EST(indx) = 0;
            [maxv indxm] = max(MESS_hat);
            EST_hat = indxm-1;
            indx = find(abs(PS_hat(1,:)-0.5) < 1e-8);
            EST_hat(indx) = 0;
            
            ER(tit,IT) = sum(xor(EST,s))/m;
            P_byz_byz(tit,IT) = mean(P_byz(1:Nbyz));
            P_byz_honest(tit,IT) = mean(P_byz(Nbyz+1:end));
            
        else
            [ lambda_u ] = calc_lambda_u( omega_u, lambda_d, -1 );
            [ nu_u ] = calc_nu_u( lambda_u, R, eps, Pmal(1), -1 );
            [ tau_r ] = calc_tau_r( phi_r, nu_u );
            [ tau_l ] = calc_tau_l( phi_l, nu_u );
            [ phi_r ] = calc_phi_r( rho, tau_r, 0.5 );
            [ phi_l ] = calc_phi_l( rho, tau_l );
            [ nu_d ] = calc_nu_d( phi_r, phi_l, nu_u );
            [ lambda_d ] = calc_lambda_d( nu_d, R, eps, Pmal(1), -1 );
            [ omega_d ] = calc_omega_d( lambda_d );
            
            MESS(1,:) = exp(sum(log(nu_u+1e-20),2)+log(phi_r+1e-20)+log(phi_l+1e-20));
            MESS(2,:) = exp(sum(log(1-nu_u+1e-20),2)+log(1-phi_r+1e-20)+log(1-phi_l+1e-20));
            MESS1_byz = exp(sum(log(lambda_d+1e-20),2)+log(omega_u));
            MESS2_byz = exp(sum(log(1-lambda_d+1e-20),2)+log(1-omega_u));
            
            PS = MESS./repmat(sum(MESS),2,1);
            P_byz = MESS1_byz./(MESS1_byz+MESS2_byz);
            
            PS_ALL(tit,IT).vals = PS;
            P_byz_ALL(tit,IT).vals = P_byz;
            
            %PS_ALL2(tit,IT).val = PS;
            [maxv indxm] = max(MESS);
            EST = indxm-1;
            indx = find(abs(PS(1,:)-0.5) < 1e-8);
            EST(indx) = 0;
            
            ER(tit,IT) = sum(xor(EST,s))/m;
            P_byz_byz(tit,IT) = mean(P_byz(1:Nbyz));
            P_byz_honest(tit,IT) = mean(P_byz(Nbyz+1:end));
            
        end;
    end;
    if ER(Num_iteraz_turbo,IT) > 0
        check = 1;
    end;
    %Verifica se la sequenza alternativa (EST2) ha maggiore probabilita'
    
    if flag_byz_behavior == 0 & Pmal(1) == 1 & Pmal(2) == 1
        EST1 = EST;
        EST2 = ~EST1;
    else
        EST1 = EST;
        EST2 = EST_hat;
    end;
    indx_byz_1 = find(P_byz > 0.5);
    indx_hon_1 = setxor(1:length(P_byz),indx_byz_1);
    indx_byz_2 = find(P_byz <= 0.5);
    indx_hon_2 = setxor(1:length(P_byz),indx_byz_1);
    
    %Calcolo della probabilita' a priori
    diff_s = xor(EST1(2:end),EST1(1:end-1));
    pS1 = prod((1-rho).^(~diff_s).*rho.^(diff_s));
    
    %Calcolo della probabilita' a priori
    diff_s = xor(EST2(2:end),EST2(1:end-1));
    pS2 = prod((1-rho).^(~diff_s).*rho.^(diff_s));
    
    Pr_1_log = log(pS1);
    Pr_2_log = log(pS2);

    for q = 1:m
        [ pri_si ] = calcola_pr_s_synch( R(:,q), EST1(q), EST2(q), 0, eps, Pmal, flag_byz_behavior );
        Pr_eff_log = Pr_eff_log + sum(log(pri_si(1:Nbyz)));
        [ pri_si ] = calcola_pr_s_synch( R(:,q), s(q), s_hat(q), 1, eps, Pmal, flag_byz_behavior );
        Pr_eff_log = Pr_eff_log + sum(log(pri_si(Nbyz+1:end)));
        [ pri_si ] = calcola_pr_s_synch( R(:,q), EST(q), EST_hat(q), 0, eps, Pmal, flag_byz_behavior );
        Pr_dec_log = Pr_dec_log + sum(log(pri_si(1:Nbyz)));
        [ pri_si ] = calcola_pr_s_synch( R(:,q), EST(q), EST_hat(q), 1, eps, Pmal, flag_byz_behavior );
        Pr_dec_log = Pr_dec_log + sum(log(pri_si(Nbyz+1:end)));
    end;
        
        
        delta = eps*(1-Pmal(1))+(1-eps)*Pmal(1);
    
    MAT_DIFF1 = xor(R,repmat(EST1,n,1));
    metrica1 = log(pS1) + sum(log((1-alpha)*prod(((1-eps).^(~MAT_DIFF1)).*(eps.^(MAT_DIFF1)),2)+alpha*prod(((1-delta).^(~MAT_DIFF1)).*(delta.^(MAT_DIFF1)),2)));
    
    MAT_DIFF2 = xor(R,repmat(EST2,n,1));
    metrica2 = log(pS2) + sum(log((1-alpha)*prod(((1-eps).^(~MAT_DIFF2)).*(eps.^(MAT_DIFF2)),2)+alpha*prod(((1-delta).^(~MAT_DIFF2)).*(delta.^(MAT_DIFF2)),2)));
    if metrica1 > metrica2
        EST = EST1;
    else
        P_byz_byz(tit,IT) = 1-mean(P_byz(1:Nbyz));
        P_byz_honest(tit,IT) = 1-mean(P_byz(Nbyz+1:end));
        EST = EST2;
    end;
    ER(Num_iteraz_turbo,IT) = sum(xor(EST,s))/m;
    
    if ER(Num_iteraz_turbo,IT) > 0
        check = 1;
    end;
    
    %Calcolo del bound sulla decisione ottima
    Pr_eff_log = 0;
    Pr_dec_log = 0;
    
    %Calcolo della probabilita' a priori
    EST1 = s;
    indx0 = find(EST1(1:end-1) == 0);
    indx1 = find(EST1(1:end-1) == 1);
    diff_s = xor(EST1(2:end),EST1(1:end-1));
    pS1 = 1;
    if ~isempty(indx0)
        pS1 = prod((1-rho).^(diff_s(indx0)).*rho.^(~diff_s(indx0)));
    end;
    if ~isempty(indx1)
        pS1 = pS1*prod((1-rho).^(diff_s(indx1)).*rho.^(~diff_s(indx1)));
    end;
    Pr_eff_log = Pr_eff_log + log(pS1);
    
    EST1 = s_hat;
    indx0 = find(EST1(1:end-1) == 0);
    indx1 = find(EST1(1:end-1) == 1);
    diff_s = xor(EST1(2:end),EST1(1:end-1));
    pS1 = 1;
    if ~isempty(indx0)
        pS1 = prod((1-rho_hat).^(diff_s(indx0)).*rho_hat.^(~diff_s(indx0)));
    end;
    if ~isempty(indx1)
        pS1 = pS1*prod((1-rho_hat).^(diff_s(indx1)).*rho_hat.^(~diff_s(indx1)));
    end;
    Pr_eff_log = Pr_eff_log + log(pS1);
    
    EST1 = EST;
    indx0 = find(EST1(1:end-1) == 0);
    indx1 = find(EST1(1:end-1) == 1);
    diff_s = xor(EST1(2:end),EST1(1:end-1));
    pS1 = 1;
    if ~isempty(indx0)
        pS1 = prod((1-rho).^(diff_s(indx0)).*rho.^(~diff_s(indx0)));
    end;
    if ~isempty(indx1)
        pS1 = pS1*prod((1-rho).^(diff_s(indx1)).*rho.^(~diff_s(indx1)));
    end;
    Pr_dec_log = Pr_dec_log + log(pS1);
    
    EST1 = EST_hat;
    indx0 = find(EST1(1:end-1) == 0);
    indx1 = find(EST1(1:end-1) == 1);
    diff_s = xor(EST1(2:end),EST1(1:end-1));
    pS1 = 1;
    if ~isempty(indx0)
        pS1 = prod((1-rho_hat).^(diff_s(indx0)).*rho_hat.^(~diff_s(indx0)));
    end;
    if ~isempty(indx1)
        pS1 = pS1*prod((1-rho_hat).^(diff_s(indx1)).*rho_hat.^(~diff_s(indx1)));
    end;
    Pr_dec_log = Pr_dec_log + log(pS1);
    
    %Calcolo delle probabilita' a posteriori
    if ER(Num_iteraz_turbo,IT) > 0
        for q = 1:m
            [ pri_si ] = calcola_pr_s_synch( R(:,q), s(q), s_hat(q), 0, eps, Pmal, flag_byz_behavior );
            Pr_eff_log = Pr_eff_log + sum(log(pri_si(1:Nbyz)));
            [ pri_si ] = calcola_pr_s_synch( R(:,q), s(q), s_hat(q), 1, eps, Pmal, flag_byz_behavior );
            Pr_eff_log = Pr_eff_log + sum(log(pri_si(Nbyz+1:end)));
            [ pri_si ] = calcola_pr_s_synch( R(:,q), EST(q), EST_hat(q), 0, eps, Pmal, flag_byz_behavior );
            Pr_dec_log = Pr_dec_log + sum(log(pri_si(1:Nbyz)));
            [ pri_si ] = calcola_pr_s_synch( R(:,q), EST(q), EST_hat(q), 1, eps, Pmal, flag_byz_behavior );
            Pr_dec_log = Pr_dec_log + sum(log(pri_si(Nbyz+1:end)));
        end;
        if Pr_dec_log >= Pr_eff_log
            ER_OTT_bound(IT) = ER(Num_iteraz_turbo,IT);
        end;
    end;
    %Calcolo con decisione a maggioranza
    MAJ_C = sum(R);
    EST_MAJ = zeros(1,m);
    EST_MAJ(MAJ_C > n/2) = 1;
    ER_MAJ(IT) = sum(xor(EST_MAJ,s))/m;
    
    %Calcolo con decisione a maggioranza con rimozione
    REP = zeros(n,1);
    for g = 1:n
        REP(g) = sum(xor(R(g,:),EST_MAJ));
    end;
    [sortv indxs] = sort(REP);
    indx_honest = indxs(1:n-Nbyz);
    %indx_honest = Nbyz+1:m;
    Rh = R(indx_honest,:);
    MAJ_C = sum(Rh);
    EST_MAJ = zeros(1,m);
    EST_MAJ(MAJ_C > length(indx_honest)/2) = 1;
    ER_MAJ_r(IT) = sum(xor(EST_MAJ,s))/m;
    if rem(IT,100) == 0
        fprintf('Error rate con Turbo decoding = %f\n',mean(ER(Num_iteraz_turbo,1:IT)));
        fprintf('Lower bound error rate con optimal decoding = %f\n',mean(ER_OTT_bound(1:IT)));
    end;
end;
fprintf('Errr rate with majority rule = %f\n',mean(ER_MAJ));
fprintf('Errr rate with majority rule after removals = %f\n',mean(ER_MAJ_r));
for l = 1:Num_iteraz_turbo
    fprintf('Error rate with Message Passing after %d iteration = %f - Lower bound optimal = %f\n',l,mean(ER(l,:)),mean(ER_OTT_bound));
    fprintf('Probability of byzantine conditioned on byzantine after %d iteration = %f\n',l,mean(P_byz_byz(l,:)));
    fprintf('Probability of byzantine conditioned on honest after %d iteration = %f\n',l,mean(P_byz_honest(l,:)));
end;

%fprintf('Errr rate con algoritmo BCJR prima della rimozione = %f\n',mean(ER));
%fprintf('Errr rate con algoritmo BCJR dopo la rimozione = %f\n',mean(ER_r));


