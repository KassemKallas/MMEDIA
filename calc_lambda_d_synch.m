function [ lambda_d ] = calc_lambda_d_synch( R, nu_d, nu_d_hat, eps, Pmal, flag_byz_behavior )

m = size(nu_d,1); %Numero osservazioni
n = size(nu_d,2); %Numero nodi
lambda_d = zeros(n,m);
indx_ch = 1:m;

for i = indx_ch
    NUM = 0;
    DEN = 0;
    for si = 0:1
        for si_hat = 0:1
            T0 = (calcola_pr_s_synch( R(:,i), si, si_hat, 0, eps, Pmal, flag_byz_behavior )).';
            T0 = T0.*(nu_d(i,:)*(1-si)+si*(1-nu_d(i,:))).*(nu_d_hat(i,:)*(1-si_hat)+(1-nu_d_hat(i,:))*si_hat);
            NUM = NUM + T0;
            DEN = DEN + T0;
            T1 = (calcola_pr_s_synch( R(:,i), si, si_hat, 1, eps, Pmal, flag_byz_behavior )).';
            T1 = T1.*(nu_d(i,:)*(1-si)+si*(1-nu_d(i,:))).*(nu_d_hat(i,:)*(1-si_hat)+(1-nu_d_hat(i,:))*si_hat);
            DEN = DEN + T1;
        end;
        
    end;
    lambda_d(:,i) = NUM./DEN;
end

