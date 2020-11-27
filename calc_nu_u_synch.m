function [ nu_u ] = calc_nu_u_synch( R, lambda_u, nu_d_hat, eps, Pmal, flag_byz_behavior, flag_verso )
m = size(lambda_u,2); %Numero osservazioni
n = size(lambda_u,1); %Numero nodi
nu_u = zeros(m,n);
indx_ch = 1:m;
for i = indx_ch
    NUM = 0;
    DEN = 0;
    for hj = 0:1
        if flag_verso == 1
            for si_hat = 0:1
                T0 = calcola_pr_s_synch( R(:,i), 0, si_hat, hj, eps, Pmal, flag_byz_behavior );
                T0 = T0.*(lambda_u(:,i)*(1-hj)+hj*(1-lambda_u(:,i))).*(nu_d_hat(i,:).'*(1-si_hat)+(1-nu_d_hat(i,:).')*si_hat);
                NUM = NUM + T0;
                DEN = DEN + T0;
                T1 = calcola_pr_s_synch( R(:,i), 1, si_hat, hj, eps, Pmal, flag_byz_behavior );
                T1 = T1.*(lambda_u(:,i)*(1-hj)+hj*(1-lambda_u(:,i))).*(nu_d_hat(i,:).'*(1-si_hat)+(1-nu_d_hat(i,:).')*si_hat);
                DEN = DEN + T1;
            end;
        else
            for si = 0:1
                T0 = calcola_pr_s_synch( R(:,i), si, 0, hj, eps, Pmal, flag_byz_behavior );
                T0 = T0.*(lambda_u(:,i)*(1-hj)+hj*(1-lambda_u(:,i))).*(nu_d_hat(i,:).'*(1-si)+(1-nu_d_hat(i,:).')*si);
                NUM = NUM + T0;
                DEN = DEN + T0;
                T1 = calcola_pr_s_synch( R(:,i), si, 1, hj, eps, Pmal, flag_byz_behavior );
                T1 = T1.*(lambda_u(:,i)*(1-hj)+hj*(1-lambda_u(:,i))).*(nu_d_hat(i,:).'*(1-si)+(1-nu_d_hat(i,:).')*si);
                DEN = DEN + T1;
            end;
        end;
    end;
    nu_u(i,:) = NUM./DEN;
end

