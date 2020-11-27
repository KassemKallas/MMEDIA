function [ phi_r ] = calc_phi_r( rho, tau_r, pS10 )
if length(rho) == 1
    rho = [rho rho];
end;
phi_r = 0*tau_r;
phi_r(1) = pS10;
phi_r(2:end) = rho(1)*tau_r(1:end-1)+(1-rho(2))*(1-tau_r(1:end-1));
end
