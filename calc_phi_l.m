function [ phi_l ] = calc_phi_l( rho, tau_l )
if length(rho) == 1
    rho = [rho rho];
end;
phi_l = 0*tau_l;
phi_l(end) = rho(1)/(rho(1)+rho(2));
phi_l(1:end-1) = rho(1)*tau_l(2:end)+(1-rho(2))*(1-tau_l(2:end));
end
