function dTdt = mixMe(t,T)

global R0 Ve C m0 C_p0 C_a rho_a

% Same rho and Cp
%     dTdt = - (Ve.^3*t.^2)./((R0^3)/3+ Ve.^3*t.^3);

% Better
    dTdt = C * T.*t.^2 / (m0*C_p0 + 4*pi*C_a*rho_a*Ve^3.*t.^3);