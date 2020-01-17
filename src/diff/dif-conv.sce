Nx=200;
Nt=200;
kappa=1e-2;
exec("my_cholesky.sce")

function y=phi_0(x)
//--------------------
// condition initiale
//--------------------
if 0 <= x & x < 0.25 then
    y = 0;
elseif 0.25 <= x & x < 0.375 then
    y = 2 * (x - 0.25);
elseif 0.375 <= x & x < 0.5 then
    y = 2 * (0.5 - x);
elseif 0.5 <= x & x <= 1 then
    y = 0;
end
endfunction


function y=conv(x)
//--------------------
// fonction de convection
//--------------------
y = 0.4 * (x - 0.25);
endfunction


//--------------------
// initialiser phi et assembler les matrices M et N
//-------------------
fin=Nt;
dt = 1 / Nt;
dx = 1 / Nx;
maillage = linspace(0,1,Nx);
phi_i = [];
indice = 1;
for element=maillage,
    phi_i(indice) = phi_0(element);
    indice = indice + 1;
end

// Initialisation des matrices pour assembler M et N.
I = eye(Nx, Nx);
J = [zeros(1, Nx);eye(Nx - 1, Nx)] + [zeros(Nx, 1),eye(Nx, Nx - 1)];
J(1, $) = 1;
J($, 1) = 1;
K = [zeros(1, Nx);eye(Nx - 1, Nx)];
K(1, $) = 1;     // Prendre K' pour former l'autre partie de la matrice.

// Assemblage des matrices (M doit se modifier au fur et a mesure a cause de c(x)).
N = kappa * dt / dx^2 * J + (-1 - 2 * kappa * dt / dx^2) * I;

x = 0;
phi = phi_i;    // Initialisation de phi.
for i=1:fin
    // On recalcule M avec la nouvelle valeur de x.
    M = (1 - conv(x)^2 * dt^2 / dx^2) * I + conv(x) * dt / (2 * dx) * K + conv(x) * dt / (2 * dx) * (conv(x) * dt / dx - 1) * K';
    phi=my_cholesky(N,M*phi);
    x = x + dx;
end
scf;
fig = gcf();
a = gca();
title = a.title;
title.foreground = 9;
title.font_size = 4;
title.font_style = 5;
title.text = "$\kappa = 10^{-2}$"
plot(maillage, [phi_i phi]);
xs2png(fig, "kappa_1e-2.png")
