// définition des pas de maillage
dx=Lx/Nx;
dy=Ly/Ny;

// définition de kappa
kappa = 1e-4;

exec("dif-conv-f.sce")

// définition des maillages horizontaux et verticaux
maillage_x=linspace(0,(Nx-1)/Nx*Lx,Nx)';
maillage_y=linspace(0,(Ny-1)/Ny*Ly,Ny)';

// initialisation des vecteurs convections
cx=zeros(Nx,Ny); //composante x de la vitesse de convection
cy=zeros(Nx,Ny); //composante y de la vitesse de convection

// initialisation du vecteur initiale
phi_i=zeros(Nx,Ny); //condtion initiale

//------------------------------------------
//Remplissage des tableaux cx cy phi phi_i
//------------------------------------------

// Initialisation des tableaux cx et cy et phi_i
for i = 1 : Nx
    for j = 1 : Ny
        cx(i,j) = conv_2d((i-1)*dy, (j-1)*dx)(1);
        cy(i,j) = conv_2d((i-1)*dy, (j-1)*dx)(2);
        phi_i(i,j) = phi_0_2d((i-1)*dx, (j-1)*dy);
     end
end

// Initialisation de phi par phi_i
phi = phi_i;

// calcul du pas de temps
dt=min(calcul_dt(cx,dx),calcul_dt(cy,dy));
Nt=floor(Tf/dt);

// résolution de l'équation
for k=1:Nt
    phi = solveur_2D(phi, cx, cy, Nx, Ny, kappa, dt, dx, dy);
end
