
function dt=calcul_dt(U,dx)
  dt=dx/max(abs(U));
endfunction

function vort=solveur_1D(vort, Ux, Nx, kappa, dt, dx)
    //----------------------------------------------
    //implémente une itération temporelle
    // de l'algo 1D codé dans le script dif-conv.sce
    //   vort = champ transporté (correspond à phi dans le sujet)
    //   Ux = vitesse sur la composante x
    //   Nx = nombre de points de discrétisation spatiale en x et y.
    //   kappa = constante de diffusion dynamique
    //   dt = pas de temps
    //   dx = pas d'espace
    //-----------------------------------------------
    J = zeros(Nx,Nx) //on commence par assembler N
    J(1,$) = -1
    J($,1) = -1
    J = J + 2*diag(ones(Nx,1)) - diag(ones(Nx-1,1),-1) - diag(ones(Nx-1, 1),1)
    N = diag(ones(Nx,1)) + kappa * (dt / dx^2) * J

    A = diag(ones(Nx-1,1),1) - diag(ones(Nx-1,1),-1)
    A(1,$)=-1
    A($,1)=1
    for i = 1:Nx
      A(i,:) = A(i,:) * Ux(i) * (dt/dx) * 0.5
    end
    B = diag(ones(Nx-1,1),1) - 2 * diag(ones(Nx,1)) + diag(ones(Nx-1,1),-1)
    B(1,$)=1
    B($,1)=1
    for i = 1:Nx
      B(i,:) = B(i,:) * Ux(i)^2 * (dt/dx)^2 * 0.5
    end
    M = diag(ones(Nx,1)) - A + B
    vort  = umfpack(sparse(N), "\",M*vort)
endfunction

function vort=solveur_2D(vort, Ux, Uy, Nx, Ny, kappa, dt, dx, dy)
    //----------------------------------------------
    //implémente une itération temporelle
    // de l'algo de splitting 2D utilisant l'algo 1D
    //   vort = champ 2D transporté (correspond à phi dans le sujet)
    //   Ux = vitesse sur la composante x
    //   Uy = vitesse sur la composante y
    //   (Nx, Ny) nombre de points de discrétisation spatiale en x et y.
    //   kappa = constante de diffusion dynamique
    //   dt    = pas de temps
    //   (dx, dy) pas d'espaces dans chaque direction
    //-----------------------------------------------
    for i = 1:Ny
      vort(i,:) = (solveur_1D(vort(i,:)',Ux(i,:),Nx,kappa,dt,dx))'
    end
    for j = 1:Nx
      vort(:,j) = solveur_1D(vort(:,j),Uy(:,j)',Ny,kappa,dt,dy)
    end
endfunction
