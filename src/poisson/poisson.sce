
// Retourne la fréquence d'échantillonage de la transformée de Fourier discrète
function [freqs]=fftfreq(N, L)
    // Calculer les fréquences d'échantillonage en fonction de L et de la parité de N
if modulo(N, 2) == 0 then
    freqs = linspace(0, N / 2 - 1, N / 2);
    freqs = 2 * %i * %pi / L * cat(2, freqs, linspace(-N / 2, -1, N / 2));
else
    freqs = linspace(0, (N - 1) / 2, (N + 1) / 2);
    freqs = 2 * %i * %pi / L * cat(2, freqs, linspace(-(N - 1) / 2, -1, (N - 1) / 2));
end
endfunction


// Résolution de l'équation de Poisson en dimension 2 en utilisant la FFT
//    laplacien(psi) = f
// Entrée: f de taille (Ny,Nx) sur une domaine de taille (Ly,Lx)
// Sortie: psi, solution de l'équation
function [psi]=poisson_2d(f, Nx, Ny, Lx, Ly)
    kx = fftfreq(Nx, Lx);
    ky = fftfreq(Ny, Ly);
    psi_hat = zeros(Ny,Nx);
    f_hat = fft(f, "nonsymmetric");
    kx2 = kx^2;     // On calcule une fois pour toutes kx^2 et ky^2 pour eviter
    ky2 = ky^2;     // de les calculer a chaque passage de boucle.
    for p=1:Ny
        for q=1:Nx
            if p == 1 & q == 1 then
                psi_hat(p,q) = 0;
            else
                psi_hat(p,q) = f_hat(p, q) / (kx2(q) + ky2(p));
            end
        end
    end
    psi = real(ifft(psi_hat, "nonsymmetric"))
endfunction


// Résolution de l'équation de Poisson avec rot en dimension 2 en utilisant la FFT
//    laplacien(Ux) = -dW/dy
//    laplacien(Uy) = +dW/dx
// Entrée: champs de vorticité W de taille (Ny,Nx) sur un domaine de taille (Ly,Lx)
// Sortie: Ux et Uy, vitesses solution des équations
function [Ux,Uy]=poisson_curl_2d(W, Nx, Ny, Lx, Ly)
  kx = fftfreq(Nx,Lx)
  ky = fftfreq(Ny,Ly)
  ux_hat = zeros(Ny,Nx)
  uy_hat = zeros(Ny,Nx)
  W_hat = fft(W)
  for p=1:Ny
    for q=1:Nx
      if p==1 & q==1 then
	ux_hat(p,q) = 0
	uy_hat(p,q) = 0
      else
	ux_hat(p,q) = -ky(p) *( W_hat(p,q)/(kx(q)^2+ky(p)^2))
	uy_hat(p,q) = kx(q)*(W_hat(p,q)/(kx(q)^2 + ky(p)^2))
      end
    end
  end
    Ux = real(ifft(ux_hat, "nonsymmetric"))
    Uy = real(ifft(uy_hat, "nonsymmetric"))
endfunction

