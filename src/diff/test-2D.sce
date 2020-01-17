// Variables utiles à la résolution
Nx=100;
Ny=100;
nu=0.001
Lx=1;
Ly=1;
Tf=0.5;

function u=conv_2d(y,x)
    // cette fonction permet de calculer le vecteur de convection
   alpha=1;
   beta=1;
   u = beta.*[cos(alpha)*x-sin(alpha)*y,sin(alpha)*x+cos(alpha)*y];
endfunction

function z=phi_0_2d(y,x)
  // cette fonction permet de calculer la fonction initiale
  p_0=[0.5 0.3];
  r_0=0.2;
 if (x-p_0(1))**2+(y-p_0(2))**2>r_0**2 then
    z=0;
  else
    z=1-((x-p_0(1))**2+(y-p_0(2))**2)/r_0**2;
  end
endfunction

exec("dif-conv-2D.sce")

//---------------------
//Affichage graphique
//--------------------
// on commence par couper la fenêtre graphique en deux parties pour chacun
// des deux graphiques
fig = gcf()
xset("colormap", jetcolormap(64))

subplot(121)
// on trace la solution initiale en 3D
colorbar(min(phi_i), max(phi_i))
a = gca()
a.isoview = "on"
title = a.title;
title.foreground = 9;
title.font_size = 4;
title.font_style = 5;
title.text = "Condition initiale"
plot3d1(maillage_x, maillage_y, phi_i)

subplot(122)
// on trace la solution finale en 3D
colorbar(min(phi), max(phi))
a = gca()
a.isoview = "on"
title = a.title;
title.foreground = 9;
title.font_size = 4;
title.font_style = 5;
title.text = "Solution au temps final"
plot3d1(maillage_x, maillage_y, phi)

xs2png(fig, "poisson_2d.png")
