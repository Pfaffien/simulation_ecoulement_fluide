// Initialisation de oméga, notée w par la suite
function [w]=omega(y,x)
    w = 8*(%pi^2)*cos(2*%pi*x)*cos(2*%pi*y)
endfunction


// Initialise -dw/dy
function [f]=init_field_1(y,x)
    f = 16*(%pi^3)*cos(2*%pi*x)*sin(2*%pi*y)
endfunction


// Initialise dw/dx
function [f]=init_field_2(y,x)
    f = -16*(%pi^3)*sin(2*%pi*x)*cos(2*%pi*y)
    endfunction


// Solution de référence du problème laplacien(u_x(x,y)) = -dw/dy
function [ref]=solution_field_1(y,x)
    ref = -2*%pi*cos(2*%pi*x)*sin(2*%pi*y)
endfunction


// Solution de référence du problème laplacien(u_x(x,y)) = dw/dx
function [ref]=solution_field_2(y,x)
    ref = 2*%pi*sin(2*%pi*x)*cos(2*%pi*y)
endfunction


// Affichage de oméga, des solutions de référence selon x et y,
// des solutions obtenues selon x et y, ainsi que l'erreur commise
// selon x et y.
// w     = oméga(x,y)   -- oméga
// ref_x = u_x(x,y)     -- solution analytique selon x
// ref_y = u_y(x,y)     -- solution analytique selon y
// sol_x                -- solution du solveur selon x
// sol_y                -- solution du solveur selon y

function plot_error(w, ref_x, ref_y, sol_x, sol_y)
    fig = gcf()
    xset("colormap", jetcolormap(64))

    subplot(331)
    colorbar(min(ref_x),max(ref_x))
    a = get("current_axes");
    title = a.title;
    title.foreground = 9;
    title.font_size = 4;
    title.text = "$\textbf{u}_{x,ref}(x,y)$"
    plot3d1(Y,X,ref_x)

    subplot(332)
    colorbar(min(sol_x),max(sol_x))
    a = get("current_axes");
    title = a.title;
    title.foreground = 9;
    title.font_size = 4;
    title.text = "$\textbf{u}_{x,solveur}(x,y)$"
    plot3d1(Y,X,sol_x)

    subplot(333)
    colorbar(min(abs(ref_x - sol_x)),max(abs(ref_x - sol_x)))
    a = get("current_axes");
    title = a.title;
    title.foreground = 9;
    title.font_size = 4;
    title.text = "$|\textbf{u}_{x,ref} - \textbf{u}_{x,solveur}|(x,y)$"
    plot3d1(Y,X,abs(ref_x - sol_x))

    subplot(334)
    colorbar(min(ref_y),max(ref_y))
    a = get("current_axes");
    title = a.title;
    title.foreground = 9;
    title.font_size = 4;
    title.text = "$\textbf{u}_{y,ref}(x,y)$"
    plot3d1(Y,X,ref_y)

    subplot(335)
    colorbar(min(sol_y),max(sol_y))
    a = get("current_axes");
    title = a.title;
    title.foreground = 9;
    title.font_size = 4;
    title.text = "$\textbf{u}_{y,solveur}(x,y)$"
    plot3d1(Y,X,sol_y)

    subplot(336)
    colorbar(min(abs(ref_y - sol_y)),max(abs(ref_y - sol_y)))
    a = get("current_axes");
    title = a.title;
    title.foreground = 9;
    title.font_size = 4;
    title.text = "$|\textbf{u}_{y,ref} - \textbf{u}_{y,solveur}|(x,y)$"
    plot3d1(Y,X,abs(ref_y - sol_y))

    subplot(338)
    colorbar(min(w),max(w))
    a = get("current_axes");
    title = a.title;
    title.foreground = 9;
    title.font_size = 4;
    title.text = "$\omega(x,y)$"
    plot3d1(Y,X,w)

    xs2png(fig, "poisson_curl_error.png")
endfunction


// Fonction de test pour le solveur curl de Poisson
function test_poisson_curl(Lx, Ly, Nx, Ny)
    printf("::Testing poisson curl operator::")
    printf("\n  Domain size:    [%0.2f, %0.2f]", Lx, Ly)
    printf("\n  Discretization: [%i, %i]", Nx, Ny)

    X = linspace(0.0, Lx*(Nx-1)/Nx, Nx)
    Y = linspace(0.0, Ly*(Ny-1)/Ny, Ny)

    printf("\n\n Initializing field w(x,y).")
    w = feval(Y, X, omega)

    printf("\n Initializing reference solution x-axis ref_x(x,y).")
    ref_x = feval(Y, X, solution_field_1)

    printf("\n Initializing reference solution y-axis ref_y(x,y).")
    ref_y = feval(Y, X, solution_field_2)

    dir = get_absolute_file_path("test_poisson_curl.sce")
    file = dir+"poisson.sce"
    printf("\n\n Loading poisson_curl_2d function from file %s%s%s.", char(39), file, char(39))
    exec(file, -1)

    printf("\n\n Computing Poisson solution on x-axis sol_x(x,y) \n and on y-axis sol_y(x,y).")
    [sol_x, sol_y] = poisson_curl_2d(w, Nx, Ny, Lx, Ly)

    printf("\n Computing error on x-axis |ref_x - sol_x|(x,y).")
    err_x = abs(ref_x - sol_x)

    printf("\n Computing error on y-axis |ref_y - sol_y|(x,y).")
    err_y = abs(ref_y - sol_y)

    file = pwd()+"/poisson_curl_error.png"
    printf("\n\n Plotting everything to %s%s%s.", char(39), file, char(39))
    plot_error(w, ref_x, ref_y, sol_x, sol_y)

    printf("\n\n")
    mErr_x = max(err_x)
    mErr_y = max(err_y)
    max_error = 1e-12

    if (mErr_x > max_error) then
        printf("   Maximal error on x-axis is %.10ef, TEST FAILURE (max_error=%.10ef).\n", mErr_x, max_error)
        exit(1)
    elseif (mErr_y > max_error) then
        printf("   Maximal error on y-axis is %.10ef, TEST FAILURE (max_error=%.10ef).\n", mErr_x, max_error)
        exit(1)
    else
        printf("   Maximal error is only %.10ef on x-axis and %.10ef on y-axis, TEST SUCCESS.\n", mErr_x, mErr_y)
        exit(0)
    end
endfunction


// Taille du domaine
Lx = 1.0
Ly = 1.0

// Discretisation du domaine
Nx = 64
Ny = 32

test_poisson_curl(Lx, Ly, Nx, Ny)
