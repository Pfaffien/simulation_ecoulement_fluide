
// Domain size and discretization
Lx = 1.0
Ly = 1.0
Nx = 128
Ny = 128
dx = Lx/Nx
dy = Ly/Ny
X = linspace(0.0, Lx*(Nx-1)/Nx, Nx)
Y = linspace(0.0, Ly*(Ny-1)/Ny, Ny)
dt = 0.09

// Simulation parameters
T     = 1.50
nu    = 0.5e-4
rho   = 100.0
delta = 0.05

function u = f(t,xy)
  ux = Ux(xy(1),xy(2))
  uy = Uy(xy(1),xy(2))
  u = [ux,uy]
endfunction


// Initialize vorticity
function [W] = init_vorticity(y,x)
    // initialize vorticity W(x,y)
    d_y = [];
    if y <= 0.5 then
        d_y = rho / cosh(rho * (y - 0.25))^2;
    else
        d_y = -rho / cosh(rho * (0.75 - y))^2;
    end
    W = 2 * %pi * delta * cos(2 * %pi * x) - d_y
endfunction


// Plot and dump fields (at every steps of the simulation so that we can generate a video)
function plot_fields(W, Ux, Uy, iteration)
    fig = scf(0)
    clf()
    fig.color_map = [jetcolormap(64); whitecolormap(1)]
    fig.background = 65;
    cmm = [0,64]

    subplot(131)
    a = gca()
    a.isoview = "on"
    title = a.title
    title.foreground = 9
    title.font_size = 4
    title.font_size = 5
    title.text = "$\omega(x,y)$"
    colorbar(min(W), max(W), cmm)
    Sgrayplot(X, Y, W', colminmax=cmm)
    xlabel("x")
    ylabel("y")

    subplot(132)
    a = gca()
    a.isoview = "on"
    title = a.title
    title.foreground = 9
    title.font_size = 4
    title.font_size = 5
    title.text = "$U_x(x,y)$"
    colorbar(min(Ux), max(Ux), cmm)
    Sgrayplot(X, Y, Ux', colminmax=cmm)
    xlabel("x")
    ylabel("y")

    subplot(133)
    a = gca()
    a.isoview = "on"
    title = a.title
    title.foreground = 9
    title.font_size = 4
    title.font_size = 5
    title.text = "$U_y(x,y)$"
    colorbar(min(Uy), max(Uy), cmm)
    Sgrayplot(X, Y, Uy', colminmax=cmm)
    xlabel("x")
    ylabel("y")

    figname = sprintf("ite_%04d.png", iteration)
    xs2png(fig, figname)
endfunction


// Plot isocontours (only at t=0.80 and t=1.20)
function plot_isocontours(W, figname)
    if (t~=0.80) & (t~=1.20) then
        return
    end

    fig = scf(1)
    clf()
    xset("fpf",' ')

    // display the isocontours
    index = -70:10:70
    //index = -36:6:36
    contourf(X*Nx/Lx,Y*Ny/Ly,W,index)
    fchamp(f,0,1:8:Nx-1,1:8:Ny-1)

    figname = sprintf("isocontours_%f.png", t)
    xs2png(fig, figname)
endfunction


// Load the Poisson solver and the advection-diffusion solver
dir  = get_absolute_file_path("simu.sce")
file = dir+"poisson/poisson.sce"
exec(file, -1)
file = dir+"diff/dif-conv-f.sce"
exec(file, -1)

// Figure setup (fig0 = fields, fig1 = isocontours)
figure(0, "position", [0,0,1400,400])
figure(1, "position", [0,0,800,800])

// Initialize vorticity and loop untill final time using adaptive timestep
t = 0.0
ite = 0
W = feval(Y, X, init_vorticity)
Ux = ones(Ny,Nx)
Uy = zeros(Ny,Nx)
while t<T
    // compute velocity from vorticity
    [Ux,Uy] = poisson_curl_2d(W,Nx,Ny,Lx,Ly)

    dt = min(calcul_dt(Ux,dx),calcul_dt(Uy,dy))

    // compute new timestep from stability criteria
    if (t<0.80) & (t+dt>0.80) then
        dt = 0.80-t
    elseif (t<1.20) & (t+dt>1.20) then
        dt = 1.20-t
    elseif (t+dt>T) then
        dt = T-t
    end

    printf("\niteration %i, from t=%f to t=%f", ite, t, t+dt)
    plot_fields(W,Ux,Uy,ite)
    plot_isocontours(W,t)

    // advection-diffusion on vorticity
    W = solveur_2D(W,Ux,Uy,Nx,Ny,nu,dt,dx,dy)
    // update t and ite
    ite = ite+1
    t=t+dt

end
plot_fields(W,Ux,Uy,ite)
plot_isocontours(W,t)

printf("\nDone in %i iterations!\n", ite)
exit(0)
