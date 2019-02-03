Physics = Dict(
n      => 3,                  % stress exponent for power law rheology
Vbc    => 66.4437;            % boundary velocity difference
Lx     => 0.86038;            % domain size x
Ly     => Lx;                 % domain size y
T0     => 49.3269/n;          % initial temperature
r      => 0.0737;             % radius of initial perturbation
Tamp   => 0.1*T0;             % amplitude of initial perturbation
% Numerics
nx     = 94;                 % number of cells x
ny     = nx;                 % number of cells y
nt     = 54;                 % number time steps
nout   = 100;                % check residual each nout iteration
noutp  = 2;                  % display graphics each nout time step
niter  = 1e5;                % max nonlinear iterations
epsi   = 1e-5;               % non-linear tolerance
tetp   = 0.5;                % reduction of PT steps for pressure
tetv   = 0.5;                % reduction of PT steps for velocity
tetT   = 0.5;                % reduction of physical time step for temperature
rel    = 1e-1;               % relaxation factor for non-linear viscosity
Vdamp  = 4.0;                % velocity damping for momentum equations
eta_b  = 1.0;                % numerical compressibility
% Pre-processing
dampx  = 1*(1-Vdamp/nx);     % velocity damping for x-momentum equation
dampy  = 1*(1-Vdamp/ny);     % velocity damping for y-momentum equation
mpow   = -(1-1/n)/2;         % exponent for strain rate dependent viscosity
dx     = Lx/nx;              % grid step in x
dy     = Ly/ny;              % grid step in y
time   = 0;
% Mesh
xn = 0:dx:Lx; xc = dx/2:dx:Lx-dx/2;
yn = 0:dy:Ly; yc = dy/2:dy:Ly-dy/2;
[xc2,  yc2] = ndgrid(xc,yc);
[xvx2,yvx2] = ndgrid(xn,yc);
[xvy2,yvy2] = ndgrid(xc,yn);
% Intial fields
Vx         =  Vbc*xvx2/Lx;
Vy         = -Vbc*yvy2/Ly;
T          =  zeros(nx  ,ny  );
P          =  zeros(nx  ,ny  );
etac       =   ones(nx  ,ny  );
qx         =  zeros(nx+1,ny  );
qy         =  zeros(nx  ,ny+1);
dVxdtauVx  =  zeros(nx-1,ny  );
dVydtauVy  =  zeros(nx  ,ny-1);
dVxdtauVx0 =  zeros(nx-1,ny  );
dVydtauVy0 =  zeros(nx  ,ny-1);
T((xc2.^2+yc2.^2)<r^2) = Tamp;                                             % initial temperature pertubation
dtT        = tetT*1/4.1*min(dx,dy)^2;                                      % explicit timestep for 2D diffusion
E = 0; W = 0;
% Action
for it = 1:nt % ------ Physical timesteps
    To    = T;                                                             % temperature from previous step (for backward-Euler integration)
    time  = time + dtT;                                                    % update physical time
    errs  = []; tic
    for iter = 1:niter % ------ Pseudo-Transient cycles
        err   = [Vx(:); Vy(:); P(:); T(:); etac(:)];
        dVxdtauVx0      = dVxdtauVx + dampx.*dVxdtauVx0;                   % used for damping x momentum residuals
        dVydtauVy0      = dVydtauVy + dampy.*dVydtauVy0;                   % used for damping y momentum residuals
        % ------ Kinematics
        Vx_exp          = [ Vx(:,1), Vx, Vx(:,end) ];                      % expand array using BC's - Free slip
        Vy_exp          = [ Vy(1,:); Vy; Vy(end,:) ];
        divV            = diff(Vx,1,1)/dx + diff(Vy,1,2)/dy;
        Exxc            = diff(Vx,1,1)/dx - 1/2*divV;
        Eyyc            = diff(Vy,1,2)/dy - 1/2*divV;
        Exyv            = 0.5*(diff(Vx_exp,1,2)/dy + diff(Vy_exp,1,1)/dx);
        Exyc            = 0.25*(Exyv(1:end-1,1:end-1) + Exyv(2:end,1:end-1) + Exyv(1:end-1,2:end) + Exyv(2:end,2:end));
        Eii2            = 0.5*(Exxc.^2 + Eyyc.^2) + Exyc.^2;               % strain rate invariant
        % ------ Rheology
        etac_phys       = Eii2.^mpow.*exp( -T.*(1./(1+T./T0)) );           % physical viscosity
        etac            = exp(rel*log(etac_phys) + (1-rel)*log(etac));     % numerical shear viscosity
        etav            = zeros(size(etac,1) + [1 1]);                     % expand viscosity fom cell centroids to vertices
        etav(2:end-1,2:end-1) = 0.25*(etac(1:end-1,1:end-1) + etac(2:end,2:end) + etac(1:end-1,2:end) + etac(2:end,1:end-1));
        etav(:,[1 end]) = etav(:,[2 end-1]);
        etav([1 end],:) = etav([2 end-1],:);
        % ------ Pseudo-Time steps
        dtauP           = tetp*  4.1/min(nx,ny)*etac*(1.0+eta_b);
        dtauVx          = tetv*1/4.1*(min(dx,dy)^2./( 0.5*(etac(2:end,:) + etac(1:end-1,:)) ))/(1+eta_b);
        dtauVy          = tetv*1/4.1*(min(dx,dy)^2./( 0.5*(etac(:,2:end) + etac(:,1:end-1)) ))/(1+eta_b);
        dtauT           = tetT*1/4.1*min(dx,dy)^2;
        % ------ Fluxes
        qx(2:end-1,:)   = -diff(T,1,1)/dx;
        qy(:,2:end-1)   = -diff(T,1,2)/dy;
        Sxx             = -P + 2*etac.*(Exxc + eta_b*divV);
        Syy             = -P + 2*etac.*(Eyyc + eta_b*divV);
        Txy             =      2*etav.*Exyv;
        Hs              = 4*etac.*Eii2;
        % ------ Residuals
        dVxdtauVx       =               diff(Txy(2:end-1,:),1,2)/dy + diff(Sxx,1,1)/dx;
        dVydtauVy       =               diff(Txy(:,2:end-1),1,1)/dx + diff(Syy,1,2)/dy;
        dPdtauP         =            -  divV;
        dTdtauT         = (To-T)/dtT - (diff(qx,1,1)/dx + diff(qy,1,2)/dy) + Hs;
        % ------ Updates
        Vx(2:end-1,:)   = Vx(2:end-1,:) + dtauVx.*(dVxdtauVx + dampx.*dVxdtauVx0); % update with damping
        Vy(:,2:end-1)   = Vy(:,2:end-1) + dtauVy.*(dVydtauVy + dampy.*dVydtauVy0); % update with damping
        P               = P             + dtauP .*dPdtauP;
        T               = T             + dtauT .*dTdtauT;
        if mod(iter,nout)==0 % -------------------------------------- Check
            fu   = [dVxdtauVx(:); dVydtauVy(:)]; fp = dPdtauP(:); fT = dTdtauT(:);
            errs = [errs; norm(fu)/length(fu), norm(fp(:))/length(fp(:)), norm(fT(:))/length(fT(:))];
            if max(errs(end,:))<epsi, break; end
            fprintf('iter = %d\nf_{u} = %1.3e\nf_{p} = %1.3e\nf_{T} = %1.3e\n',iter, errs(end,1), errs(end,2), errs(end,3) )
        end%mod
    end%niter
    % Post-Processing
    cpu = toc; Gbs = 8*(3*2)*nx*ny*iter/cpu/1e9;
    fprintf('> it=%2.d, iter=%1.2e, err=%1.2e, max(T)=%1.2e, GB/s=%1.3f \n',[it,iter/1000,max(errs(end,1)),(max(T(:))),Gbs])
    % Fields away from box corners
    Tvis   = T;    Tvis  ((xc2.^2+yc2.^2)<(r+0.5*r)^2) = 0; Tvis  (((xc2-Lx).^2+(yc2-Ly).^2)<(r+0.5*r)^2) = 0;
    Eiivis = Eii2; Eiivis((xc2.^2+yc2.^2)<(r+0.5*r)^2) = 0; Eiivis(((xc2-Lx).^2+(yc2-Ly).^2)<(r+0.5*r)^2) = 0;
    % Conservation test
    rhoc   = 1;
    Edot   = rhoc*(T(:)-To(:))/dtT;
    E      = E + dx*dy*dtT*sum(Edot);
    % Volume integral
    Wdot   = 4.*Eii2.*etac;
    W      = W + dx*dy*dtT*sum(Wdot(:));
    evol(:,it) = [min(Tvis(:)) max(Tvis(:)) min(Eiivis(:)) max(Eiivis(:)) min(Hs(:)) max(Hs(:)) min(etac(:)) max(etac(:)) iter time dtT cpu];
    if mod(it,noutp)==0 || it==1
        FS = 20;
        figure(1),clf,colormap('jet'),set(gcf,'Color','white')
        plot(1:nout:iter, log10(errs(:,1)'/errs(1,1)),'b.-',1:nout:iter, log10(errs(:,3)'/errs(1,3)),'r.-')
        xlabel('Iterations', 'interpreter', 'latex', 'FontSize', FS), ylabel('$\| f_{u} \|_{\mathrm{rel}}$, $\| f_{T} \|_{\mathrm{rel}}$', 'interpreter', 'latex', 'FontSize', FS)
        leg=legend('$\| f_{u} \|_{\mathrm{rel}}$', '$\| f_{T} \|_{\mathrm{rel}}$','Location', 'SouthWest'); legend boxoff, set(leg,'interpreter', 'latex'), set(gca, 'FontSize', FS,'Linewidth',1.6 )
        drawnow
        
        figure(2),clf,colormap('jet'),set(gcf,'Color','white')
        subplot(211), plot((1:it)*dtT*1e3, sqrt(evol(4,:))/(Vbc/Lx)', 'b.'),axis([0 3 0 16])
        ylabel('$\max$ $\dot{\epsilon}_{II}$ / $\dot{\epsilon}_{BG}$', 'interpreter', 'latex', 'FontSize', FS), set(gca, 'FontSize', FS,'Linewidth',1.6 )
        subplot(212), plot((1:it)*dtT*1e3, evol(2,:), 'b.'),axis([0 3 0.02 4])
        ylabel('$\max$ $T$', 'interpreter', 'latex', 'FontSize', 20)
        xlabel('$t \times 10^{-3}$', 'interpreter', 'latex', 'FontSize', 20), set(gca, 'FontSize', FS,'Linewidth',1.6 )
        drawnow
        
        figure(3),hold on,colormap('jet'),set(gcf,'Color','white')
        hold on, plot((it)*dtT*1e3, E, 'b.', (it)*dtT*1e3, W, 'dr'),box on
        ylabel('$E, W$', 'interpreter', 'latex', 'FontSize', FS),xlabel('$t \times 10^{-3}$', 'interpreter', 'latex', 'FontSize', FS), axis([0 3 0 1])
        leg=legend('E, DI.', 'W, DI.','Location', 'SouthEast'); legend boxoff, set(leg,'interpreter', 'latex'), set(gca, 'FontSize', FS,'Linewidth',1.6 )
        drawnow
        
        figure(4),clf,colormap('jet'),set(gcf,'Color','white')
        subplot(121), imagesc(xc,yc,flipud(   T')),colorbar,axis image,
        title ('$T$', 'interpreter', 'latex', 'FontSize', FS),xlabel('$x$', 'interpreter', 'latex', 'FontSize', FS),ylabel('$y$', 'interpreter', 'latex', 'FontSize', FS), set(gca, 'FontSize', FS,'Linewidth',1.6 )
        subplot(122), imagesc(xc,yc,flipud(log10(sqrt(Eii2))')),colorbar,axis image
        title ('$\dot{\epsilon}_{II}$ / $\dot{\epsilon}_{BG}$', 'interpreter', 'latex', 'FontSize', FS),xlabel('$x$', 'interpreter', 'latex', 'FontSize', FS), set(gca, 'FontSize', FS,'Linewidth',1.6 )
        drawnow
    end
end



