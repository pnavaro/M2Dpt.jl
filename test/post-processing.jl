#    % Post-Processing
#    cpu = toc; Gbs = 8*(3*2)*nx*ny*iter/cpu/1e9;
#    fprintf('> it=%2.d, iter=%1.2e, err=%1.2e, max(T)=%1.2e, GB/s=%1.3f \n',[it,iter/1000,max(errs(end,1)),(max(T(:))),Gbs])
#    % Fields away from box corners
#    Tvis   = T;    Tvis  ((xc2.^2+yc2.^2)<(r+0.5*r)^2) = 0; Tvis  (((xc2-Lx).^2+(yc2-Ly).^2)<(r+0.5*r)^2) = 0;
#    Eiivis = Eii2; Eiivis((xc2.^2+yc2.^2)<(r+0.5*r)^2) = 0; Eiivis(((xc2-Lx).^2+(yc2-Ly).^2)<(r+0.5*r)^2) = 0;
#    % Conservation test
#    rhoc   = 1;
#    Edot   = rhoc*(T(:)-To(:))/dtT;
#    E      = E + dx*dy*dtT*sum(Edot);
#    % Volume integral
#    Wdot   = 4.*Eii2.*etac;
#    W      = W + dx*dy*dtT*sum(Wdot(:));
#    evol(:,it) = [min(Tvis(:)) max(Tvis(:)) min(Eiivis(:)) max(Eiivis(:)) min(Hs(:)) max(Hs(:)) min(etac(:)) max(etac(:)) iter time dtT cpu];
#    if mod(it,noutp)==0 || it==1
#        FS = 20;
#        figure(1),clf,colormap('jet'),set(gcf,'Color','white')
#        plot(1:nout:iter, log10(errs(:,1)'/errs(1,1)),'b.-',1:nout:iter, log10(errs(:,3)'/errs(1,3)),'r.-')
#        xlabel('Iterations', 'interpreter', 'latex', 'FontSize', FS), ylabel('$\| f_{u} \|_{\mathrm{rel}}$, $\| f_{T} \|_{\mathrm{rel}}$', 'interpreter', 'latex', 'FontSize', FS)
#        leg=legend('$\| f_{u} \|_{\mathrm{rel}}$', '$\| f_{T} \|_{\mathrm{rel}}$','Location', 'SouthWest'); legend boxoff, set(leg,'interpreter', 'latex'), set(gca, 'FontSize', FS,'Linewidth',1.6 )
#        drawnow
#        
#        figure(2),clf,colormap('jet'),set(gcf,'Color','white')
#        subplot(211), plot((1:it)*dtT*1e3, sqrt(evol(4,:))/(Vbc/Lx)', 'b.'),axis([0 3 0 16])
#        ylabel('$\max$ $\dot{\epsilon}_{II}$ / $\dot{\epsilon}_{BG}$', 'interpreter', 'latex', 'FontSize', FS), set(gca, 'FontSize', FS,'Linewidth',1.6 )
#        subplot(212), plot((1:it)*dtT*1e3, evol(2,:), 'b.'),axis([0 3 0.02 4])
#        ylabel('$\max$ $T$', 'interpreter', 'latex', 'FontSize', 20)
#        xlabel('$t \times 10^{-3}$', 'interpreter', 'latex', 'FontSize', 20), set(gca, 'FontSize', FS,'Linewidth',1.6 )
#        drawnow
#        
#        figure(3),hold on,colormap('jet'),set(gcf,'Color','white')
#        hold on, plot((it)*dtT*1e3, E, 'b.', (it)*dtT*1e3, W, 'dr'),box on
#        ylabel('$E, W$', 'interpreter', 'latex', 'FontSize', FS),xlabel('$t \times 10^{-3}$', 'interpreter', 'latex', 'FontSize', FS), axis([0 3 0 1])
#        leg=legend('E, DI.', 'W, DI.','Location', 'SouthEast'); legend boxoff, set(leg,'interpreter', 'latex'), set(gca, 'FontSize', FS,'Linewidth',1.6 )
#        drawnow
#        
#        figure(4),clf,colormap('jet'),set(gcf,'Color','white')
#        subplot(121), imagesc(xc,yc,flipud(   T')),colorbar,axis image,
#        title ('$T$', 'interpreter', 'latex', 'FontSize', FS),xlabel('$x$', 'interpreter', 'latex', 'FontSize', FS),ylabel('$y$', 'interpreter', 'latex', 'FontSize', FS), set(gca, 'FontSize', FS,'Linewidth',1.6 )
#        subplot(122), imagesc(xc,yc,flipud(log10(sqrt(Eii2))')),colorbar,axis image
#        title ('$\dot{\epsilon}_{II}$ / $\dot{\epsilon}_{BG}$', 'interpreter', 'latex', 'FontSize', FS),xlabel('$x$', 'interpreter', 'latex', 'FontSize', FS), set(gca, 'FontSize', FS,'Linewidth',1.6 )
#        drawnow



