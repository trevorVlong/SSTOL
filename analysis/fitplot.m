function [done] = fitplot(fig,plotnum,spm,spn,xl,yl,xlbl,ylbl,ttl,alfavec,cmat,DCJmod)
%visplot plots line data from the fitting of coefficient curves
%   Detailed explanation goes here

    figure(fig)
    cax = [0 25];
    subplot(spm,spn,plotnum);
    cmap = jet(length(DCJmod));
    alfamat = alfavec.*ones(size(cmat));
    
    for ii = 1:2:length(DCJmod)
        plot(alfamat(:,ii),cmat(:,ii),'Color',cmap(ii,:));
        hold on
    end
    colormap('jet')
    caxis(cax);    
    grid on
    colorbar
    title(ttl)
    xlabel(xlbl)
    ylabel(ylbl)
    xlim(xl)
    ylim(yl)
    axis fill
    %axis square
    done = 1;
end

