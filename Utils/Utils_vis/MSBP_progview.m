function MSBP_progview(obj,fignum,plot_range,pltCurve,upto)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Displays a triframe cross-sectional view of a 3D volume, and may also
% plot a figure associated with a cost curve (optional)
%
% July 06, 2020 by Shwetadwip Chowdhury
% 
% inputs:   
%           obj:            3D matrix 
%           fignum:         number assigned to this figure
%           plot_range:     color axis with which to display values
%           pltCurve:       curve that may get plot (optional)
%           upto:           boundary up to which curve will be plot
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(fignum);
    subplot(2,2,1);
    imagesc(real(squeeze(obj(:,:,end/2)))); axis equal; axis tight;
    clim(plot_range); colormap gray; colorbar; title('x,y');
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    
    subplot(2,2,2);
    imagesc(real(squeeze(obj(:,end/2,:)))); axis equal; axis tight;
    clim(plot_range); colormap gray; colorbar; title('x,z');
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    
    subplot(2,2,3);
    imagesc(real(squeeze(obj(end/2,:,:)))); axis equal; axis tight;
    clim(plot_range); colormap gray; colorbar; title('y,z');
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    
    subplot(2,2,4);
    plot(log10(1+pltCurve),':o'); title('cost function'); axis tight;
    
    if nargin == 5
        subplot(2,2,4);
        start = 1;
        % if upto>50
        %     start = upto-50;
        % end
        plot(log10(1+pltCurve(start:upto)),':o'); title('cost function'); axis tight;
    end
end