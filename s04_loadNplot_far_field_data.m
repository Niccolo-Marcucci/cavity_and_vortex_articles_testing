clear, close all
load far_field_data
% 
% useful_ux = abs(ux)<0.2;
% useful_uy = abs(uy)<0.2;

[ux,uy]=meshgrid(ux,uy);
% ux = ux(useful_ux,useful_uy);
% uy = uy(useful_ux,useful_uy);
% E_phi = E_phi(useful_ux,useful_uy);
% E_theta = E_theta(useful_ux,useful_uy);


% cos(theta) = uz
theta = real( acos( sqrt(1 - ux.^2 - uy.^2)));
cos_phi = ux./sin(theta);
sin_phi = uy./sin(theta);

Ex = E_theta.*cos_phi + E_phi.*sin_phi;
Ey = E_theta.*sin_phi + E_phi.*cos_phi;

EL = +sqrt(2)/2*Ex + sqrt(2)/2*Ey*exp(+1i*pi/2);
ER = -sqrt(2)/2*Ex + sqrt(2)/2*Ey*exp(+1i*pi/2);


% figure
% plot_surf(ux,uy,real(E_phi),"something about e_phi");
% plot_surf(ux,uy,real(E_theta),"something about e_theta");
% plot_surf(ux,uy,abs(Ey).^2   +abs(Ex).^2,'intensity from Ex and Ey');
% plot_surf(ux,uy,abs(E_phi).^2+abs(E_theta).^2,'intensity from E_phi and E_theta');

figure
subplot(2,2,1)
plot_surf(ux,uy,abs(ER).^2,"Right circular polarization intensity");
subplot(2,2,2)
plot_surf(ux,uy,abs(EL).^2,"Left circular polarization intensity");
subplot(2,2,3)
plot_surf(ux,uy,angle(ER),"Right circular polarization phase");
subplot(2,2,4)
plot_surf(ux,uy,angle(EL),"Left circular polarization phase");

%%
S3 = 1i*(Ex.*conj(Ey)-Ey.*conj(Ex));
chi = 0.5*asin( real(S3)./(abs(Ex).^2+abs(Ey).^2));

figure
subplot(1,2,1)
plot_surf(ux,uy,real(S3),"S3 Stokes parameter");
subplot(1,2,2)
plot_surf(ux,uy,abs(tan(chi)),"eccentricity |tan\chi|");
    

function plot_surf(ux,uy,Quantity,picture_title)
    s=surf(ux,uy,Quantity);
    s.EdgeColor='none';
    colormap('jet')
    colorbar
    view(2)
    xlabel("ux");
    ylabel('uy');
    axis('square')
    title(picture_title)
end