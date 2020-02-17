clear%,close all
folder="SIM02_spiral_outcoupler/far_field_data/";
%betterSource_charge
for top_charge=1
load(strcat(folder,"spiral_outcoupler","_charge",string(top_charge),"_539"));
% load(strcat(folder,"SPIRALE_PEC_nonaz"))CloserSource_smalleRadius
% useful_ux = abs(ux)<0.2;
% useful_uy = abs(uy)<0.2;
% ux = ux(useful_ux,useful_uy);
% uy = uy(useful_ux,useful_uy);
% E_phi = E_phi(useful_ux,useful_uy);
% E_theta = E_theta(useful_ux,useful_uy);

E_theta2 = zeros(size(E_phi));
E_phi2 = E_theta2;
 for lambda = [ 524,531, 539]
load(strcat(folder,"spiral_outcoupler","_charge",string(top_charge),"_",string(lambda)));
E_theta2 = E_theta2 + E_theta;
E_phi2 = E_phi2 + E_phi;
 end

 E_theta=E_theta2/3;
 E_phi=E_phi2/3;
 
 [ux,uy]=meshgrid(ux,uy);
 
% cos(thneta) = uz
theta = real( acos( sqrt(1 - ux.^2 - uy.^2)));
cos_phi = ux./sin(theta);
sin_phi = uy./sin(theta);

Ex = + E_phi.*sin_phi + E_theta.*cos_phi;%
Ey = - E_phi.*cos_phi + E_theta.*sin_phi;%

EL = sqrt(2)/2*( Ex + Ey*exp(+1i*pi/2));
ER = sqrt(2)/2*(-Ex + Ey*exp(+1i*pi/2));

S3 = 1i*(Ex.*conj(Ey)-Ey.*conj(Ex));
chi = 0.5*asin( real(S3)./(abs(Ex).^2+abs(Ey).^2));

% figure
% plot_surf(ux,uy,sin_phi,"something about e_phi");
% plot_surf(ux,uy,theta,"something about e_theta");

% figure
% subplot(1,2,2)
% plot_surf(ux,uy,abs(Ey).^2   +abs(Ex).^2,'Intensity from Ex and Ey');
% subplot(1,2,1)
% plot_surf(ux,uy,abs(E_phi).^2 + abs(E_theta).^2,'Intensity from E_\phi and E_\theta');
% %%
figure
subplot(2,3,1)
plot_surf(ux,uy,abs(ER).^2,"Right circular polarization intensity");
subplot(2,3,2)
plot_surf(ux,uy,abs(EL).^2,"Left circular polarization intensity");
subplot(2,3,4)
plot_surf(ux,uy,angle(ER),"Right circular polarization phase");
subplot(2,3,5)
plot_surf(ux,uy,angle(EL),"Left circular polarization phase");

subplot(2,3,3)
plot_surf(ux,uy,real(S3),"S3 Stokes parameter");
subplot(2,3,6)
plot_surf(ux,uy,abs(tan(chi)),"eccentricity |tan\chi|");
%  end
end



function plot_surf(ux,uy,Quantity,picture_title)
    s=surf(ux,uy,Quantity);
    s.EdgeColor='none';
    colormap('hot')
    colorbar
    view(2)
    xlabel("ux");
    ylabel('uy');
    axis('square')
    title(picture_title)
end