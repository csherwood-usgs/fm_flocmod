% flocmod_agregation_statistics
for iv1=1:nv_mud
   for iv2=1:nv_mud
      % Shear (eqn 9)
      f_coll_prob_sh(iv1,iv2)=1.0/6.0*(f_diam(iv1)+f_diam(iv2))^3.0;
      % Differential settling
      f_coll_prob_ds(iv1,iv2)=0.250*pi*(f_diam(iv1)+f_diam(iv2))^2.0 ...
         *grav/mu*abs((f_rho(iv1)-rhoref)*f_diam(iv1)^2.0 ...
         -(f_rho(iv2)-rhoref)*f_diam(iv2)^2.0);      
   end
end
