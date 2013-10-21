%**************************************************************************
%  Shear fragmentation : LOSS : f_l2 (stet...f13)
%**************************************************************************
% TODO - this is easy to vectorize
for iv1=1:nv_mud
   if (f_diam(iv1)>dfragmax)
      % shear fragmentation
      f_l3(iv1)=f_l3(iv1)+(1.0-f_ero_frac)*f_beta*f_diam(iv1)*((f_diam(iv1)-f_dp0)/f_dp0)^(3.0-f_nf);
      % shear erosion
      if ((f_mass(iv1)-f_mass(f_ero_iv)*f_ero_nbfrag) > f_mass(f_ero_iv))
         f_l3(iv1)=f_l3(iv1)+f_ero_frac*f_beta*f_diam(iv1)*((f_diam(iv1)-f_dp0)/f_dp0)^(3.0-f_nf);
      end
   end
end
