% fm_shear_frag_gain
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
% Shear fragmentation : GAIN : f_g3
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

for iv1=1:nv_mud
   for iv2=iv1:nv_mud
      
      if(iv1==1)
         f_masslo = 0.0;
      else
         f_masslo = f_mass(iv1-1);
      end
      
      if (f_diam(iv2)>dfragmax)
         % binary fragmentation   
         if (f_mass(iv2)/f_nb_frag > f_masslo ...
               && f_mass(iv2)/f_nb_frag <= f_mass(iv1))
            
            if (iv1 == 1)
               f_weight=1.0;
            else
               f_weight=(f_mass(iv2)/f_nb_frag-f_masslo)/(f_mass(iv1)-f_masslo);
            end
         elseif (f_mass(iv2)/f_nb_frag > f_mass(iv1) ...
               && f_mass(iv2)/f_nb_frag < f_mass(iv1+1))
            f_weight=1.0-(f_mass(iv2)/f_nb_frag-f_mass(iv1))/(f_mass(iv1+1)-f_mass(iv1));
         else            
            f_weight=0.0;
         end
      else
         f_weight=0.0;
      end
      
      f_g3(iv2,iv1)=f_g3(iv2,iv1)+(1.0-f_ero_frac)*(1.0-f_ater)*f_weight*f_beta ...
         *f_diam(iv2)*((f_diam(iv2)-f_dp0)/f_dp0)^(3.0-f_nf)           ...
         *f_mass(iv2)/f_mass(iv1);
      
      % ternary fragmentation
      if (f_diam(iv2)>dfragmax)
         if (f_mass(iv2)/(2.0*f_nb_frag) > f_masslo ...
               && f_mass(iv2)/(2.0*f_nb_frag) <= f_mass(iv1))
            if (iv1 == 1)
               f_weight=1.0;
            else
               f_weight=(f_mass(iv2)/(2.0*f_nb_frag)-f_mass(iv1-1))/(f_mass(iv1)-f_masslo);
            end
         elseif (f_mass(iv2)/(2.0*f_nb_frag) > f_mass(iv1) ...
               && f_mass(iv2)/(2.0*f_nb_frag) < f_mass(iv1+1))
            f_weight=1.0-(f_mass(iv2)/(2.0*f_nb_frag)-f_mass(iv1))/(f_mass(iv1+1)-f_mass(iv1));            
         else
            f_weight=0.0;
            
         end
         % update for ternary fragments
         f_g3(iv2,iv1)=f_g3(iv2,iv1)+(1.0-f_ero_frac)*(f_ater)*f_weight*f_beta ...
            *f_diam(iv2)*((f_diam(iv2)-f_dp0)/f_dp0)^(3.0-f_nf)           ...
            *f_mass(iv2)/f_mass(iv1);
         
         % Floc erosion
         if ((f_mass(iv2)-f_mass(f_ero_iv)*f_ero_nbfrag) > f_mass(f_ero_iv))
            if (((f_mass(iv2)-f_mass(f_ero_iv)*f_ero_nbfrag) >f_masslo) ...
                  && (f_mass(iv2)-f_mass(f_ero_iv)*f_ero_nbfrag) <= f_mass(iv1))              
               if (iv1 == 1)
                  f_weight=1.0;
               else
                  f_weight=(f_mass(iv2)-f_mass(f_ero_iv)*f_ero_nbfrag-f_masslo)/(f_mass(iv1)-f_masslo);
               end
            elseif ((f_mass(iv2)-f_mass(f_ero_iv)*f_ero_nbfrag) > f_mass(iv1) ...
                  && (f_mass(iv2)-f_mass(f_ero_iv)*f_ero_nbfrag) < f_mass(iv1+1))
               f_weight=1.0-(f_mass(iv2)-f_mass(f_ero_iv)*f_ero_nbfrag-f_mass(iv1))/(f_mass(iv1+1)-f_mass(iv1));       
            else
               f_weight=0.0;
            end
            
            % update for eroded floc masses            
            f_g3(iv2,iv1)=f_g3(iv2,iv1)+f_ero_frac*f_weight*f_beta                    ...
               *f_diam(iv2)*((f_diam(iv2)-f_dp0)/f_dp0)^(3.0-f_nf)           ...
               *(f_mass(iv2)-f_mass(f_ero_iv)*f_ero_nbfrag)/f_mass(iv1);
            
            if (iv1 == f_ero_iv)               
               f_g3(iv2,iv1)=f_g3(iv2,iv1)+f_ero_frac*f_beta                           ...
                  *f_diam(iv2)*((f_diam(iv2)-f_dp0)/f_dp0)^(3.0-f_nf)           ...                                      ...
                  *f_ero_nbfrag*f_mass(f_ero_iv)/f_mass(iv1);
            end
         end
      end % condition on dfragmax
   end
end
clear f_weight
