% fm_aggregation_gain
%********************************************************************************
% agregation : GAIN : f_g1_sh and f_g1_ds
%********************************************************************************

% the Fortran indexing starting at 0 is hard to replicate in Matlab,
% so I have used diffmass and f_masslo to work around.
diffmass = diff([f_mass; 0]);
for iv1=1:nv_mud
   for iv2=1:nv_mud
      for iv3=iv2:nv_mud
         if(iv1==1)
            f_masslo = 0.0;
         else
            f_masslo = f_mass(iv1-1);
         end
         if((f_mass(iv2)+f_mass(iv3)) > f_masslo ...
               && ((f_mass(iv2)+f_mass(iv3)) <= f_mass(iv1)))
            
            %f_weight=(f_mass(iv2)+f_mass(iv3)-f_mass(iv1-1))/(f_mass(iv1)-f_mass(iv1-1));
            f_weight=(f_mass(iv2)+f_mass(iv3)-f_masslo)/(diffmass(iv1));
            
         elseif ((f_mass(iv2)+f_mass(iv3)) > f_mass(iv1) ...
               && ((f_mass(iv2)+f_mass(iv3)) < f_mass(iv1+1)))
            
            if (iv1 == nv_mud)
               f_weight=1.0;
            else
               %f_weight=1.0-(f_mass(iv2)+f_mass(iv3)-f_mass(iv1))/(f_mass(iv1+1)-f_mass(iv1));
               f_weight=1.0-(f_mass(iv2)+f_mass(iv3)-f_mass(iv1))/(diffmass(iv1+1));
            end
            
         else
            f_weight=0.0;
         end
         
         f_g1_sh(iv2,iv3,iv1)=f_weight*f_alpha*f_coll_prob_sh(iv2,iv3)*(f_mass(iv2)+f_mass(iv3))/f_mass(iv1);
         f_g1_ds(iv2,iv3,iv1)=f_weight*f_alpha*f_coll_prob_ds(iv2,iv3)*(f_mass(iv2)+f_mass(iv3))/f_mass(iv1);
         
      end
   end
end
clear f_weight