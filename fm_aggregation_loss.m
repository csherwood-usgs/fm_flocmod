% fm_aggregation_loss
%**************************************************************************
%   Shear agregation : LOSS : f_l1
%**************************************************************************
for iv1=1:nv_mud
   for iv2=1:nv_mud     
      if(iv2 == iv1)
         mult=2.0;
      else
         mult=1.0;
      end
      f_l1_sh(iv2,iv1)=mult*f_alpha*f_coll_prob_sh(iv2,iv1);
      f_l1_ds(iv2,iv1)=mult*f_alpha*f_coll_prob_ds(iv2,iv1);      
   end
end
clear mult