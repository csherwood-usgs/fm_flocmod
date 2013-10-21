% fm_comp_fsd
% This processes NNin and returns NNout
NNout = zeros(size(NNin));
tmp_g1=0.0;
tmp_g3=0.0;
tmp_g4=0.0;
tmp_l1=0.0;
tmp_l3=0.0;
tmp_l4=0.0;
f_g1_tmp=zeros(nv_mud,nv_mud,nv_mud);
f_l1_tmp=zeros(nv_mud,nv_mud);

if (l_COLLFRAG)
   fm_collfrag
end

for iv1=1:nv_mud
   for iv2=1:nv_mud
      for iv3=1:nv_mud
         if (l_ASH)
            f_g1_tmp(iv2,iv3,iv1)=f_g1_tmp(iv2,iv3,iv1)+f_g1_sh(iv2,iv3,iv1)*Gval;
         end
         if (l_ADS)
            f_g1_tmp(iv2,iv3,iv1)=f_g1_tmp(iv2,iv3,iv1)+f_g1_ds(iv2,iv3,iv1);
         end
         
         tmp_g1=tmp_g1+(NNin(iv3)*(f_g1_tmp(iv2,iv3,iv1))*NNin(iv2));
         
         if (l_COLLFRAG)
            tmp_g4=tmp_g4+(NNin(iv3)*(f_g4(iv2,iv3,iv1)*Gval)*NNin(iv2));
         end
      end
      
      tmp_g3=tmp_g3+f_g3(iv2,iv1)*NNin(iv2)*Gval^1.5;
      
      if (l_ASH)
         f_l1_tmp(iv2,iv1)=f_l1_tmp(iv2,iv1)+f_l1_sh(iv2,iv1)*Gval;
      end
      if (l_ADS)
         f_l1_tmp(iv2,iv1)=f_l1_tmp(iv2,iv1)+f_l1_ds(iv2,iv1)*Gval;
      end
      
      tmp_l1=tmp_l1+(f_l1_tmp(iv2,iv1))*NNin(iv2);
      
      if (l_COLLFRAG)
         tmp_l4=tmp_l4+(f_l4(iv2,iv1)*Gval)*NNin(iv2);
      end
   end
   
   tmp_l1=tmp_l1*NNin(iv1);
   tmp_l4=tmp_l4*NNin(iv1);
   
   tmp_l3=f_l3(iv1)*Gval^1.50*NNin(iv1);
   
   NNout(iv1)=NNin(iv1)+f_dt*(tmp_g1+tmp_g3+tmp_g4-(tmp_l1+tmp_l3+tmp_l4));
   
   tmp_g1=0.0;
   tmp_g3=0.0;
   tmp_g4=0.0;
   tmp_l1=0.0;
   tmp_l3=0.0;
   tmp_l4=0.0;
end
