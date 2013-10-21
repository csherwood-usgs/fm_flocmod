% fm_collfrag
f_fp=0.10
f_fy=1e-10
f_cfcst=3.0/16.0
f_g4=zeros(nv_mud,nv_mud,nv_mud);
f_l4=zeros(nv_mud,nv_mud);

for iv1=1:nv_mud
   for iv2=1:nv_mud
      for iv3=iv2:nv_mud
         % fragmentation after collision probability based on Gval for particles iv2 and iv3
         % gcolfrag=(collision induced shear) / (floc strength)
         gcolfragmin=2.0*(Gval*(f_diam(iv2)+f_diam(iv3))).^2.0*f_mass(iv2)*f_mass(iv3)  ...
            /(pi*f_fy*f_fp*f_diam(iv3).^2.0*(f_mass(iv2)+f_mass(iv3))         ...
            *((f_rho(iv3)-rhoref)/rhoref).^(2.0/(3.0-f_nf)));
         
         gcolfragmax=2.0*(Gval*(f_diam(iv2)+f_diam(iv3))).^2.0*f_mass(iv2)*f_mass(iv3)  ...
            /(pi*f_fy*f_fp*f_diam(iv2).^2.0*(f_mass(iv2)+f_mass(iv3))         ...
            *((f_rho(iv2)-rhoref)/rhoref).^(2.0/(3.0-f_nf)));
         
         % first case : iv3 not eroded, iv2 eroded forming 2 particles : iv3+f_cfcst*iv2 / iv2-f_cfcst*iv2
         if (gcolfragmin<1.0 && gcolfragmax>=10)
            if (((f_mass(iv3)+f_cfcst*f_mass(iv2))>f_mass(iv1-1)) &&  ...
                  ((f_mass(iv3)+f_cfcst*f_mass(iv2))<=f_mass(iv1)))             
               f_weight=((f_mass(iv3)+f_cfcst*f_mass(iv2)-f_mass(iv1-1))/(f_mass(iv1)-f_mass(iv1-1)));
            elseif (f_mass(iv3)+f_cfcst*f_mass(iv2)>f_mass(iv1)  && ...
                  f_mass(iv3)+f_cfcst*f_mass(iv2)<f_mass(iv1+1));
               if (iv1.eq.nv_mud)
                  f_weight=1.0;
               else
                  f_weight=1.0-((f_mass(iv3)+f_cfcst*f_mass(iv2)-f_mass(iv1))/(f_mass(iv1+1)-f_mass(iv1)));
               end   
            else
               f_weight=0.0;
            end
            f_g4(iv2,iv3,iv1)=f_g4(iv2,iv3,iv1)+f_weight*(f_coll_prob_sh(iv2,iv3))   ...
               *(f_mass(iv3)+f_cfcst*f_mass(iv2))/f_mass(iv1);
            
            if (f_mass(iv2)-f_cfcst*f_mass(iv2)>f_mass(iv1-1)   && ...
                  f_mass(iv2)-f_cfcst*f_mass(iv2)<=f_mass(iv1));
               f_weight=((f_mass(iv2)-f_cfcst*f_mass(iv2)-f_mass(iv1-1))/(f_mass(iv1)-f_mass(iv1-1)));
            elseif (f_mass(iv2)-f_cfcst*f_mass(iv2)>f_mass(iv1)  &&  ...
                  f_mass(iv2)-f_cfcst*f_mass(iv2)<f_mass(iv1+1));               
               if (iv1.eq.nv_mud)
                  f_weight=1.0;
               else
                  f_weight=1.0-((f_mass(iv2)-f_cfcst*f_mass(iv2)-f_mass(iv1))/(f_mass(iv1+1)-f_mass(iv1)));
               end
            else
               f_weight=0.0;
            end
            f_g4(iv2,iv3,iv1)=f_g4(iv2,iv3,iv1)+f_weight*(f_coll_prob_sh(iv2,iv3))   ...
               *(f_mass(iv2)-f_cfcst*f_mass(iv2))/f_mass(iv1);   
         % second case : iv3 eroded and iv2 eroded forming 3 particles : iv3-f_cfcst*iv3 / iv2-f_cfcst*iv2 / f_cfcst*iv3+f_cfcst*iv2
         elseif (gcolfragmin>=1.0 && gcolfragmax>=10)   % iv2 and iv3 eroded forming new (third) particle
            if (f_cfcst*f_mass(iv2)+f_cfcst*f_mass(iv3)>f_mass(iv1-1) &&  ...
                  f_cfcst*f_mass(iv2)+f_cfcst*f_mass(iv3)<=f_mass(iv1))
               f_weight=((f_cfcst*f_mass(iv2)+f_cfcst*f_mass(iv3)-f_mass(iv1-1))/(f_mass(iv1)-f_mass(iv1-1)));
            elseif (f_cfcst*f_mass(iv2)+f_cfcst*f_mass(iv3)>f_mass(iv1) &&  ...
                  f_cfcst*f_mass(iv2)+f_cfcst*f_mass(iv3)<f_mass(iv1+1));
               if (iv1==nv_mud)
                  f_weight=1.0;
               else
                  f_weight=1.0-((f_cfcst*f_mass(iv2)+f_cfcst*f_mass(iv3)-f_mass(iv1))/(f_mass(iv1+1)-f_mass(iv1)));
               end
            else
               f_weight=0.0;
            end
            f_g4(iv2,iv3,iv1)=f_g4(iv2,iv3,iv1)+f_weight*(f_coll_prob_sh(iv2,iv3))   ...
               *(f_cfcst*f_mass(iv2)+f_cfcst*f_mass(iv3))/f_mass(iv1);
            if ((1.0-f_cfcst)*f_mass(iv2)>f_mass(iv1-1) &&  ...
                  (1.0-f_cfcst)*f_mass(iv2)<=f_mass(iv1))
               f_weight=((1.0-f_cfcst)*f_mass(iv2)-f_mass(iv1-1))/(f_mass(iv1)-f_mass(iv1-1)); 
            elseif ((1.0-f_cfcst)*f_mass(iv2)>f_mass(iv1) &&  ...
                  (1.0-f_cfcst)*f_mass(iv2)<f_mass(iv1+1))   
               if (iv1==nv_mud)
                  f_weight=1.0;
               else
                  f_weight=1.0-(((1.0-f_cfcst)*f_mass(iv2)-f_mass(iv1))/(f_mass(iv1+1)-f_mass(iv1)));
               end            
            else
               f_weight=0.0;
            end           
            f_g4(iv2,iv3,iv1)=f_g4(iv2,iv3,iv1)+f_weight*(f_coll_prob_sh(iv2,iv3))   ...
               *((1.0-f_cfcst)*f_mass(iv2))/f_mass(iv1);
            if ((1.0-f_cfcst)*f_mass(iv3)>f_mass(iv1-1) &&  ...
                  (1.0-f_cfcst)*f_mass(iv3)<=f_mass(iv1))              
               f_weight=((1.0-f_cfcst)*f_mass(iv3)-f_mass(iv1-1))/(f_mass(iv1)-f_mass(iv1-1));
            elseif ((1.0-f_cfcst)*f_mass(iv3)>f_mass(iv1) &&  ...
                  (1.0-f_cfcst)*f_mass(iv3)<f_mass(iv1+1));
               if (iv1.eq.nv_mud)
                  f_weight=1.0;
               else
                  f_weight=1.0-(((1.0-f_cfcst)*f_mass(iv3)-f_mass(iv1))/(f_mass(iv1+1)-f_mass(iv1)));
               end;               
            else
               f_weight=0.0;
            end
            f_g4(iv2,iv3,iv1)=f_g4(iv2,iv3,iv1)+f_weight*(f_coll_prob_sh(iv2,iv3))   ...
               *((1.0-f_cfcst)*f_mass(iv3))/f_mass(iv1);
         end % end collision test case
      end
   end
end

for iv1=1:nv_mud
   for iv2=1:nv_mud
      
      gcolfragiv1=2.0*(Gval*(f_diam(iv1)+f_diam(iv2))).^2.0*f_mass(iv1)*f_mass(iv2)  ...
         /(pi*f_fy*f_fp*f_diam(iv1).^2.0*(f_mass(iv1)+f_mass(iv2))         ...
         *((f_rho(iv1)-rhoref)/rhoref).^(2.0/(3.0-f_nf)));
      
      gcolfragiv2=2.0*(Gval*(f_diam(iv1)+f_diam(iv2))).^2.0*f_mass(iv1)*f_mass(iv2)  ...
         /(pi*f_fy*f_fp*f_diam(iv2).^2.0*(f_mass(iv1)+f_mass(iv2))         ...
         *((f_rho(iv2)-rhoref)/rhoref).^(2.0/(3.0-f_nf)));
      
      mult=1.0;
      if (iv1.eq.iv2); mult=2.0; end
      if (iv1.eq.MAX(iv1,iv2) && gcolfragiv1>=1.0)
         f_l4(iv2,iv1)=f_l4(iv2,iv1)+mult*(f_coll_prob_sh(iv1,iv2));
      elseif (iv1.eq.MIN(iv1,iv2) && gcolfragiv2>=1.0)
         f_l4(iv2,iv1)=f_l4(iv2,iv1)+mult*(f_coll_prob_sh(iv1,iv2));
      end
   end
end

f_g4(1:nv_mud,1:nv_mud,1:nv_mud)=f_g4(1:nv_mud,1:nv_mud,1:nv_mud)*f_collfragparam;
f_l4(1:nv_mud,1:nv_mud)=f_l4(1:nv_mud,1:nv_mud)*f_collfragparam;
