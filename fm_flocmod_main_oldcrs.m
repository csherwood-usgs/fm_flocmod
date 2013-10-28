% fm_flocmod_main
clear
mu = 0.0010
grav = 9.81
rhoref = 1030;
rhosp = 2650;
f_dp0 = 4e-6;
f_nf = 1.9
l_ADS=0
l_ASH=1
l_COLLFRAG=0
f_dmax=0.001500
f_nb_frag=2.
f_alpha=0.35
f_beta=0.15
f_ater=0.
f_ero_frac=0.0
f_ero_nbfrag=2.
f_ero_iv=1
f_mneg_param=0.000
f_collfragparam=0.01
f_test=1
dfragmax=0.00003
epsilon = 1e-8;
% min concentration below which flocculation processes are not calculated
f_clim=0.001

dt = 1.0
tstart = 0.0
tend   = 750.0 * 60.0

% size classes
f_diam = 1e-6 * ...
   [4.0, 6.1, 9.3, 14.2, 21.8, 33.2, 50.7, 77.5, 118.3, 180.6, 275.8, 421.2, 643.2, 982.3, 1500.0]';

nv_mud = length(f_diam);
% initial concentrations
% TODO - cv_wat should (nv_mud,t)
cv_wat = zeros(nv_mud,1);
cv_wat(5)=0.093;

f_vol = (pi/6.0)*f_diam.^3;
f_rho = rhoref+(rhosp-rhoref)*(f_dp0./f_diam).^(3.0-f_nf);
f_mass = zeros(nv_mud+1,1);
f_mass(1:nv_mud) = f_vol.*(f_rho-rhoref);
f_mass(nv_mud+1) = f_mass(nv_mud)*2.0+1.0;
if (f_diam(1) == f_dp0)
   f_mass(1)=f_vol(1)*rhosp;
end

f_ws = grav*(f_rho-rhoref).*f_diam.^2.0/(18.*0.001);

fm_print_init
% dimension arrays
f_coll_prob_sh=zeros(nv_mud,nv_mud);
f_coll_prob_ds=zeros(nv_mud,nv_mud);
f_g1_sh = zeros(nv_mud,nv_mud,nv_mud);
f_g1_ds = zeros(nv_mud,nv_mud,nv_mud);
f_g3 = zeros(nv_mud,nv_mud);
f_l3 = zeros(nv_mud);
f_g4=zeros(nv_mud,nv_mud,nv_mud);
f_l4=zeros(nv_mud,nv_mud);

% floc kernals
fm_flocmod_aggregation_statistics
fm_aggregation_gain
fm_shear_frag_gain
fm_aggregation_loss
fm_shear_frag_loss

fm_kernal_stats

% I think we can delete fmass(nv_mud+1) now
f_mass = f_mass(1:nv_mud,1);

fid = fopen('fm.dat','w');
t = tstart;
nt = 0;
while (t<tend)
   nt=nt+1;
   dtmin=dt;
   
   f_dt=dt;
   dttemp=0.0;
   
   cv_tmp=cv_wat; % concentration of all mud classes in one grid cell
   cvtotmud=sum(cv_tmp);
   % TODO - fix calculation of G
   fm_Gval
   f_gval = Gval;
   if( mod(nt,600) == 0 )
      fprintf(1,'t, G, cvtotmud: %f %f %f\n',t,Gval,cvtotmud)
   end
   NNin=cv_tmp./f_mass;
   
   if( any (NNin<0.0) )
      fprintf(1,'***************************************\n')
      fprintf(1,'CAUTION, negative mass at t = %f\n', t)
      fprintf(1,'***************************************\n')
   end
   
   if (cvtotmud > f_clim)
      while (dttemp <= dt)
         %     print*, 'f_dt:',f_dt
         
         fm_comp_fsd % NNin -> NNout
         % fm_mass_control
         ineg = find(NNout<0.0);
         mneg = sum( -NNout(ineg).*f_mass(ineg) );
         
         %     fprintf(1, 'mneg',mneg
         if (mneg > f_mneg_param)
            while (mneg > f_mneg_param)
               f_dt=min(f_dt/2.0,dt-dttemp);
               fm_comp_fsd % NNin -> NNout
               % fm_mass_control
               ineg = find(NNout<0.0);
               mneg = sum( -NNout(ineg).*f_mass(ineg) );
            end
            %         if (f_dt<1.0)
            %           fprintf(1, 'apres : Gval,f_dt',Gval, f_dt,dttemp
            %  end
         else
            
            if (f_dt<dt)
               while (mneg <f_mneg_param)
                  
                  if (dttemp+f_dt == dt)
                     fm_comp_fsd % NNin -> NNout
                     break
                  else
                     dt1=f_dt;
                     f_dt=min(2.0*f_dt,dt-dttemp);
                     fm_comp_fsd % NNin -> NNout
                     % fm_mass_control
                     ineg = find(NNout<0.0);
                     mneg = sum( -NNout(ineg).*f_mass(ineg) );
                     if (mneg > f_mneg_param)
                        f_dt=dt1;
                        fm_comp_fsd % NNin -> NNout
                        break
                     end
                  end
               end
            end
         end
         dtmin = min(dtmin,f_dt);
         dttemp = dttemp+f_dt;
         NNin = NNout; % update new Floc size distribution
         
         fm_mass_distribute % redistribute negative masses in NNin (if any) over positive classes,
         % depends on f_mneg_param
         
         if (abs(sum(NNin.*f_mass)-cvtotmud)>epsilon*100.0)
            fprintf(1, 'CAUTION flocculation routine not conservative!\n')
            fprintf(1, 'time = %f\n',t);
            fprintf(1, 'f_dt= %f\n',f_dt);
            fprintf(1, 'before : cvtotmud= %f\n',cvtotmud);
            fprintf(1, 'after  : cvtotmud= %f\n',sum(NNin.*f_mass));
            fprintf(1, 'absolute difference  : cvtotmud= %f\n',abs(cvtotmud-sum(NNin.*f_mass)));
            fprintf(1, 'absolute difference reference  : espilon= %f\n',epsilon);
            fprintf(1, 'before redistribution %f\n', sum(NNout.*f_mass));
            fprintf(1, 'after redistribution %f\n', sum(NNin.*f_mass));
            error('Simultation stopped')
         end
         
         if (dttemp == dt) break; end
      end % loop on full dt
   end % only if cvtotmud > f_clim
   
   if (abs( sum( NNin.*f_mass )-cvtotmud) > epsilon*10.0)
      fprintf(1, 'CAUTION flocculation routine not conservative!\n');
      fprintf(1, 'time = %g\n',t);
      fprintf(1, 'before : cvtotmud= %f\n',cvtotmud)
      fprintf(1, 'after  : cvtotmud= %f\n',sum( NNin.*f_mass ))
      fprintf(1, 'absolute difference  : cvtotmud= %f\n',...
         abs(cvtotmud-sum( NNin.*f_mass )))
      fprintf(1, 'absolute difference reference  : espilon= %f\n',epsilon);
      error('Simultation stopped')
   end
   
   % update mass concentration for all mud classes
   cv_wat  = NNin.*f_mass;
   
   % compute floc distribution statistics before output
   f_csum=0.0;
   f_ld50=1;
   f_ld10=1;
   f_ld90=1;
   
   f_davg = sum(NNin.*f_mass.*f_diam)./(sum(NNin.*f_mass)+eps);
   f_dtmin = dtmin;
   
   for iv1=1:nv_mud
      f_csum=f_csum + NNin(iv1)*f_mass(iv1)/((sum(NNin.*f_mass))+eps);
      if (f_csum > 0.1 && f_ld10)
         f_d10 =f_diam(iv1);
         f_ld10 = 0;
      end
      
      if (f_csum > 0.5 && f_ld50)
         f_d50 = f_diam(iv1);
         f_ld50=0;
      end
      
      if (f_csum > 0.9 && f_ld90)
         f_d90=f_diam(iv1);
         f_ld90=0;
      end
   end
   
   fprintf(fid,'%f %f %f\n',t, Gval, f_d50*1e6);
   t = t+dt;
end
fclose(fid)
fprintf(1,'END flocmod_main\n')
