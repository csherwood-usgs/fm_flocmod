% fm_print_init - Print out initial values
  fprintf(1,'\n');
  fprintf(1,'***********************\n')
  fprintf(1,'    FLOCMOD\n')
  fprintf(1,'***********************\n')
  fprintf(1,'class  diameter (um)  volume (m3)  density (kg/m3)  mass (kg) Ws (m/s)\n')
  for iv=1:nv_mud
     fprintf(1,'% 3d % 7.1f % 14g % 6.1f % 14g % 12f\n',...
        iv,f_diam(iv)*1e6,f_vol(iv),f_rho(iv),f_mass(iv),f_ws(iv))
  end
  fprintf(1,'\n')
  fprintf(1,' *** PARAMETERS ***\n')
  fprintf(1,'Primary particle size (f_dp0)                                : %f\n',f_dp0)
  fprintf(1,'Fractal dimension (f_nf)                                     : %f\n',f_nf)
  fprintf(1,'Flocculation efficiency (f_alpha)                            : %f\n',f_alpha)
  fprintf(1,'Floc break up parameter (f_beta)                             : %f\n',f_beta)
  fprintf(1,'Nb of fragments (f_nb_frag)                                  : %f\n',f_nb_frag)
  fprintf(1,'Ternary fragmentation (f_ater)                               : %f\n',f_ater)
  fprintf(1,'Floc erosion (pct of mass) (f_ero_frac)                      : %f\n',f_ero_frac)
  fprintf(1,'Nb of fragments by erosion (f_ero_nbfrag)                    : %f\n',f_ero_nbfrag)
  fprintf(1,'fragment class (f_ero_iv)                                    : %f\n',f_ero_iv)
  fprintf(1,'negative mass tolerated before redistribution (f_mneg_param) : %f\n',f_mneg_param)
  fprintf(1,'Boolean for differential settling aggregation (L_ADS)        : %d\n',l_ADS)
  fprintf(1,'Boolean for shear aggregation (L_ASH)                        : %d\n',l_ASH)
  fprintf(1,'Boolean for collision fragmenation (L_COLLFRAG)              : %d\n',l_COLLFRAG)
  fprintf(1,'Collision fragmentation parameter (f_collfragparam)          : %f\n',f_collfragparam)
  fprintf(1,'\n')
  fprintf(1,'*** END FLOCMOD INIT *** \n')    

  if ~(l_ADS+l_ASH)
     fprintf(1,'CAUTION : incompatible flocculation kernel options : \n')
     fprintf(1,'*****************************************************\n')
     fprintf(1,'l_ADS=%d\n',l_ADS)
     fprintf(1,'l_ASH=%d\n',l_ASH)
     error('simulation stopped')
  end
