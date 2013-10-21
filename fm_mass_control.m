% fm_mass_control - Compute mass in every class after flocculation and returns negative mass if any
mneg=0.0;
% TODO - Vectorize
for iv1=1:nv_mud
   if (NN(iv1)<0.0)
      mneg=mneg-NN(iv1)*f_mass(iv1);
   end
end