% fm_Gval - Returns Gval(t) for test case
Gval=0.0;
if (t < 7201.0)
   Gval=1.0;
elseif (t < 8401.0)
   Gval=2.0;
elseif (t < 9601.0)
   Gval=3.0;
elseif (t < 10801.0)
   Gval=4.0;
elseif (t < 12601.0)
   Gval=12.0;
elseif (t < 13801.0)
   Gval=4.0;
elseif (t < 15001.0)
   Gval=3.0;
elseif (t < 16201.0)
   Gval=2.0;
elseif (t < 21601.0)
   Gval=1.0;
elseif (t < 25201.0)
   Gval=0.0;
elseif (t < 30601.0)
   Gval=1.0;
elseif (t < 31801.0)
   Gval=2.0;
elseif (t < 33001.0)
   Gval=3.0;
elseif (t < 34201.0)
   Gval=4.0;
elseif (t < 36001.0)
   Gval=12.0;
elseif (t < 37201.0)
   Gval=4.0;
elseif (t < 38401.0)
   Gval=3.0;
elseif (t < 39601.0)
   Gval=2.0;
elseif (t < 45001.0)
   Gval=1.0;
elseif (t < 48601.0)
   Gval=0.0;
elseif (t < 54001.0)
   Gval=1.0;
elseif (t < 55201.0)
   Gval=2.0;
elseif (t < 56401.0)
   Gval=3.0;
elseif (t < 57601.0)
   Gval=4.0;
else
   Gval=12.0;
end