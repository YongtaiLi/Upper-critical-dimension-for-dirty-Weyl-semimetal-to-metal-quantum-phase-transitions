import numpy as np
import scipy
import math

# Sometimes there's (non-lethal) "RuntimeWarning". If bothered by that, try this below:
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)

import time # monitoring the time
import os # better manipulating the output files
import gc # garbage collector

class Core_Routines:
  ################################################################
  # This class contains main routines used in this project. 
  #
  #                                             --Yongtai Li
  ################################################################

  def KPM_Anderson_ADOS(self, main, final_system): # The default kernel is the Jackson kernel
    import kwant
    import tinyarray

    time_per_incre = np.zeros(len(main.list_of_Nm)) # An array containing the time costed with each increment of Nm

    DOS_by_Nm = np.zeros((main.enegrid, main.Nm_power_diff)) 
    omega = np.linspace(-main.enerange, main.enerange, main.enegrid)

    for Nm_power in range(int(math.log(main.Nm_init, 2)), int(math.log(main.Nm_fin, 2)) + 1):
      time_DOS_start = time.time()

      if Nm_power == int(math.log(main.Nm_init, 2)):
        spectrum = kwant.kpm.SpectralDensity(final_system, num_moments = main.Nm_init,
                                             num_vectors = main.random_vecs, eps = 0.05,
                                             rng = main.kpm_seed)
      else:
        Nm_to_add = 2**(Nm_power - 1)
        spectrum.add_moments(num_moments = Nm_to_add)

      integrate = spectrum.integrate() # Note that this quantity is complex (though the imaginary part is very small)
      energies = spectrum(omega)
      DOS_Nm = np.array(energies.real)
      DOS_Nm = DOS_Nm / integrate.real

      ####################################################################################
      # If an energy value is out of the actual spectrum of a system, then the KPM-calculated DOS
      # will be "nan" there. But physically they should be zero instead of "nan".
      # So, we set the "nan" to, for numerical sake, an infinitesimal number.
      ####################################################################################
      DOS_Nm[np.isnan(DOS_Nm)] = 1.0E-10

      ###############################
      # Check for the severe errors
      ###############################
      main.KPM_error_tag = 0
      DOS_flat = np.ravel(DOS_Nm)

      for nwn in range(len(DOS_flat)):
        if DOS_flat[nwn] < 0.0E0:
          if DOS_flat[nwn] > -1.0E-8: # minor error, typically happens when Nm is too large
            DOS_flat[nwn] = abs(DOS_flat[nwn])

            main.KPM_error_tag = 1
            main.moment_error = spectrum.num_moments
          else: # Severe error, and one would have to redo this calculation
            # write this error DOS into a file
            name_error_DOS = os.path.join(main.errordir, "DOS_error" + "_W" + str(main.dis_strength).replace(".","") + \
                    "_Ndis" + str(main.Ndis) + "_attempt" + str(main.disorder_attempt) + ".dat")
            with open(name_error_DOS, "w") as error_DOS_file:
              for i in range(len(omega)):
                error_DOS_file.write(str("{:15.8e}".format(omega[i])) + "   " + str( "{:15.8e}".format(DOS_flat[i])) + "\n")

            main.KPM_error_tag = 2
            main.name_error_DOS = name_error_DOS
            main.moment_error = spectrum.num_moments
            break # break the "for nwn" loop, to the current iteration of "for Nm_power" loop

      if main.KPM_error_tag == 2:
        break # break the "for Nm_power" loop, and the entire disorder realization will be "thrown out of window"
      else:
        DOS_by_Nm[:, (Nm_power - int(math.log(main.Nm_init, 2)))] += DOS_Nm[:] # Load DOS_Nm to DOS_by_Nm

        time_DOS_end = time.time()
        time_per_incre[(Nm_power - int(math.log(main.Nm_init, 2)))] += (time_DOS_end - time_DOS_start)
    # end of "for Nm_power" 

    if main.KPM_error_tag == 2: # the worst case scenario
      return main.KPM_error_tag, main.name_error_DOS, main.moment_error
    elif main.KPM_error_tag == 1:
      main.DOS_by_Nm, main.omega = DOS_by_Nm, omega
      main.time_per_incre = time_per_incre
      return main.omega, main.DOS_by_Nm, main.KPM_error_tag, main.time_per_incre, main.moment_error
    else: # The best case scenario
      main.DOS_by_Nm, main.omega = DOS_by_Nm, omega
      main.time_per_incre = time_per_incre
      return main.omega, main.DOS_by_Nm, main.KPM_error_tag, main.time_per_incre
