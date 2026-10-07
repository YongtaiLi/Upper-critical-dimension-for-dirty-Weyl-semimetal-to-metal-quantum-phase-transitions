import os

class Postroutines:
  #################################################################################################################
  #
  #   This file contains functions that print out data and calculate some observables.
  #                                                                                           --Yongtai Li
  #
  #################################################################################################################

  def write_DOS(self, main):
    name_ADOS = os.path.join(main.outdir, "ADOS_" + main.label + "_W" + str(main.dis_strength).replace(".", "") + ".dat")

    with open(name_ADOS, "w") as ADOS_file:
      ADOS_file.write("## nwn")
      for Nm_power_red in range(main.tot_ADOS.shape[1]):
        ADOS_file.write("       Nm = " + str(main.list_of_Nm[Nm_power_red]))
      ADOS_file.write("\n")

      for i in range(len(main.omega)): # Rows are for frequencies
        ADOS_file.write(" {:11.4e}".format(main.omega[i]))
  
        for Nm_power_red in range(main.tot_ADOS.shape[1]): # columns are for Nm
          ADOS_file.write("  {:15.8e}".format(main.tot_ADOS[i, Nm_power_red]))
        ADOS_file.write("\n")

#############################################################################################
#############################################################################################

  def reset_temp_data(self, main):
    # Remove stale temporary results before any MPI rank starts working.
    import shutil

    for name in os.listdir(main.tempdir):
      path_tempfiles = os.path.join(main.tempdir, name)
      if os.path.isfile(path_tempfiles) or os.path.islink(path_tempfiles):
        os.remove(path_tempfiles)
      elif os.path.isdir(path_tempfiles):
        shutil.rmtree(path_tempfiles)

#############################################################################################
#############################################################################################

  def write_temp_data(self, nrun, nbin, main): # Write temporary datafiles within folder temp_datafiles/
    import numpy as np
    import shutil

    # Create (sub-) directories for each bin
    tempdir_per_bin = main.tempdir + "/bin_" + str("{:03d}".format(nrun))
    os.makedirs(tempdir_per_bin, exist_ok = True)

    # Each rank owns distinct nrun values.  Clear this bin only when its first
    # realization is saved, so later nbin files are not deleted.
    if nbin == 0:
      for name in os.listdir(tempdir_per_bin):
        path_in_bin = os.path.join(tempdir_per_bin, name)
        if os.path.isfile(path_in_bin) or os.path.islink(path_in_bin):
          os.remove(path_in_bin)
        elif os.path.isdir(path_in_bin):
          shutil.rmtree(path_in_bin)

    # Clearing up each bin before writing new datafiles in this bin.
    # Usually it wouldn't hurt to not do this (that is, to save the temporary files after every disorder realizations)
    # But it would look very messy if the "bin_size" is too large (since in that case there would be "bin_size" number of temporary files in a sub-directory).
    if main.bin_size > 10:
      for name in os.listdir(tempdir_per_bin):
        path_in_bin = os.path.join(tempdir_per_bin, name)
        if os.path.isfile(path_in_bin) or os.path.islink(path_in_bin):
          os.remove(path_in_bin)
        elif os.path.isdir(path_in_bin):
          shutil.rmtree(path_in_bin)
  
    # Write the temporary files within the "tempdir_per_bin"
    name_ADOS_temp = os.path.join(tempdir_per_bin, "ADOS_" + main.label + "_W" + str(main.dis_strength).replace(".", "") + \
            "_nrun" + str("{:03d}".format(nrun)) + "_nbin" + str("{:02d}".format(nbin)) +".dat")

    with open(name_ADOS_temp, "w") as ADOS_temp_file:
      ADOS_temp_file.write("## No. of realization %s completed in this bin and written below. They are not normalized/averaged." %(nbin + 1) + "\n")
      ADOS_temp_file.write("## nwn")
      for Nm_power_red in range(main.tot_ADOS_bins.shape[1]): # Number of columns, for all the Nm
        ADOS_temp_file.write("       Nm = " + str(main.list_of_Nm[Nm_power_red]))
      ADOS_temp_file.write("\n")

      for i in range(len(main.omega)): # Rows are for frequencies
        ADOS_temp_file.write(" {:11.4e}".format(main.omega[i]))

        for Nm_power_red in range(main.tot_ADOS_runs.shape[1]): # columns are for Nm
          ADOS_temp_file.write("  {:15.8e}".format(main.tot_ADOS_bins[i, Nm_power_red]))
        ADOS_temp_file.write("\n")
