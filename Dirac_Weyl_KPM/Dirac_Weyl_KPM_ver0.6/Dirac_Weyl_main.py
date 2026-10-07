import numpy as np
import scipy
import math

import kwant
import tinyarray

import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)

import time # monitoring the time
import os # better manipulating the output files
import gc # garbage collector

import sys

class Disordered_Dirac_Weyl:
  ###############################################################
  # This is the main routine the takes care of level statistics, ADOS, and IPR calculations.
  #
  #                                                                 --Yongtai Li
  ###############################################################

  def __init__(self, params):
    # data_and_output
    self.label = params["Data_and_output"]["label"]

    # Lattice
    self.dimension = int(params["Lattice"]["dimension"])
    self.spin = int(params["Lattice"]["spin"])
    self.boundary_condition = str(params["Lattice"]["bound_cond"])
    self.linear_size = int(params["Lattice"]["linear_size"])
    self.chem_pot = np.float64(params["Lattice"]["mu"])

    # Disorder
    self.dis_strength = np.float64(params["Disorder"]["W"])
    self.distribution = params["Disorder"]["distr"]

    # Disorder_averaging
    self.bin_size = int(params["Disorder_averaging"]["bin_size"])
    self.runs = int(params["Disorder_averaging"]["runs"])
    # Number of disorder realizations is Ndis = "bin_size" * "runs"

    # Symmetry_breaking_terms
    self.Hidden_order = params["Symmetry_breaking_terms"]["Hidden_order"] # a boolean value
    if (self.Hidden_order == True):
      self.alpha_Hidden = np.float64(params["Symmetry_breaking_terms"]["alpha_Hidden"])
    elif (self.Hidden_order == False):
      pass
    else:
      raise ValueError("Please specify Hidden_order, whether it's 'True' or 'False'.")
    
    # KPM_ADOS
    self.KPM_ADOS = params["KPM_ADOS_details"]["KPM_ADOS"] # a boolean value
    if (self.KPM_ADOS == True):  
      self.Nm_init = int(params["KPM_ADOS_details"]["Nm_init"])
      self.Nm_fin = int(params["KPM_ADOS_details"]["Nm_fin"])
      self.random_vecs = int(params["KPM_ADOS_details"]["random_vecs"])

      # the starting/ending number of moments must be a POWER of BASE 2 and interger exponent
      if (float(math.log(self.Nm_init, 2)).is_integer()==False) or (float(math.log(self.Nm_fin, 2)).is_integer()==False):
        raise ValueError("The initial/final numbers of moments must be a POWER of BASE 2 and interger exponent.")
      if self.Nm_init > self.Nm_fin:
        raise ValueError("The initial number of moments canNOT be GREATER than the final one.")

      self.Nm_power_diff = int(math.log(self.Nm_fin, 2) - math.log(self.Nm_init, 2) + 1) # Difference in the exponents of Nm

      self.enegrid = int(params["KPM_ADOS_details"]["enegrid"])
      self.enerange = np.float64(params["KPM_ADOS_details"]["enerange"])

      # a list of all the number of moments considered in our calculation (mostly for output writing purposes)
      list_of_Nm = []
      for Nm_power_red in range(self.Nm_power_diff):
        if Nm_power_red == 0:
          list_of_Nm.append(self.Nm_init)
        else:
          last_element = list_of_Nm[Nm_power_red - 1]
          list_of_Nm.append(last_element + 2**int(math.log(last_element, 2)))
      self.list_of_Nm = list_of_Nm

    elif (self.KPM_ADOS == False):
      pass
    else:
      raise ValueError("Please specify KPM_ADOS, whether it's 'True' or 'False'.")

    # output file writting
    self.outdir = "datafiles"
    os.makedirs(self.outdir, exist_ok = True)

    # temporary file writting
    self.tempdir = "tempfiles"
    os.makedirs(self.tempdir, exist_ok = True)

    # erronous datafile writting
    self.errordir = "errorfiles"
    os.makedirs(self.errordir, exist_ok = True)
  
  def Main(self):
    start_time = time.time()

    #############################################
    # Writing output file in this local machine
    #############################################
    #name_out = str("YT" + "_W" + str(self.dis_strength).replace(".","") + ".out")
    name_out = os.path.join(self.outdir, "YT_" + self.label + "_W" + str(self.dis_strength).replace(".", "") + ".out")
    output_file = open(name_out, "w", buffering = 1)

    output_file.write("Datafiles are named with: " + str(self.label) + "\n\n")

    output_file.write("Dimension of the lattice: " +str(self.dimension) + "D" + "\n")
    output_file.write("Spin of the fermions = " + str(self.spin) + "/2 " + "\n")
    output_file.write("Boundary condition: " + str(self.boundary_condition) + "\n")
    output_file.write("Linear size, L = " + str(self.linear_size) + "\n\n")

    output_file.write("Disorder strength, W = " + str(self.dis_strength) + "\n")
    output_file.write("Disorder distribution = " + self.distribution + "\n")
    output_file.write("Chemical potential, mu = " + str(self.chem_pot) + "\n\n")

    output_file.write("Hidden-order conditions: " + str(self.Hidden_order) + "\n")
    if (self.Hidden_order == True):
      output_file.write("Hidden-order coupling constant, alpha = " + str(self.alpha_Hidden) + "\n\n")

    output_file.write("KPM-ADOS calculation? " + str(self.KPM_ADOS) + "\n")
    if (self.KPM_ADOS == True):
      output_file.write("  Energy range = " + str(self.enerange) + "\n")
      output_file.write("  Energy grid = " + str(self.enegrid) + "\n\n")

      output_file.write("  Number of initial moments = " + str(self.Nm_init) + "\n")
      output_file.write("  Number of final moments = " + str(self.Nm_fin) + "\n")
      
      output_file.write("  Number of random vectors = " + str(self.random_vecs) + "\n")
      if (self.random_vecs > 20) or (self.random_vecs < 9):
        output_file.write("\n" + "  #########################################" + "\n")
        output_file.write("  # " + "\n")
        output_file.write("  # WARNING: Your number of random vectors is out of the recommended range. (9-20)" + "\n")
        output_file.write("  # Take your results with a pinch of salt!." + "\n")
        output_file.write("  # " + "\n")
        output_file.write("  #########################################" + "\n\n")

      output_file.write("  All the numbers of moments, Nm = " + str(self.list_of_Nm) + "\n\n")
      total_time_per_incre = np.zeros(len(self.list_of_Nm))
    else:
      pass

    output_file.write("Disorder configurations to average = "+str(self.runs * self.bin_size) + \
            " (Every %s of them is in a bin)"%self.bin_size +"\n\n")

    ###################################
    #  Initialization
    ###################################
    from Hamiltonians import Dirac_Weyl_Hamiltonians
    Hamiltonians = Dirac_Weyl_Hamiltonians()

    from postroutines import Postroutines # The class "Postroutines" is used as a helper class
    self.data_output = Postroutines()     # which contains several methods to be used below

    if (self.KPM_ADOS == True):
      # prep for ADOS
      self.tot_ADOS_bins = np.zeros((self.enegrid, self.Nm_power_diff))
      self.tot_ADOS_runs = np.zeros((self.enegrid, self.Nm_power_diff))

      from core_routines import Core_Routines
      self.KPM_calculations = Core_Routines()

    #####################################
    # The main loop
    #####################################
    nrun = 0
    while nrun < self.runs: # Let "nrun" be the index of each single bin
      nbin = 0
      while nbin < self.bin_size: # Let "nbin" be the index of disorder realization in a single bin
        Ndis = nrun * self.bin_size + nbin # The index of disorder realization, from 0 to (runs*bin_size - 1)

        ###################################
        # Initializing clean and full Hamiltonians in every realization.
        # The clean Hamiltonian will not be retained throughout the entire calculation, 
        # it will be deleted right after having the full Hamiltonian for the sake of memory saving.
        ###################################
        if (self.Hidden_order == True) and (self.dimension != 3):
          raise ValueError("Currently, we do not support Hidden order on other dimensions except 3D !")

        if (self.dimension == 2): 
          if (self.spin == 1):
            self.clean_Hamiltonian = Hamiltonians.Dirac_2D(self.linear_size, self.chem_pot, self.boundary_condition)
            self.size = self.clean_Hamiltonian.shape[0]; self.number_of_sites = int(self.size / 2)
          else:
            raise ValueError("in 2D, other spins except 1/2 is not accepted!")

        elif (self.dimension == 3):
          if ((self.Hidden_order == True) and (self.spin != 1)):
            raise ValueError("Make sure to set spin = 1 before switching on Hidden_order = True")

          if (self.spin == 1):
            if (self.Hidden_order == False):
              self.clean_Hamiltonian = Hamiltonians.Dirac_1_over_2(self.linear_size, self.chem_pot, self.boundary_condition)
              self.size = self.clean_Hamiltonian.shape[0]; self.number_of_sites = int(self.size / 2)
            else:
              self.clean_Hamiltonian = Hamiltonians.Dirac_3D_Hidden_order(self.linear_size, self.chem_pot, self.alpha_Hidden, 
                                                                          self.boundary_condition) # there's an alternate version, "Dirac_3D_Hidden_order_ver2"
              self.size = self.clean_Hamiltonian.shape[0]; self.number_of_sites = int(self.size / 4) # with hidden order, the number of spin-orbital is 4
          elif (self.spin == 3):
            self.clean_Hamiltonian = Hamiltonians.Dirac_3_over_2(self.linear_size, self.chem_pot, self.boundary_condition)
            self.size = self.clean_Hamiltonian.shape[0]; self.number_of_sites = int(self.size / 4)
          elif (self.spin == 5):
            self.clean_Hamiltonian = Hamiltonians.Dirac_5_over_2(self.linear_size, self.chem_pot, self.boundary_condition)
            self.size = self.clean_Hamiltonian.shape[0]; self.number_of_sites = int(self.size / 6)
          else:
            raise ValueError("Please specify the correct spin number!")

        elif (self.dimension == 4):
          if (self.spin == 1):
            self.clean_Hamiltonian = Hamiltonians.Dirac_4D(self.linear_size, self.chem_pot, self.boundary_condition)
            self.size = self.clean_Hamiltonian.shape[0]; self.number_of_sites = int(self.size / 4) # in 4D systems, the number of spin orbitals is 4
          else:
            raise ValueError("in 4D, other spins except 1/2 is not accepted!")

        elif (self.dimension == 5):
          if (self.spin == 1):
            self.clean_Hamiltonian = Hamiltonians.Dirac_5D(self.linear_size, self.chem_pot, self.boundary_condition)
            self.size = self.clean_Hamiltonian.shape[0]; self.number_of_sites = int(self.size / 4) # in 5D systems, the number of spin orbitals is 4
          else:
            raise ValueError("in 5D, other spins except 1/2 is not accepted!")

        elif (self.dimension == 6):
          if (self.spin == 1):
            self.clean_Hamiltonian = Hamiltonians.Dirac_6D(self.linear_size, self.chem_pot, self.boundary_condition)
            self.size = self.clean_Hamiltonian.shape[0]; self.number_of_sites = int(self.size / 8) # in 6D systems, the number of spin orbitals is 8
          else:
            raise ValueError("in 6D, other spins except 1/2 is not accepted!")

        else:
          raise ValueError("Please specify the correct dimension number!")

        if (Ndis==0):
          output_file.write("The total number of sites = " + str(self.number_of_sites) + "\n")
          output_file.write("The Hamiltonian dimension = " + str(self.size) + "\n")
          output_file.write("The job begins at %s" %time.asctime(time.localtime()) + "\n")
          output_file.write("###########################################################" + "\n\n")

        start_time_Ndis = time.time()

        # Make full Hamiltonian
        disorder_potential = Hamiltonians.make_random(W = self.dis_strength, dimension = self.dimension, Hidden_order = self.Hidden_order,
                                                      distribution = self.distribution, size = self.number_of_sites, spin = self.spin)

        full_Hamiltonian = self.clean_Hamiltonian.copy()
        full_Hamiltonian = full_Hamiltonian.tocsr() + disorder_potential

        # delete the clean Hamiltonian and disorder potential to yield some memory
        del self.clean_Hamiltonian
        del disorder_potential
        gc.collect()

        self.Ndis = Ndis

        if (self.KPM_ADOS == True):
          self.KPM_calculations.KPM_Anderson_ADOS(self, final_system = full_Hamiltonian)
          
          ###############################
          # Check for the severe errors
          ###############################
          if self.KPM_error_tag == 2:
            output_file.write("\n" + "## Warning: SEVERE error" + "\n")
            output_file.write("## This occurs at Nm = " + str(self.moment_error) + "\n")
            output_file.write("## DOS printed as " + self.name_error_DOS + "\n")
            output_file.write("## We re-do this Ndis iteration." + "\n\n")
          else:
            if self.KPM_error_tag == 1:
              output_file.write("\n" + "## Warning: minor error: a negative infinitesimal appeared" + "\n")
              output_file.write("## This occurs when raising Nm up to: " + str(self.moment_error) + "\n")
              output_file.write("## The minor numerical error is rectified." + "\n")

            # ADOS: Save in the "runs" bin
            self.tot_ADOS_bins += self.DOS_by_Nm # the arrays "self.tot_ADOS_bins" and "self.DOS_by_Nm" have the same dimension

            # Save the temporary data
            self.data_output.write_temp_data(nrun, nbin, self)
            total_time_per_incre += self.time_per_incre # Accumulation of time consumption resolved in Nm

            # delete the DOS calculated for this particular disorder realization
            del self.DOS_by_Nm
            gc.collect()
        # end of "if (self.KPM_ADOS == "True")"
        
        ############
        # If possible, some other routines to do for each disorder realization
        ############
        if (self.KPM_error_tag != 2): # only when the error tag indicates no severe error can we safely proceed.
          end_time_Ndis = time.time()
          output_file.write("The data in Ndis = %s " %Ndis + "saved. Time consumption: %s s"
                              %(end_time_Ndis - start_time_Ndis)+"\n")
          nbin = nbin + 1
        
        # delete "full_Hamiltonian"
        del full_Hamiltonian
        gc.collect()
      # end of "while nbin"

      if (self.KPM_ADOS == True):
        self.tot_ADOS_runs += self.tot_ADOS_bins # The array "self.tot_ADOS_runs" is added by "self.tot_ADOS_bins" after each nrun
        self.tot_ADOS_bins[:, :] = 0.0e0 # initializing the "self.tot_ADOS_bins" after one iteration of nrun

      nrun = nrun + 1
    # end of "while nrun"

    if (self.KPM_ADOS == True):
      # finalizing ADOS
      self.tot_ADOS = self.tot_ADOS_runs / (self.runs * self.bin_size)
      output_file.write("\n" + "The total ADOS successfully obtained. \n")

      self.data_output.write_DOS(self)

    end_time = time.time()

    #######################################
    # Writing ending messages
    #######################################
    output_file.write("###########################################################" + "\n")
    output_file.write("\n" + "Job Complete! \n")
    output_file.write("\n" + "Total time consumption in this calculation: %s s" %(end_time - start_time) + "\n\n")
    output_file.write("within which, " + "\n\n")
    
    if (self.KPM_ADOS == True):
      output_file.write("  In KPM-ADOS calculations:")
      for i in range(self.Nm_power_diff):
        if (i==0):
          output_file.write("  An increment up to Nm = %s " %self.list_of_Nm[i] +
                          "costed %s s" %total_time_per_incre[i] + "\n")
        else:
          output_file.write("                             An increment up to Nm = %s " %self.list_of_Nm[i] + 
                          "costed %s s" %total_time_per_incre[i] + "\n")

    output_file.write("\n" + "Job ends at %s" %time.asctime(time.localtime()) + "\n")
    output_file.close()

    exit()

############################################################################################
############################################################################################

def readin(input_file_name):
  import yaml

  with open(input_file_name, "r") as input_file:
    params = yaml.safe_load(input_file)

  return params

if __name__ == "__main__":
  input_file_name = "input.yaml" # sys.argv[1]

  params = readin(input_file_name)
 
  simulation = Disordered_Dirac_Weyl(params)
  simulation.Main()

