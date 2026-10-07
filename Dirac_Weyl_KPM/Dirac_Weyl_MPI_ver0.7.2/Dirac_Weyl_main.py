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
    self.base_seed = int(params["Disorder_averaging"].get("base_seed", 271828))
    # Number of disorder realizations is Ndis = "bin_size" * "runs"

    if self.bin_size < 1 or self.runs < 1:
      raise ValueError("Both 'bin_size' and 'runs' must be positive integers.")
    
    if self.base_seed < 0:
      raise ValueError("'base_seed' must be a non-negative integer.")
    
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
      self.omega = np.linspace(-self.enerange, self.enerange, self.enegrid)
      self.max_retries = int(params["KPM_ADOS_details"].get("max_retries", 10))

      if self.max_retries < 0:
        raise ValueError("'max_retries' must be non-negative.")

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
  
  ################################
  # Construct one clean Hamiltonian and record its lattice dimensions.
  ################################
  def _build_clean_hamiltonian(self, Hamiltonians): # The "_" in front of the method name means this method is mainly for internal use. It does not have actual impact.
    # Also, "Hamiltonians = Dirac_Weyl_Hamiltonians()" is a class. So this method and a few methods below take a class as an input
    if self.Hidden_order and self.dimension != 3:
      raise ValueError("Currently, we do not support Hidden order outside 3D!")

    if self.dimension == 2:
      if self.spin != 1:
        raise ValueError("In 2D, only spin 1/2 is accepted!")
      if self.Hidden_order == True:
        raise ValueError("In 2D, we don't support hidden-order yet!")
      clean_Ham = Hamiltonians.Dirac_2D(self.linear_size, self.chem_pot, self.boundary_condition)
      orbitals_per_site = 2

    elif self.dimension == 3:
      if (self.Hidden_order == True) and (self.spin != 1):
        raise ValueError("Set spin = 1 before switching on Hidden_order.")
      
      if (self.spin == 1) and (self.Hidden_order == False):
        clean_Ham = Hamiltonians.Dirac_1_over_2(self.linear_size, self.chem_pot, self.boundary_condition)
        orbitals_per_site = 2
      elif (self.spin == 1) and (self.Hidden_order == True):
        clean_Ham = Hamiltonians.Dirac_3D_Hidden_order(self.linear_size, self.chem_pot, self.alpha_Hidden, self.boundary_condition)
        orbitals_per_site = 4
      elif self.spin == 3 and (self.Hidden_order == False):
        clean_Ham = Hamiltonians.Dirac_3_over_2(self.linear_size, self.chem_pot, self.boundary_condition)
        orbitals_per_site = 4
      elif self.spin == 5 and (self.Hidden_order == False):
        clean_Ham = Hamiltonians.Dirac_5_over_2(self.linear_size, self.chem_pot, self.boundary_condition)
        orbitals_per_site = 6
      else:
        raise ValueError("Please specify the correct spin number!")
    elif self.dimension in (4, 5, 6):
      if self.spin != 1:
        raise ValueError("In dimensions 4--6, only spin 1/2 is accepted!")
      if self.Hidden_order == True:
        raise ValueError("In dimensions 4--6, we don't support hidden-order yet!")
      
      # A dictionary that connects dimension number (d > 3) and the methods that generate Hamiltonians
      # We can do this shortcut only because 1. we only consider spin 1/2 in d > 3 (not spin 3/2, 5/2, etc.) and 2. no symmetry-breaking orders enter in d > 3
      builders = {4: Hamiltonians.Dirac_4D,
                  5: Hamiltonians.Dirac_5D,
                  6: Hamiltonians.Dirac_6D}
      
      clean_Ham = builders[self.dimension](self.linear_size, self.chem_pot, self.boundary_condition)
      orbitals_per_site = 4 if self.dimension in (4, 5) else 8
    else:
      raise ValueError("Please specify the correct dimension number!")

    self.size = clean_Ham.shape[0]
    self.number_of_sites = int(self.size / orbitals_per_site)
    
    return clean_Ham

  ###########################
  # Return deterministic disorder RNG and Kwant seed for one attempt.
  ###########################
  def _random_streams(self, nrun, nbin, attempt):
    ##############################################################
    # For MPI-based random processes (both random disordered potential and random vectors in KPM), we first create a class that mixes sources of entropy 
    # in a reproducible way. In other words, we first processes the vanilla interger seeds (here provided by "self.base_seed") into high-quality, 
    # randomized initial states. 
    #
    # Then, from those classes, we obtain the random-quantity generators as "disorder_rng" for the disordered potentials and "kpm_seed" for random vectors in KPM. 
    ##############################################################
    disorder_sequence = np.random.SeedSequence([self.base_seed, nrun, nbin, attempt, 0]) 
    kpm_sequence = np.random.SeedSequence([self.base_seed, nrun, nbin, attempt, 1])
    disorder_rng = np.random.default_rng(disorder_sequence) # A random number generator for disordered potentials (NOT the disorder itself!)
    kpm_seed = int(kpm_sequence.generate_state(1, dtype=np.uint32)[0])
    
    return disorder_rng, kpm_seed

  ###########################
  # Calculate one complete disorder bin; this is the MPI work unit.
  #
  # Here we arrange the main functionalities (such as KPM-ADOS) in this method, which will be called in every one iteration on disorder realizations. 
  ###########################
  def _calculate_one_run(self, nrun, output_file, Hamiltonians):
    if self.KPM_ADOS == True:
      bin_ADOS = np.zeros((self.enegrid, self.Nm_power_diff))
      bin_time_per_incre = np.zeros(self.Nm_power_diff) # a list of float points, they measure the time consumed per increment of Nm
    else:
      bin_ADOS = None
      bin_time_per_incre = None

    messages = []
    nbin = 0; attempt = 0 # "nbin" is capped at self.bin_size, while "attempt" is (or should be) capped at "self.max_retries"

    while nbin < self.bin_size:
      # We've allowed the program to retry several times if a severe error (i.e., big negative number in spectral weight) pops up. 
      # If the times of "attempt"s exceed the "self.max_retries", then something is fundamentally wrong with KPM and/or the lattice.
      if attempt > self.max_retries:
        raise RuntimeError("nrun = %d, nbin = %d exceeded max_retries = %d. " % (nrun, nbin, self.max_retries))

      Ndis = nrun * self.bin_size + nbin # the index of disorder realization (in general)

      start_time_Ndis = time.time()

      clean_Hamiltonian = self._build_clean_hamiltonian(Hamiltonians)
      disorder_rng, self.kpm_seed = self._random_streams(nrun, nbin, attempt)
      
      self.disorder_attempt = attempt

      if self.rank == 0:
        if (nbin == 0) and (nrun == 0) and (attempt == 0):
          output_file.write("The total number of sites = %d \n" % self.number_of_sites)
          output_file.write("The Hamiltonian dimension = %d \n\n" % self.size)
          output_file.write("###########################################################\n\n")
          # "output_file" should be equal to "None" for all other ranks except rank-0, for which it actually refers to the output file. 

      # Obtain the disorder potential (a diagonal matrix having the same dimension as that of the clean Hamiltonian)
      disorder_potential = Hamiltonians.make_random(W = self.dis_strength, dimension = self.dimension, Hidden_order = self.Hidden_order, 
                                                    distribution = self.distribution, size = self.number_of_sites, spin = self.spin, rng = disorder_rng)
      
      # Create the full Hamiltonian (as the clean part added to the disordered part)
      full_Hamiltonian = clean_Hamiltonian.tocsr() + disorder_potential
      
      del clean_Hamiltonian, disorder_potential
      gc.collect()

      severe_error = False # The severe error tag is True if self.KPM_error_tag == 2 
                           # There might be some severe error coming from other routines (maybe Lanczos or so), we don't have to worry about it now, though.
      if self.KPM_ADOS == True:
        self.KPM_calculations.KPM_Anderson_ADOS(self, final_system = full_Hamiltonian)
        severe_error = (self.KPM_error_tag == 2)

        if severe_error == True:
          messages.append("## Severe KPM error: Negative definite appears for Ndis = %d, attempt = %d, Nm = %s; "
                          "DOS written to %s. Retrying with a new deterministic stream." % (Ndis, attempt, self.moment_error, self.name_error_DOS))
        
        else: # we accumulate data as no severe error shows up.
          if self.KPM_error_tag == 1:
            messages.append("## Minor KPM error: negative-infinitesimal appears for Ndis = %d "
                            "at Nm = %s; rectified." % (Ndis, self.moment_error))
          
          bin_ADOS += self.DOS_by_Nm # The accumulated ADOS by the end of one particular Ndis-th disorder realization is done. 
          bin_time_per_incre += self.time_per_incre # to accumulate the time consumed for this disorder realization
          
          self.tot_ADOS_bins = bin_ADOS
          self.data_output.write_temp_data(nrun, nbin, self)
          
          del self.DOS_by_Nm
      
      ##############################################
      # If possible, some other routines carries for each disorder realization (such as Lanczos algorithm or ED or something else ...)
      ##############################################

      del full_Hamiltonian
      gc.collect()

      if severe_error == True:
        attempt += 1 # Note that we do NOT accumulate any data when "severe_error" is triggered in this realization. Everything is thrown out of window then.
        continue
      else:
        # real-time update of writing messages. This is NOT written in the output file "YT_....out",
        # since it's strongly discouraged to ask multiple MPI ranks to open and write in the same file.
        # Instead, we let the real-time update be "printed out" in the standard output (if there exists one). 
        message = ("The data in Ndis = %d saved by MPI rank %d. Time consumption: %s s" % (Ndis, self.rank, time.time() - start_time_Ndis))
        messages.append(message)
        print(message, flush = True)

        nbin += 1; attempt = 0
    # End of while nbin < self.bin_size
    
    return bin_ADOS, bin_time_per_incre, {"nrun": nrun, "rank": self.rank, "size": self.size, "number_of_sites": self.number_of_sites, "messages": messages}
    # The last one is a summary dictionary.     

  ###################################################
  # This method writes the header in the output file titled "YT_... .out"
  # Note, this method is called only by the 0th worker (i.e., when "self.rank == 0")
  ###################################################
  def _write_header(self, output_file, mpi_size):
    output_file.write("Datafiles are named with: " + str(self.label) + "\n\n")

    output_file.write("Dimension of the lattice: " + str(self.dimension) + "D\n")
    output_file.write("Spin of the fermions = " + str(self.spin) + "/2\n")
    output_file.write("Boundary condition: " + str(self.boundary_condition) + "\n")
    output_file.write("Linear size, L = " + str(self.linear_size) + "\n\n")
    
    output_file.write("Disorder strength, W = " + str(self.dis_strength) + "\n")
    output_file.write("Disorder distribution = " + self.distribution + "\n")
    output_file.write("Chemical potential, mu = " + str(self.chem_pot) + "\n")
    output_file.write("Base random seed = " + str(self.base_seed) + "\n\n")
    
    output_file.write("Hidden-order conditions: " + str(self.Hidden_order) + "\n")
    if self.Hidden_order == True:
      output_file.write("  Hidden-order coupling constant, alpha = " + str(self.alpha_Hidden) + "\n\n")
    else:
      output_file.write("\n")
    
    output_file.write("KPM-ADOS calculation? " + str(self.KPM_ADOS) + "\n")
    if self.KPM_ADOS == True:
      output_file.write("  Energy range = " + str(self.enerange) + "\n")
      output_file.write("  Energy grid = " + str(self.enegrid) + "\n")
      output_file.write("  Number of initial moments = " + str(self.Nm_init) + "\n")
      output_file.write("  Number of final moments = " + str(self.Nm_fin) + "\n")
      output_file.write("  Number of random vectors = " + str(self.random_vecs) + "\n")
      output_file.write("  All the numbers of moments, Nm = " + str(self.list_of_Nm) + "\n\n")

      if (self.random_vecs > 20) or (self.random_vecs < 9):
        output_file.write(" ##  WARNING: random_vecs is outside the recommended range 9--20.\n\n")
    else:
      output_file.write("\n")
    
    output_file.write("Disorder configurations to average = %d (Every %d of them is in a bin)\n\n" % (self.runs * self.bin_size, self.bin_size))
    
    output_file.write("MPI ranks, mpi_size = %d \n" % mpi_size)
    if mpi_size == 1:
      output_file.write("This job is effectively a one-CPU serial job \n\n")
    else:
      output_file.write("Work assignment: rank r calculates nrun = r, r + MPI_size, r + 2 * MPI_size ...\n\n")

  ##############################################
  # The main program of this code. 
  ##############################################
  def Main(self):
    # Distribute independent nrun bins over MPI ranks and reduce the ADOS.
    from mpi4py import MPI
    from Hamiltonians import Dirac_Weyl_Hamiltonians
    from postroutines import Postroutines

    comm = MPI.COMM_WORLD
    self.rank = comm.Get_rank() # The rank of each processor, an interger. Different processors handling this script has different "self.rank".
    mpi_size = comm.Get_size()  # The number of processors participating this job. 

    start_time = time.time()

    output_file = None
    
    # It is only the rank-0 processor that writes the output file. The barrier prevents other
    # ranks (i.e., other processors) from writing their distinct bin directories during that cleanup.
    # It is dangerous to allow multiple MPI ranks to write on the same file, that may causes interleave or corrupt lines. 
    if self.rank == 0:
      name_out = os.path.join(self.outdir,"YT_" + self.label + "_W" + str(self.dis_strength).replace(".", "") + ".out")
      output_file = open(name_out, "w", buffering = 1)
      self._write_header(output_file, mpi_size)
    
      output_file.write("The job begins at %s\n" % time.asctime(time.localtime(start_time)))
      output_file.write("MPI calculation is in progress.\n")
      output_file.write("###########################################################\n\n")
    # comm.Barrier()
    
    Hamiltonians = Dirac_Weyl_Hamiltonians()
    
    self.data_output = Postroutines()
    
    if self.KPM_ADOS == True:
      from core_routines import Core_Routines
      self.KPM_calculations = Core_Routines()
      self.tot_ADOS_bins = np.zeros((self.enegrid, self.Nm_power_diff))
      self.tot_ADOS_runs = np.zeros((self.enegrid, self.Nm_power_diff)) # Kept for compatibility with the output helper.

    # It is only the rank-0 that clears shared temporary output.
    if self.rank == 0:
      self.data_output.reset_temp_data(self)
    comm.Barrier() # Ask those faster-running ranks to wait for the slower ranks. 
    
    # The quantities with prefix "local_" are used locally within each processor
    if self.KPM_ADOS == True:
      local_ADOS = np.zeros((self.enegrid, self.Nm_power_diff))
      local_time_per_incre = np.zeros(self.Nm_power_diff)
    local_summaries = []

    #################################################
    # Embarrassingly parallel region: no communication occurs inside this loop.
    #################################################
    for nrun in range(self.rank, self.runs, mpi_size):
      bin_ADOS, bin_time, summary = self._calculate_one_run(nrun, output_file, Hamiltonians)
      if self.KPM_ADOS == True:
        local_ADOS += bin_ADOS
        local_time_per_incre += bin_time
      local_summaries.append(summary)
      print("MPI rank %d completed nrun= %d" % (self.rank, nrun), flush=True)

    gathered_summaries = comm.gather(local_summaries, root = 0)

    # Reduce data from other ranks to rank 0.
    if self.KPM_ADOS == True:
      total_ADOS = np.zeros_like(local_ADOS) if self.rank == 0 else None
      total_time_per_incre = (np.zeros_like(local_time_per_incre) if self.rank == 0 else None)
      comm.Reduce(local_ADOS, total_ADOS, op=MPI.SUM, root=0)
      comm.Reduce(local_time_per_incre, total_time_per_incre, op=MPI.SUM, root=0)

    if self.rank != 0:
      return None

    summaries = sorted([item for rank_items in gathered_summaries for item in rank_items], 
                       key=lambda item: item["nrun"])
    if len(summaries) != self.runs:
      raise RuntimeError("MPI completed %d of %d nrun bins." % (len(summaries), self.runs))

    # Write the arranged version of messages from every single disorder realization into the output file "YT....out"
    for summary in summaries:
      for message in summary["messages"]:
        output_file.write(message + "\n")

    if self.KPM_ADOS == True:
      # Do the disorder average (the total DOS divided by the number of disorder realizations)
      self.tot_ADOS_runs = total_ADOS # The total ADOS (not yet divided by the number of disorder realizations)
      self.tot_ADOS = total_ADOS / (self.runs * self.bin_size) # The ADOS
      
      self.data_output.write_DOS(self)
      output_file.write("\nThe total ADOS successfully obtained.\n")

    end_time = time.time()
    output_file.write("###########################################################\n")
    output_file.write("\nJob Complete!\n")
    output_file.write("\nMPI wall-clock time: %s s\n\n" % (end_time - start_time))
    output_file.write("...within which: \n")
    if self.KPM_ADOS:
      output_file.write("  Accumulated KPM-ADOS time across all realizations:\n")
      for i in range(self.Nm_power_diff):
        output_file.write(
          "    Increment up to Nm = %s costed %s s\n"
          % (self.list_of_Nm[i], total_time_per_incre[i]))
    output_file.write("\nJob ends at %s\n" % time.asctime(time.localtime()))
    output_file.close()

    return self.tot_ADOS if self.KPM_ADOS else None

############################################################################################
############################################################################################

def readin(input_file_name):
  import yaml

  with open(input_file_name, "r") as input_file:
    params = yaml.safe_load(input_file)

  return params

if __name__ == "__main__":
  from mpi4py import MPI
  import traceback

  input_file_name = sys.argv[1] if len(sys.argv) > 1 else "input.yaml"

  try:
    params = readin(input_file_name)
    simulation = Disordered_Dirac_Weyl(params)
    simulation.Main()
  except Exception:
    print("Unhandled exception on MPI rank %d:" % MPI.COMM_WORLD.Get_rank(),
          file=sys.stderr, flush=True)
    traceback.print_exc()

    # An immediate coordinated stop prevents other ranks from hanging inside
    # a gather, Reduce, or Barrier after one rank has failed.
    MPI.COMM_WORLD.Abort(1)

