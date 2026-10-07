import numpy as np
import scipy

class Dirac_Weyl_Hamiltonians:
  #####################################################
  #
  # This class contains Hamiltonians for a variety of TB-based Dirac/Weyl systems, for 2D, 3D, 4D, 5D, and 6D systems, 
  # for (in 3D systems) spins-1/2, 3/2, 5/2, etc (functions/methods named "Dirac_1_over_2", "Dirac_3_over_2", and "Dirac_5_over_2", respectively)
  # for Dirac systems on 2D, 4D, 5D, and 6D hypercubic lattices (functions/methods named "Dirac_2D", "Dirac_4D", "Dirac_5D", and "Dirac_6D", respectively),
  # and (in 3D) Dirac systems with a symmetry-breaking term that induces a hidden order (functions/methods named "Dirac_3D_Hidden_order"). 
  #
  # The lattices may have either periodic boundary conditions (boundary_condition == "PBC") or open boundary conditions (boundary_condition == "OBC").
  #
  #                                                                             --Yongtai Li
  #####################################################

  def Dirac_1_over_2(self, L, mu, boundary_condition):
    from scipy.sparse import diags, kron, eye, csr_matrix
    import gc

    # Define the Pauli matrices sigma_0, sigma_x, sigma_y, and sigma_z
    sigma_0 = eye(2, format="csr") # 2-by-2 identity matrix
    sigma_x = csr_matrix(np.array([[0, 1], [1,  0]], dtype = np.complex128))
    sigma_y = csr_matrix(np.array([[0, -1j], [1j, 0]], dtype = np.complex128))
    sigma_z = csr_matrix(np.array([[1, 0], [0, -1]], dtype = np.complex128))

    # Define the bonds (i.e., Hamiltonians that resolved in real-space, in direction-wise)
    v = 1.0e0 # To reproduce Andras Szabo's results, set v = 0.50e0
    Hx = diags([(v / (2j)) * np.ones(L - 1), (-v / (2j)) * np.ones(L - 1)], [1, -1], format = "lil", dtype = np.complex128) # [1 for LL, -1 for UR]
    Hy = diags([(v / (2j)) * np.ones(L - 1), (-v / (2j)) * np.ones(L - 1)], [1, -1], format = "lil", dtype = np.complex128)
    Hz = diags([(v / (2j)) * np.ones(L - 1), (-v / (2j)) * np.ones(L - 1)], [1, -1], format = "lil", dtype = np.complex128)

    if (boundary_condition == "PBC"):
      Hx[0, L - 1] = -v / (2j); Hx[L - 1, 0] = v / (2j) # Periodic boundary conditions (PBC)
      Hy[0, L - 1] = -v / (2j); Hy[L - 1, 0] = v / (2j)
      Hz[0, L - 1] = -v / (2j); Hz[L - 1, 0] = v / (2j)

    Hx = Hx.tocsr(); Hy = Hy.tocsr(); Hz = Hz.tocsr()

    # Define the real-space identity matrices
    Ix = eye(L, format = "csr")
    Iy = eye(L, format = "csr")
    Iz = eye(L, format = "csr")

    # assemble the sub-terms to make the total Hamiltonian. 
    # The order of tensor products: |z> X |y> X |x> X |spin>
    Hamiltonian = (kron(Iz, kron(Iy, kron(Hx, sigma_x))) 
                   + kron(Iz, kron(Hy, kron(Ix, sigma_y)))
                   + kron(Hz, kron(Iy, kron(Ix, sigma_z))))

    # Define the on-site chemical potential. 
    #I_0 = eye(Hamiltonian.shape[0], format = "csr")
    I_0 = kron(Iz, kron(Iy, Ix))
    Hamiltonian = Hamiltonian - kron(np.float64(mu)*I_0, sigma_0)

    del Hx, Hy, Hz, Ix, Iy, Iz
    gc.collect()

    return Hamiltonian

  def Dirac_3_over_2(self, L, mu, boundary_condition):
    ####################################################
    # The Hamiltonian represents a Lorentz-symmetric Dirac system 
    # (for higher-spin fermions they may lose Lorentz symmetry depending on how one constructs the Hamiltonian).
    ####################################################
    from scipy.sparse import diags, kron, eye, csr_matrix
    import gc

    # Define the Lorentz-symmetric Gamma matrices (They are analytically derived prior to writing this code)
    Gamma_0 = eye(4, format = "csr") # 4-by-4 identity matrix
    Gamma_x = csr_matrix(np.array([[0, 0, 0, 1], 
                                   [0, 0, 1, 0], 
                                   [0, 1, 0, 0], 
                                   [1, 0, 0, 0]]), dtype = np.complex128)

    Gamma_y = csr_matrix(np.array([[0, 0, 0, +1j], 
                                   [0, 0, -1j, 0], 
                                   [0, +1j, 0, 0], 
                                   [-1j, 0, 0, 0]]), dtype = np.complex128)
    
    Gamma_z = diags([1, -1, 1, -1], format = "csr", dtype = np.complex128) 

    # Define the bonds (i.e., Hamiltonians that resolved in real-space, in direction-wise)
    v = 1.0e0
    Hx = diags([(v / (2j)) * np.ones(L - 1), (-v / (2j)) * np.ones(L - 1)], [1, -1], format = "lil", dtype = np.complex128) # [1 for LL, -1 for UR]
    Hy = diags([(v / (2j)) * np.ones(L - 1), (-v / (2j)) * np.ones(L - 1)], [1, -1], format = "lil", dtype = np.complex128)
    Hz = diags([(v / (2j)) * np.ones(L - 1), (-v / (2j)) * np.ones(L - 1)], [1, -1], format = "lil", dtype = np.complex128)

    if (boundary_condition == "PBC"):
      Hx[0, L - 1] = -v / (2j); Hx[L - 1, 0] = v / (2j) # Periodic boundary conditions (PBC)
      Hy[0, L - 1] = -v / (2j); Hy[L - 1, 0] = v / (2j)
      Hz[0, L - 1] = -v / (2j); Hz[L - 1, 0] = v / (2j)

    Hx = Hx.tocsr(); Hy = Hy.tocsr(); Hz = Hz.tocsr()

    # Define the real-space identity matrices
    Ix = eye(L, format = "csr")
    Iy = eye(L, format = "csr")
    Iz = eye(L, format = "csr")

    # assemble the sub-terms to make the total Hamiltonian.
    # The order of tensor products: |z> X |y> X |x> X |spin>
    Hamiltonian = (kron(Iz, kron(Iy, kron(Hx, Gamma_x)))
                   + kron(Iz, kron(Hy, kron(Ix, Gamma_y)))
                   + kron(Hz, kron(Iy, kron(Ix, Gamma_z)))) 

    # Define the on-site chemical potential.
    I_0 = kron(Iz, kron(Iy, Ix))
    Hamiltonian = Hamiltonian - kron(np.float64(mu)*I_0, Gamma_0)

    del Hx, Hy, Hz, Ix, Iy, Iz
    gc.collect()

    return Hamiltonian

  def Dirac_5_over_2(self, L, mu, boundary_condition):
    ####################################################
    # The Hamiltonian represents a Lorentz-symmetric Dirac system
    # (for higher-spin fermions they may lose Lorentz symmetry depending on how one constructs the Hamiltonian).
    ####################################################
    from scipy.sparse import diags, kron, eye, csr_matrix
    import gc
    
    # Define the Lorentz-symmetric Gamma matrices (They are analytically derived prior to writing this code)
    Gamma_0 = eye(6, format = "csr") # 4-by-4 identity matrix
    Gamma_x = csr_matrix(np.array([[0, 0, 0, 0, 0, 1], 
                                   [0, 0, 0, 0, 1, 0], 
                                   [0, 0, 0, 1, 0, 0], 
                                   [0, 0, 1, 0, 0, 0], 
                                   [0, 1, 0, 0, 0, 0], 
                                   [1, 0, 0, 0, 0, 0]]), dtype = np.complex128)
    
    Gamma_y = csr_matrix(np.array([[0, 0, 0, 0, 0, -1j],
                                   [0, 0, 0, 0, +1j, 0],
                                   [0, 0, 0, -1j, 0, 0],
                                   [0, 0, +1j, 0, 0, 0],
                                   [0, -1j, 0, 0, 0, 0],
                                   [+1j, 0, 0, 0, 0, 0]]), dtype = np.complex128)
    
    Gamma_z = diags([+1, -1, +1, -1, +1, -1], format = "csr", dtype = np.complex128)

    # Define the bonds (i.e., Hamiltonians that resolved in real-space, in direction-wise)
    v = 1.0e0
    Hx = diags([(v / (2j)) * np.ones(L - 1), (-v / (2j)) * np.ones(L - 1)], [1, -1], format = "lil", dtype = np.complex128) # [1 for LL, -1 for UR]
    Hy = diags([(v / (2j)) * np.ones(L - 1), (-v / (2j)) * np.ones(L - 1)], [1, -1], format = "lil", dtype = np.complex128)
    Hz = diags([(v / (2j)) * np.ones(L - 1), (-v / (2j)) * np.ones(L - 1)], [1, -1], format = "lil", dtype = np.complex128)

    if (boundary_condition == "PBC"):
      Hx[0, L - 1] = -v / (2j); Hx[L - 1, 0] = v / (2j) # Periodic boundary conditions (PBC)
      Hy[0, L - 1] = -v / (2j); Hy[L - 1, 0] = v / (2j)
      Hz[0, L - 1] = -v / (2j); Hz[L - 1, 0] = v / (2j)

    Hx = Hx.tocsr(); Hy = Hy.tocsr(); Hz = Hz.tocsr()

    # Define the real-space identity matrices
    Ix = eye(L, format = "csr")
    Iy = eye(L, format = "csr")
    Iz = eye(L, format = "csr")

    # assemble the sub-terms to make the total Hamiltonian.
    # The order of tensor products: |z> X |y> X |x> X |spin>
    Hamiltonian = (kron(Iz, kron(Iy, kron(Hx, Gamma_x)))
                   + kron(Iz, kron(Hy, kron(Ix, Gamma_y)))
                   + kron(Hz, kron(Iy, kron(Ix, Gamma_z))))

    # Define the on-site chemical potential.
    #I_0 = eye(Hamiltonian.shape[0], format = "csr")
    I_0 = kron(Iz, kron(Iy, Ix))
    Hamiltonian = Hamiltonian - kron(np.float64(mu)*I_0, Gamma_0)

    del Hx, Hy, Hz, Ix, Iy, Iz
    gc.collect()

    return Hamiltonian

  def Dirac_2D(self, L, mu, boundary_condition):
    from scipy.sparse import diags, kron, eye, csr_matrix
    import gc
    
    # Define the Pauli matrices sigma_0, sigma_x, sigma_y, and sigma_z
    sigma_0 = eye(2, format="csr") # 2-by-2 identity matrix
    sigma_x = csr_matrix(np.array([[0, 1], [1,  0]], dtype = np.complex128))
    sigma_y = csr_matrix(np.array([[0, -1j], [1j, 0]], dtype = np.complex128))
    sigma_z = csr_matrix(np.array([[1, 0], [0, -1]], dtype = np.complex128))
    
    # Define the bonds (i.e., Hamiltonians that resolved in real-space, in direction-wise)
    t = 1.0e0
    Hx = diags([(t / (2j)) * np.ones(L - 1), (-t / (2j)) * np.ones(L - 1)], [1, -1], format = "lil", dtype = np.complex128) # [1 for LL, -1 for UR]
    Hy = diags([(t / (2j)) * np.ones(L - 1), (-t / (2j)) * np.ones(L - 1)], [1, -1], format = "lil", dtype = np.complex128)

    if (boundary_condition == "PBC"):
      Hx[0, L - 1] = -t / (2j); Hx[L - 1, 0] = t / (2j) # Periodic boundary conditions (PBC)
      Hy[0, L - 1] = -t / (2j); Hy[L - 1, 0] = t / (2j)

    Hx = Hx.tocsr(); Hy = Hy.tocsr()

    # Define the real-space identity matrices
    Ix = eye(L, format = "csr"); Iy = eye(L, format = "csr")

    # assemble the sub-terms to make the total Hamiltonian.
    # The order of tensor products: |y> X |x> X |spin>
    Hamiltonian = kron(Iy, kron(Hx, sigma_x)) + kron(Hy, kron(Ix, sigma_y)) 

    # Define the on-site chemical potential.
    I_0 = kron(Iy, Ix)
    Hamiltonian = Hamiltonian - kron(np.float64(mu)*I_0, sigma_0)

    del Hx, Hy, Ix, Iy
    gc.collect()

    return Hamiltonian

  def Dirac_4D(self, L, mu, boundary_condition):
    from scipy.sparse import diags, kron, eye, csr_matrix
    import gc

    # Define the Pauli matrices sigma_0, sigma_x, sigma_y, and sigma_z
    sigma_0 = eye(2, format="csr") # 2-by-2 identity matrix
    sigma_x = csr_matrix(np.array([[0, 1], [1,  0]], dtype = np.complex128))
    sigma_y = csr_matrix(np.array([[0, -1j], [1j, 0]], dtype = np.complex128))
    sigma_z = csr_matrix(np.array([[1, 0], [0, -1]], dtype = np.complex128))
    
    # Define the 4D Gamma matrices (4-by-4), they are calculated out of Pauli matrices
    Gamma_x = kron(sigma_z, sigma_x); Gamma_y = kron(sigma_z, sigma_y)
    Gamma_z = kron(sigma_z, sigma_z); Gamma_w = kron(sigma_x, sigma_0)
    Gamma_0 = kron(sigma_0, sigma_0) # Equivalently, eye(4, format="csr")

    # Define the bonds (i.e., Hamiltonians that resolved in real-space, in direction-wise)
    t = 1.0e0
    Hx = diags([(t / (2j)) * np.ones(L - 1), (-t / (2j)) * np.ones(L - 1)], [1, -1], format = "lil", dtype = np.complex128) # [1 for LL, -1 for UR]
    Hy = diags([(t / (2j)) * np.ones(L - 1), (-t / (2j)) * np.ones(L - 1)], [1, -1], format = "lil", dtype = np.complex128)
    Hz = diags([(t / (2j)) * np.ones(L - 1), (-t / (2j)) * np.ones(L - 1)], [1, -1], format = "lil", dtype = np.complex128)
    Hw = diags([(t / (2j)) * np.ones(L - 1), (-t / (2j)) * np.ones(L - 1)], [1, -1], format = "lil", dtype = np.complex128)
    
    if (boundary_condition == "PBC"):
      Hx[0, L - 1] = -t / (2j); Hx[L - 1, 0] = t / (2j) # Periodic boundary conditions (PBC)
      Hy[0, L - 1] = -t / (2j); Hy[L - 1, 0] = t / (2j)
      Hz[0, L - 1] = -t / (2j); Hz[L - 1, 0] = t / (2j)
      Hw[0, L - 1] = -t / (2j); Hw[L - 1, 0] = t / (2j)
    
    Hx = Hx.tocsr(); Hy = Hy.tocsr(); Hz = Hz.tocsr(); Hw = Hw.tocsr()

    # Define the real-space identity matrices
    Ix = eye(L, format = "csr"); Iy = eye(L, format = "csr"); Iz = eye(L, format = "csr"); Iw = eye(L, format = "csr")

    # assemble the sub-terms to make the total Hamiltonian.
    # The order of tensor products: |w> X |z> X |y> X |x> X |spin>
    Hamiltonian = (kron(Iw, kron(Iz, kron(Iy, kron(Hx, Gamma_x))))
                   + kron(Iw, kron(Iz, kron(Hy, kron(Ix, Gamma_y))))
                   + kron(Iw, kron(Hz, kron(Iy, kron(Ix, Gamma_z))))
                   + kron(Hw, kron(Iz, kron(Iy, kron(Ix, Gamma_w)))))
    
    # Define the on-site chemical potential.
    I_0 = kron(Iw, kron(Iz, kron(Iy, Ix)))
    Hamiltonian = Hamiltonian - kron(np.float64(mu)*I_0, Gamma_0)

    del Hx, Hy, Hz, Hw, Ix, Iy, Iz, Iw
    gc.collect()

    return Hamiltonian

  def Dirac_5D(self, L, mu, boundary_condition):
    from scipy.sparse import diags, kron, eye, csr_matrix
    import gc
    
    # Define the Pauli matrices sigma_0, sigma_x, sigma_y, and sigma_z
    sigma_0 = eye(2, format="csr") # 2-by-2 identity matrix
    sigma_x = csr_matrix(np.array([[0, 1], [1,  0]], dtype = np.complex128))
    sigma_y = csr_matrix(np.array([[0, -1j], [1j, 0]], dtype = np.complex128))
    sigma_z = csr_matrix(np.array([[1, 0], [0, -1]], dtype = np.complex128))

    # Define the 5D Gamma matrices (also 4-by-4), they are calculated out of Pauli matrices
    Gamma_x = kron(sigma_z, sigma_x); Gamma_y = kron(sigma_z, sigma_y)
    Gamma_z = kron(sigma_z, sigma_z); Gamma_w = kron(sigma_x, sigma_0)
    Gamma_v = kron(sigma_y, sigma_0); Gamma_0 = kron(sigma_0, sigma_0) # Equivalently, eye(4, format="csr")

    # Define the bonds (i.e., Hamiltonians that resolved in real-space, in direction-wise)
    t = 1.0e0
    Hx = diags([(t / (2j)) * np.ones(L - 1), (-t / (2j)) * np.ones(L - 1)], [1, -1], format = "lil", dtype = np.complex128) # [1 for LL, -1 for UR]
    Hy = diags([(t / (2j)) * np.ones(L - 1), (-t / (2j)) * np.ones(L - 1)], [1, -1], format = "lil", dtype = np.complex128)
    Hz = diags([(t / (2j)) * np.ones(L - 1), (-t / (2j)) * np.ones(L - 1)], [1, -1], format = "lil", dtype = np.complex128)
    Hw = diags([(t / (2j)) * np.ones(L - 1), (-t / (2j)) * np.ones(L - 1)], [1, -1], format = "lil", dtype = np.complex128)
    Hv = diags([(t / (2j)) * np.ones(L - 1), (-t / (2j)) * np.ones(L - 1)], [1, -1], format = "lil", dtype = np.complex128)

    if (boundary_condition == "PBC"):
      Hx[0, L - 1] = -t / (2j); Hx[L - 1, 0] = t / (2j) # Periodic boundary conditions (PBC)
      Hy[0, L - 1] = -t / (2j); Hy[L - 1, 0] = t / (2j)
      Hz[0, L - 1] = -t / (2j); Hz[L - 1, 0] = t / (2j)
      Hw[0, L - 1] = -t / (2j); Hw[L - 1, 0] = t / (2j)
      Hv[0, L - 1] = -t / (2j); Hv[L - 1, 0] = t / (2j)

    Hx = Hx.tocsr(); Hy = Hy.tocsr(); Hz = Hz.tocsr(); Hw = Hw.tocsr(); Hv = Hv.tocsr()

    # Define the real-space identity matrices
    Ix = eye(L, format = "csr"); Iy = eye(L, format = "csr"); Iz = eye(L, format = "csr"); Iw = eye(L, format = "csr"); Iv = eye(L, format = "csr")

    # assemble the sub-terms to make the total Hamiltonian.
    # The order of tensor products: |v> X |w> X |z> X |y> X |x> X |spin>
    Hamiltonian = (kron(Iv, kron(Iw, kron(Iz, kron(Iy, kron(Hx, Gamma_x)))))
                   + kron(Iv, kron(Iw, kron(Iz, kron(Hy, kron(Ix, Gamma_y)))))
                   + kron(Iv, kron(Iw, kron(Hz, kron(Iy, kron(Ix, Gamma_z)))))
                   + kron(Iv, kron(Hw, kron(Iz, kron(Iy, kron(Ix, Gamma_w)))))
                   + kron(Hv, kron(Iw, kron(Iz, kron(Iy, kron(Ix, Gamma_v))))))

    # Define the on-site chemical potential.
    I_0 = kron(Iv, kron(Iw, kron(Iz, kron(Iy, Ix))))
    Hamiltonian = Hamiltonian - kron(np.float64(mu)*I_0, Gamma_0)

    del Hx, Hy, Hz, Hw, Hv, Ix, Iy, Iz, Iw, Iv
    gc.collect()

    return Hamiltonian

  def Dirac_6D(self, L, mu, boundary_condition):
    from scipy.sparse import diags, kron, eye, csr_matrix
    import gc

    # Define the Pauli matrices sigma_0, sigma_x, sigma_y, and sigma_z
    sigma_0 = eye(2, format="csr") # 2-by-2 identity matrix
    sigma_x = csr_matrix(np.array([[0, 1], [1,  0]], dtype = np.complex128))
    sigma_y = csr_matrix(np.array([[0, -1j], [1j, 0]], dtype = np.complex128))
    sigma_z = csr_matrix(np.array([[1, 0], [0, -1]], dtype = np.complex128))

    # Define the d < 6 Gamma matrices, they are calculated out of Pauli matrices
    Gamma_x = kron(sigma_z, sigma_x); Gamma_y = kron(sigma_z, sigma_y)
    Gamma_z = kron(sigma_z, sigma_z); Gamma_w = kron(sigma_x, sigma_0)
    Gamma_v = kron(sigma_y, sigma_0); Gamma_0 = kron(sigma_0, sigma_0) # Equivalently, eye(4, format="csr")

    # Define the 8-by-8 Gamma matrices out of Gamma matrices and Pauli matrices. We call them "super_Gamma" matrices
    Super_Gamma_x = kron(sigma_z, Gamma_x); Super_Gamma_y = kron(sigma_z, Gamma_y)
    Super_Gamma_z = kron(sigma_z, Gamma_z); Super_Gamma_w = kron(sigma_z, Gamma_w)
    Super_Gamma_v = kron(sigma_z, Gamma_v); Super_Gamma_u = kron(sigma_x, Gamma_0)
    Super_Gamma_0 = kron(sigma_0, Gamma_0) # Equivalently, eye(8, format="csr")

    # Define the bonds (i.e., Hamiltonians that resolved in real-space, in direction-wise)
    t = 1.0e0
    Hx = diags([(t / (2j)) * np.ones(L - 1), (-t / (2j)) * np.ones(L - 1)], [1, -1], format = "lil", dtype = np.complex128) # [1 for LL, -1 for UR]
    Hy = diags([(t / (2j)) * np.ones(L - 1), (-t / (2j)) * np.ones(L - 1)], [1, -1], format = "lil", dtype = np.complex128)
    Hz = diags([(t / (2j)) * np.ones(L - 1), (-t / (2j)) * np.ones(L - 1)], [1, -1], format = "lil", dtype = np.complex128)
    Hw = diags([(t / (2j)) * np.ones(L - 1), (-t / (2j)) * np.ones(L - 1)], [1, -1], format = "lil", dtype = np.complex128)
    Hv = diags([(t / (2j)) * np.ones(L - 1), (-t / (2j)) * np.ones(L - 1)], [1, -1], format = "lil", dtype = np.complex128)
    Hu = diags([(t / (2j)) * np.ones(L - 1), (-t / (2j)) * np.ones(L - 1)], [1, -1], format = "lil", dtype = np.complex128)

    if (boundary_condition == "PBC"):
      Hx[0, L - 1] = -t / (2j); Hx[L - 1, 0] = t / (2j) # Periodic boundary conditions (PBC)
      Hy[0, L - 1] = -t / (2j); Hy[L - 1, 0] = t / (2j)
      Hz[0, L - 1] = -t / (2j); Hz[L - 1, 0] = t / (2j)
      Hw[0, L - 1] = -t / (2j); Hw[L - 1, 0] = t / (2j)
      Hv[0, L - 1] = -t / (2j); Hv[L - 1, 0] = t / (2j)
      Hu[0, L - 1] = -t / (2j); Hu[L - 1, 0] = t / (2j)

    Hx = Hx.tocsr(); Hy = Hy.tocsr(); Hz = Hz.tocsr(); Hw = Hw.tocsr(); Hv = Hv.tocsr(); Hu = Hu.tocsr()

    # Define the real-space identity matrices
    Ix = eye(L, format = "csr"); Iy = eye(L, format = "csr"); Iz = eye(L, format = "csr"); 
    Iw = eye(L, format = "csr"); Iv = eye(L, format = "csr"); Iu = eye(L, format = "csr")

    # assemble the sub-terms to make the total Hamiltonian.
    # The order of tensor products: |u> X |v> X |w> X |z> X |y> X |x> X |spin>
    Hamiltonian = (kron(Iu, kron(Iv, kron(Iw, kron(Iz, kron(Iy, kron(Hx, Super_Gamma_x)))))) 
                   + kron(Iu, kron(Iv, kron(Iw, kron(Iz, kron(Hy, kron(Ix, Super_Gamma_y))))))
                   + kron(Iu, kron(Iv, kron(Iw, kron(Hz, kron(Iy, kron(Ix, Super_Gamma_z))))))
                   + kron(Iu, kron(Iv, kron(Hw, kron(Iz, kron(Iy, kron(Ix, Super_Gamma_w))))))
                   + kron(Iu, kron(Hv, kron(Iw, kron(Iz, kron(Iy, kron(Ix, Super_Gamma_v)))))) 
                   + kron(Hu, kron(Iv, kron(Iw, kron(Iz, kron(Iy, kron(Ix, Super_Gamma_u))))))
                   )

    # Define the on-site chemical potential.
    I_0 = kron(Iu, kron(Iv, kron(Iw, kron(Iz, kron(Iy, Ix)))))
    Hamiltonian = Hamiltonian - kron(np.float64(mu)*I_0, Super_Gamma_0)

    del Hx, Hy, Hz, Hw, Hv, Hu, Ix, Iy, Iz, Iw, Iv, Iu
    gc.collect()

    return Hamiltonian

  def Dirac_3D_Hidden_order(self, L, mu, alpha_Hidden, boundary_condition):
    #####################################################
    # Theory first discovered and proposed in Roy, Bitan, arXiv:2606.27368v1 (2026)
    # The spin orbital is extended to four although the fermions may still be s-1/2.
    # The hidden order is induced by a symmetry-breaking mass term at the end of the tight-binding like term. 
    #
    # The mass matrix, M, is taken to do matrix-product with the Gamma matrices, before tensor product to form the complete Hilbert space. 
    #####################################################
    from scipy.sparse import diags, kron, eye, csr_matrix
    import gc

    # Define the Pauli matrices sigma_0, sigma_x, sigma_y, and sigma_z
    sigma_0 = eye(2, format="csr") # 2-by-2 identity matrix
    sigma_x = csr_matrix(np.array([[0, 1], [1,  0]], dtype = np.complex128))
    sigma_y = csr_matrix(np.array([[0, -1j], [1j, 0]], dtype = np.complex128))
    sigma_z = csr_matrix(np.array([[1, 0], [0, -1]], dtype = np.complex128))

    # Define the Gamma matrices and the mass matrix, they are calculated out of Pauli matrices
    # Strictly speaking, the first matrices in each tensor product are pseudo-spin tau-matrices; 
    # nonetheless, the tau-matrices have exactly the same numerical values as the Pauli matrices. 
    Gamma_x = kron(sigma_z, sigma_x); Gamma_y = kron(sigma_z, sigma_y)
    Gamma_z = kron(sigma_z, sigma_z); Gamma_0 = kron(sigma_0, sigma_0) # Equivalently, eye(4, format="csr")
    Mass_matrix = kron(sigma_x, sigma_0) # the mass matrix M

    # Define the matrix product of M and Gamma
    Mass_Gamma_x = Mass_matrix @ Gamma_x; Mass_Gamma_y = Mass_matrix @ Gamma_y; Mass_Gamma_z = Mass_matrix @ Gamma_z
    
    t = 1.0e0
    #  Define the normal, Dirac bonds (i.e., Hamiltonians that resolved in real-space, in direction-wise)
    Hx = diags([(t / (2j)) * np.ones(L - 1), (-t / (2j)) * np.ones(L - 1)], [1, -1], format = "lil", dtype = np.complex128) # [1 for LL, -1 for UR]
    Hy = diags([(t / (2j)) * np.ones(L - 1), (-t / (2j)) * np.ones(L - 1)], [1, -1], format = "lil", dtype = np.complex128)
    Hz = diags([(t / (2j)) * np.ones(L - 1), (-t / (2j)) * np.ones(L - 1)], [1, -1], format = "lil", dtype = np.complex128)

    # Define the bonds for hidden order
    Hx_M = diags([(1j * t * alpha_Hidden / (2j)) * np.ones(L - 1), (-1j * t * alpha_Hidden / (2j)) * np.ones(L - 1)], 
                  [1, -1], format = "lil", dtype = np.complex128) # [1 for LL, -1 for UR]
    Hy_M = diags([(1j * t * alpha_Hidden / (2j)) * np.ones(L - 1), (-1j * t * alpha_Hidden / (2j)) * np.ones(L - 1)], 
                  [1, -1], format = "lil", dtype = np.complex128)
    Hz_M = diags([(1j * t * alpha_Hidden / (2j)) * np.ones(L - 1), (-1j * t * alpha_Hidden / (2j)) * np.ones(L - 1)], 
                  [1, -1], format = "lil", dtype = np.complex128)

    if (boundary_condition == "PBC"): # Periodic boundary conditions (PBC)
      Hx[0, L - 1] = -t / (2j); Hx[L - 1, 0] = t / (2j)
      Hy[0, L - 1] = -t / (2j); Hy[L - 1, 0] = t / (2j)
      Hz[0, L - 1] = -t / (2j); Hz[L - 1, 0] = t / (2j)

      Hx_M[0, L - 1] = -1j * t * alpha_Hidden / (2j); Hx_M[L - 1, 0] = 1j * t * alpha_Hidden / (2j)
      Hy_M[0, L - 1] = -1j * t * alpha_Hidden / (2j); Hy_M[L - 1, 0] = 1j * t * alpha_Hidden / (2j)
      Hz_M[0, L - 1] = -1j * t * alpha_Hidden / (2j); Hz_M[L - 1, 0] = 1j * t * alpha_Hidden / (2j)

    Hx = Hx.tocsr(); Hy = Hy.tocsr(); Hz = Hz.tocsr()
    Hx_M = Hx_M.tocsr(); Hy_M = Hy_M.tocsr(); Hz_M = Hz_M.tocsr()

    # Define the real-space identity matrices
    Ix = eye(L, format = "csr"); Iy = eye(L, format = "csr"); Iz = eye(L, format = "csr")

    # assemble the sub-terms to make the total Hamiltonian.
    # The order of tensor products: |z> X |y> X |x> X |spin>
    Hamiltonian_Dir = (kron(Iz, kron(Iy, kron(Hx, Gamma_x)))
                      + kron(Iz, kron(Hy, kron(Ix, Gamma_y)))
                      + kron(Hz, kron(Iy, kron(Ix, Gamma_z)))) # The Dirac, tight-binding like part of Hamiltonian. 

    Hamiltonian_Mass = (kron(Iz, kron(Iy, kron(Hx_M, Mass_Gamma_x)))
                        + kron(Iz, kron(Hy_M, kron(Ix, Mass_Gamma_y)))
                        + kron(Hz_M, kron(Iy, kron(Ix, Mass_Gamma_z))))

    Hamiltonian = Hamiltonian_Dir + Hamiltonian_Mass

    # Define the on-site chemical potential.
    I_0 = kron(Iz, kron(Iy, Ix))
    Hamiltonian = Hamiltonian - kron(np.float64(mu)*I_0, Gamma_0)

    del Hx, Hy, Hz, Hx_M, Hy_M, Hz_M, Ix, Iy, Iz, Hamiltonian_Dir, Hamiltonian_Mass
    gc.collect()

    return Hamiltonian

  def Dirac_3D_Hidden_order_ver2(self, L, mu, alpha_Hidden, boundary_condition):
    #####################################################
    # This Hamiltonian will be the same as the one produced in "Dirac_3D_Hidden_order"
    # It's just the symmetry-breaking term is constructed through a slightly different way.
    #
    # The mass matrix, M, is first taken tensor product with "I_0", the identity matrix in real space, and then taken 
    # matrix product with the kinetic part "Hamiltonian_Dir". 
    #####################################################
    from scipy.sparse import diags, kron, eye, csr_matrix
    import gc

    # Define the Pauli matrices sigma_0, sigma_x, sigma_y, and sigma_z
    sigma_0 = eye(2, format="csr") # 2-by-2 identity matrix
    sigma_x = csr_matrix(np.array([[0, 1], [1,  0]], dtype = np.complex128))
    sigma_y = csr_matrix(np.array([[0, -1j], [1j, 0]], dtype = np.complex128))
    sigma_z = csr_matrix(np.array([[1, 0], [0, -1]], dtype = np.complex128))

    # Define the Gamma matrices and the mass matrix, they are calculated out of Pauli matrices
    # Strictly speaking, the first matrices in each tensor product are pseudo-spin tau-matrices;
    # nonetheless, the tau-matrices have exactly the same numerical values as the Pauli matrices.
    Gamma_x = kron(sigma_z, sigma_x); Gamma_y = kron(sigma_z, sigma_y)
    Gamma_z = kron(sigma_z, sigma_z); Gamma_0 = kron(sigma_0, sigma_0) # Equivalently, eye(4, format="csr")
    Mass_matrix = kron(sigma_x, sigma_0) # the mass matrix M

    t = 1.0e0
    #  Define the normal, Dirac bonds (i.e., Hamiltonians that resolved in real-space, in direction-wise)
    Hx = diags([(t / (2j)) * np.ones(L - 1), (-t / (2j)) * np.ones(L - 1)], [1, -1], format = "lil", dtype = np.complex128) # [1 for LL, -1 for UR]
    Hy = diags([(t / (2j)) * np.ones(L - 1), (-t / (2j)) * np.ones(L - 1)], [1, -1], format = "lil", dtype = np.complex128)
    Hz = diags([(t / (2j)) * np.ones(L - 1), (-t / (2j)) * np.ones(L - 1)], [1, -1], format = "lil", dtype = np.complex128)

    if (boundary_condition == "PBC"): # Periodic boundary conditions (PBC)
      Hx[0, L - 1] = -t / (2j); Hx[L - 1, 0] = t / (2j)
      Hy[0, L - 1] = -t / (2j); Hy[L - 1, 0] = t / (2j)
      Hz[0, L - 1] = -t / (2j); Hz[L - 1, 0] = t / (2j)

    Hx = Hx.tocsr(); Hy = Hy.tocsr(); Hz = Hz.tocsr()

    # Define the real-space identity matrices
    Ix = eye(L, format = "csr"); Iy = eye(L, format = "csr"); Iz = eye(L, format = "csr")
    I_0 = kron(Iz, kron(Iy, Ix))

    # assemble the sub-terms to make the total Hamiltonian.
    # The order of tensor products: |z> X |y> X |x> X |spin>
    Hamiltonian_Dir = (kron(Iz, kron(Iy, kron(Hx, Gamma_x)))
                      + kron(Iz, kron(Hy, kron(Ix, Gamma_y)))
                      + kron(Hz, kron(Iy, kron(Ix, Gamma_z)))) # The Dirac, tight-binding like part of Hamiltonian.
    
    Hamiltonian = Hamiltonian_Dir + (1j) * t * alpha_Hidden * (kron(I_0, Mass_matrix)) @ Hamiltonian_Dir

    del Hx, Hy, Hz, I_0, Ix, Iy, Iz, Hamiltonian_Dir
    gc.collect()

    return Hamiltonian

  def make_random(self, W, distribution, size, dimension, spin, Hidden_order):
    from scipy.sparse import diags, kron, eye, csr_matrix
    
    # Get a list of random number
    if distribution == "box":
      dis_pot_array = np.random.uniform(-W / 2.0E0, W / 2.0E0, size = size)
    elif distribution == "binary":
      dis_pot_array = np.random.choice([-W, W], size = size)
    elif distribution == "chiral":
      if (Hidden_order == False):
        raise ValueError("The chiral disorder distribution only works when Hidden_order = True !")
      dis_pot_array = np.random.uniform(-W / 2.0E0, W / 2.0E0, size = size)
    else:
      raise ValueError("Unknown distribution: %s" % distribution)

    # Promote the list of random number to disorder potential. The disorder potential is a diagonal matrix in the same dimension with the clean Hamiltonian. 
    if (dimension == 2):
      sigma_0 = eye(2, format = "csr")
      dis_pot = diags(dis_pot_array, format = "csr", dtype = np.float64)
      dis_pot_matrix = kron(dis_pot, sigma_0) # it's a diagonal, (2N)-by-(2N) matrix

    elif (dimension == 3):
      if (spin == 1): # for spin-1/2 systems
        if (Hidden_order == False): 
          # Create identity matrix in spin space for disordered potential
          sigma_0 = eye(2, format = "csr")
          dis_pot = diags(dis_pot_array, format = "csr", dtype = np.float64)
          dis_pot_matrix = kron(dis_pot, sigma_0) # it's a diagonal, (2N)-by-(2N) matrix
        else: # when Hidden_order == True
          if (distribution == "box"):
            Gamma_0 = eye(4, format = "csr")
            dis_pot = diags(dis_pot_array, format = "csr", dtype = np.float64)
            dis_pot_matrix = kron(dis_pot, Gamma_0) # it's a diagonal, (4N)-by-(4N) matrix
          elif (distribution == "chiral"):
            sigma_0 = eye(2, format = "csr"); sigma_z = csr_matrix(np.array([[1, 0], [0, -1]], dtype = np.complex128))
            Gamma_5 = kron(sigma_z, sigma_0)
            dis_pot = diags(dis_pot_array, format = "csr", dtype = np.float64)
            dis_pot_matrix = kron(dis_pot, Gamma_5) # it's a diagonal, (4N)-by-(4N) matrix

      elif (spin == 3): # for spin-3/2 systems
        # Create identity matrix in spin space for disordered potential
        Gamma_0 = eye(4, format = "csr")
        dis_pot = diags(dis_pot_array, format = "csr", dtype = np.float64)
        dis_pot_matrix = kron(dis_pot, Gamma_0) # it's a diagonal, (4N)-by-(4N) matrix

      elif (spin == 5): # for spin-5/2 systems
        # Create identity matrix in spin space for disordered potential
        Gamma_0 = eye(6, format = "csr")
        dis_pot = diags(dis_pot_array, format = "csr", dtype = np.float64)
        dis_pot_matrix = kron(dis_pot, Gamma_0) # it's a diagonal, (6N)-by-(6N) matrix

    elif (dimension == 4): 
      Gamma_0 = eye(4, format = "csr")
      dis_pot = diags(dis_pot_array, format = "csr", dtype = np.float64)
      dis_pot_matrix = kron(dis_pot, Gamma_0) # it's a diagonal, (4N)-by-(4N) matrix

    elif (dimension == 5):
      Gamma_0 = eye(4, format = "csr")
      dis_pot = diags(dis_pot_array, format = "csr", dtype = np.float64)
      dis_pot_matrix = kron(dis_pot, Gamma_0) # it's a diagonal, (4N)-by-(4N) matrix

    elif (dimension == 6):
      Super_Gamma_0 = eye(8, format = "csr")
      dis_pot = diags(dis_pot_array, format = "csr", dtype = np.float64)
      dis_pot_matrix = kron(dis_pot, Super_Gamma_0) # it's a diagonal, (8N)-by-(8N) matrix

    return dis_pot_matrix

