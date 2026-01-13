"""
X_CIC_cronos.py 
==> @ Created by J.M. Abril-Hernández, 2025.

This script reads the following files: 'configuration.json' (containing the core data),
random samples from the library, 'Absolute_min.txt', and 'Map3D.txt'. It then computes
the 1‑sigma interval and writes all solutions within it to 'Cloud.txt'. This file is
subsequently processed by the 'profile' function, which generates full profiles with
model magnitudes, including the chronology. These solutions are stored in 'Plot.txt',
ready for plotting with Gnuplot.

Options for the confidence region (declared in the .json file):
 --interesting-parameters -> false => coefficient 3.53 (default), true => coefficient 2.30
 
 Note that the sampling-step is dynamically determined to sample fewer than 6000 solutions from Cloud.txt

"""
import numpy as np
import sys
import json

# Reading the configuration file JSON
with open("configuration.json", "r") as archivo:
    datos = json.load(archivo)

# Assigning values to the variables
Core_data = datos["Core_data"]
OP_LN_w = datos["OP_LN_w"]  
kr = datos["kr"]
Tmr = datos["Tmr"]
sgt = datos["sgt"]
peso = datos["peso"]
interesting_parameters = datos["interesting_parameters"]



# --- Read Absolute_min.txt ---
with open('Absolute_min.txt', 'r') as file:
    lineas = [line.strip().split() for line in file if line.strip()]
# First line 
min_A0, min_w, min_sw, minChi, DoF = (
    float(lineas[0][0]), float(lineas[0][1]), float(lineas[0][2]),
    float(lineas[0][3]), float(lineas[0][4])
)

# --- Compute Δχ for 1-sigma based on the selected coefficient ---
coef = 2.30 if interesting_parameters else 3.53
D_Chi_1_sigma = (minChi**2 + coef/(DoF))**0.5 - minChi
# --- 1-sigma upper threshold ---
Chi_1_sigma = minChi + D_Chi_1_sigma
print("1-sigma threshold:")
print(f"minChi = {minChi:.6f}")
print(f"minChi_1sg = {Chi_1_sigma:.6f}")
print(f"Coefficient used: {coef} ({'interesting' if interesting_parameters else 'conservative'})")
# --- Process Mapa4D.txt and write single output file ---
count_in = 0
count_out = 0

with open('Cloud.txt', 'w', encoding='utf-8') as fout, open('Map3D.txt', 'r') as fin:
    for line in fin:
        if not line.strip():
            continue
        parts = line.strip().split()
        if len(parts) != 4:
            # Skip malformed lines instead of stopping
            continue
        A_solver = float(parts[0])
        w_solver = float(parts[1])
        sw_solver = float(parts[2])
        vChi = float(parts[3])
        if interesting_parameters:
            # Requires that A_solver == min_A0 and sA_solver == min_sA
            if vChi <= Chi_1_sigma and A_solver == min_A0:
                fout.write(line)
                count_in += 1
            else:
                count_out += 1
        else:
            # Default case (3 parameters): vChi <= Chi_1_sigma
            if vChi <= Chi_1_sigma:
                fout.write(line)
                count_in += 1
            else:
                count_out += 1
# --- Summary ---
print(f"Solvers within 1-sigma (Cloud): {count_in}")
print(f"Outside 1-sigma: {count_out}")
print(f"Total solvers processed: {count_in + count_out}")

# Dinamical fixation of sampling step
sampling_step = int(count_in / 6000 + 1)
print(f"[info] sampling_step {sampling_step}")

# Below the solvers stored in 'Cloud.txt' are processed by the 'profile' fuction.

m_i = []   # mass depth at the bottom of the slice i
Aexp_i = []  # 210Pb_exc (Bq/kg) in slice i
sgAexp_i = []  # Error in 210Pb_exc (Bq/kg) in slice i
with open(Core_data, 'r') as file:
    for line in file:
        if line.strip():
            parts = line.strip().split()
            m_i.append(float(parts[0]))
            Aexp_i.append(float(parts[1]))
            sgAexp_i.append(float(parts[2]))
# Reading the canonical representative sample of randomly sorted values following a normal typified distribution
N = len(m_i)
file2 = f"./aleat_S1/aleat_S1_{N}.txt"
print(file2)
z_1 = []
with open(file2, 'r') as file_aleat:
    for line in file_aleat:
        if line.strip():
            z_1.append(float(line.strip()))
z_1 = np.array(z_1)

ldPb= 0.03118 # 1/yr. decay constant for 210Pb
# Computing the mass depth scale referred to the mid-point of each sediment slice
mi_m = []
for k in range(N):
    if k == 0:
        dm2=m_i[0]
        mi_m.append(dm2/2)
    else:
        dm2 = m_i[k]-m_i[k-1]
        mi_m.append(m_i[k-1]+dm2/2)
# print(mi_m)
def profile(A0_solver, w_solver, sw_solver):
    """
    Given (A0_solver, w_solver, sw_solver, ), generates N pairs (A0i, wi) and finds their best ordering
    that minimizes the distance to the experimental profile, including radioactive decay
    and an optional time-mark constraint.
    Returns 
    Sol_A, Sol_w, Sol_T, Sol_Ath, Sol_flux
    """
    Sol_T = [] # Solution chronology (cumulative age by slice)
    Sol_A = [] # Ordered initial activities
    Sol_w = [] # Ordered SARs
    Sol_Ath = [] # Theoretical profile (solution including decay)
    Sol_flux = [] # Stores solution for the 210-Pb_exc flux captured by each sediment slice.
    slices = N
    chif = 0.0
    tant = 0.0
    # --- w ---
    if OP_LN_w:
        # Parameters already in log-space; sA_solver is relative in log-space
        wi = w_solver * (1.0 + sw_solver * z_1)
        wi = np.exp(wi)
    else:
        # Standard (physical) space with relative deviation
        wi = w_solver * (1.0 + sw_solver * z_1)
    # --- A0 ---
    A0i = np.full(N, A0_solver)

    # wi = np.clip(wi, 0.15 * w_central, None) # guard against non-physical small values. [Activate this line when needed]
    # A0i = np.clip(A0i, 0.1 * A0_central, None) # guard against non-physical small values. [Activate this line when needed]
    for k in range(N):
        chiant = 1.0E20
        jmin = 0
        if k == 0:
            dm2 = m_i[0]
        else:
            dm2 = m_i[k] - m_i[k - 1]
        for s in range(slices):
            Dt = dm2 / wi[s]
            Act = A0i[s] * np.exp(-ldPb * tant) * (1 - np.exp(-ldPb * Dt)) / (ldPb * Dt)
            chi = (Act - Aexp_i[k])**2
            if chi < chiant:
                jmin = s
                chiant = chi
        Sol_A.append(A0i[jmin])
        Sol_w.append(wi[jmin])
        Sol_T.append(tant + dm2 / wi[jmin])
        tant = tant + dm2 / wi[jmin]
        A0i = np.delete(A0i, jmin)
        wi = np.delete(wi, jmin)
        slices -= 1
        if k == 0:
            Dt = Sol_T[0]
            tanterior = 0
        else:
            Dt = Sol_T[k] - Sol_T[k - 1]
            tanterior = Sol_T[k - 1]
        Act = Sol_A[k] * np.exp(-ldPb * tanterior) * (1 - np.exp(-ldPb * Dt)) / (ldPb * Dt)
        chif += ((Act - Aexp_i[k]) / sgAexp_i[k])**2
        Sol_Ath.append(Act)
        Sol_flux.append(Sol_A[k]*Sol_w[k]*10)
    # End k-loop
    # Objective function definition (optional time-mark)
    chif = chif + peso * ((Tmr - Sol_T[kr]) / sgt)**2 # see more in publications
    chi_df = (chif / DoF)**0.5
    return Sol_A, Sol_w, Sol_T, Sol_Ath, Sol_flux
# end definition

# --- Plot.txt ---
with open('Plot.txt', 'w', encoding= 'utf-8') as file_out1:
    with open('Cloud.txt', 'r') as file1:
        itcount_plot = 0
        for line in file1:
            if line.strip():
                itcount_plot += 1
                # Application of sampling_step
                if sampling_step > 1 and (itcount_plot % sampling_step) != 0:
                    continue
                parts = line.strip().split()
                A0_solver = float(parts[0])
                w_solver = float(parts[1])
                sw_solver = float(parts[2])
                #sw_solver = float(parts[3])
                Sol_A, Sol_w, Sol_T, Sol_Ath, Sol_flux = profile(A0_solver, w_solver, sw_solver)
                for x1, x2, x3, x4, x5, x6, x7 in zip(m_i, mi_m, Sol_A, Sol_w, Sol_T, Sol_Ath, Sol_flux):
                    file_out1.write(f"{x1:.3f}\t{x2:.3f}\t{x3:.2f}\t{x4:.4f}\t{x5:.2f}\t{x6:.2f}\t{x7:.2f}\n")
                file_out1.write(f"\n")

with open('Solution.txt', 'w', encoding= 'utf-8') as file_out4:
    with open('Absolute_min.txt', 'r') as file4:
        for line in file4:
            if line.strip():
                parts = line.strip().split()
                A0_solver = float(parts[0])
                w_solver = float(parts[1])
                sw_solver = float(parts[2])
                #sw_solver = float(parts[3])
                # Additional information in this file is not needed here
                Sol_A, Sol_w, Sol_T, Sol_Ath, Sol_flux = profile(A0_solver, w_solver, sw_solver)
                for x1, x2, x3, x4, x5, x6, x7 in zip(m_i, mi_m, Sol_A, Sol_w, Sol_T, Sol_Ath, Sol_flux):
                    file_out4.write(f"{x1:.3f}\t{x2:.3f}\t{x3:.2f}\t{x4:.4f}\t{x5:.2f}\t{x6:.2f}\t{x7:.2f}\n")
                break
sys.exit()
