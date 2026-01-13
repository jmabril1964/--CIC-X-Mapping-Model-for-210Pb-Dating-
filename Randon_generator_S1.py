"""
This coge generates the library of random numbers required for
X-CIC, a X-mapping model for 210Pb-dating assuming a constant initial activity concentration
Created by J.M. Abril-Hernández, 2025.
"""

import numpy as np
from scipy.stats import norm
import random

def generar_normal_canonica(N):
    """Generate the canonical representative sample for a standard normal distribution."""
    probabilidades = [(i + 0.5) / N for i in range(N)]
    valores = [norm.ppf(p) for p in probabilidades]
    return np.array(valores)

# Generate the library of random samples
for N in range(5, 100):
    filename = f"aleat_S1_{N}.txt"   
    z_0 = generar_normal_canonica(N)
    z_1 = z_0.copy()
    random.shuffle(z_1)

    with open(filename, 'w', encoding='utf-8') as archivo:
        for x in z_1:
            archivo.write(f"{x:.4f}\n")

print("Done!")

