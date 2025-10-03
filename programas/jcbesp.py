import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import h, c, k
import os

def planck(wl, T):
    a = 2 * h * c**2 / wl**5
    b = h * c / (wl * k * T)
    return a / (np.exp(b) - 1 + 1e-20)

current_dir = os.path.dirname(os.path.abspath(__file__))
data_file = os.path.join(current_dir, 'jcb.dat')

r_in = r_out = d_puntos = None
estrellas = []
disco = []

try:
    with open(data_file, 'r') as f:
        for line in f:
            if line.startswith('#DSIZE'):
                partes = line.strip().split()
                r_in = float(partes[1])
                r_out = float(partes[2])
                d_puntos = float(partes[3])
            elif line.startswith('#STAR'):
                partes = line.strip().split()
                estrellas.append({
                    'R': float(partes[1]),
                    'T': float(partes[2])
                })
            elif not line.startswith('#'):
                partes = line.strip().split()
                if len(partes) == 4:
                    x, y, zona, T = map(float, partes)
                    r = np.sqrt(x**2 + y**2)
                    if r_in <= r <= r_out:
                        disco.append(T)

except FileNotFoundError:
    print(f"Error: Archivo {data_file} no encontrado")
    exit()

if not estrellas or not disco:
    raise ValueError("Datos incompletos o corruptos en el archivo")

wl = np.logspace(-9, -3, 1000)
# Distancia al sistema en metros: 2.908e7 UA * 1.5e11 m/UA = 4.362e18 m (aprox. 14.1 parsecs)
d = 2.908e7 * 1.5e11  # metros

# Convertir d_puntos de km a m
d_puntos_m = d_puntos * 1000
area_elemento = d_puntos_m**2  # ahora en m²

# =======================================================================
# 1. Espectro de estrellas y disco binario
# =======================================================================
espectros_estrellas = []
for i, estrella in enumerate(estrellas):
    R_m = estrella['R'] * 1000  # radio de km a m
    B = planck(wl, estrella['T'])
    flujo = np.pi * B * (R_m**2 / d**2)  # flujo en W/m²/m
    espectros_estrellas.append(flujo)

espectro_disco = np.zeros_like(wl)
for T_elemento in disco:
    B = planck(wl, T_elemento)
    espectro_disco += B * (area_elemento / d**2)

# Espectro combinado de las estrellas (sin disco)
espectro_estrellas_combinado = sum(espectros_estrellas)

# =======================================================================
# 2. Graficación con estilo unificado
# =======================================================================
plt.figure(figsize=(14, 7))

# Sistema binario completo (estrellas + disco) - Tono más oscuro
plt.plot(wl*1e9, espectro_estrellas_combinado + espectro_disco, 
         color='darkblue', lw=3, 
         label='Sistema binario completo (estrellas + disco)')

# Espectro combinado de las estrellas (sin disco) - Mismo color, línea punteada
plt.plot(wl*1e9, espectro_estrellas_combinado, 
         color='blue', linestyle='--', lw=2,
         label='Estrellas combinadas (sin disco)')

# Disco solo - Mismo color, línea punteada con puntos
plt.plot(wl*1e9, espectro_disco, 
         color='blue', linestyle='-.', lw=2, 
         label='Disco (sin estrellas)')

plt.xscale('log')
plt.yscale('log')
plt.xlabel('Longitud de onda (nm)', fontsize=12)
plt.ylabel('Flujo espectral (W/m²/nm)', fontsize=12)
plt.title(f'Comparación de espectros - Sistema binario\nDisco: {r_in/1.496e8:.1f}-{r_out/1.496e8:.1f} UA', 
          fontsize=14, pad=20)

# Ajustar los límites de los ejes según lo solicitado
plt.xlim(10**1, 10**7)  # Longitud de onda entre 10^1 y 10^7 nm
plt.ylim(10**(-14), 10**(-3))  # Flujo espectral entre 10^-14 y 10^-3 W/m²/nm

plt.legend(fontsize=10, loc='upper right')
plt.grid(True, which='both', alpha=0.3)
plt.tight_layout()

plt.savefig('jcb_esp.png', dpi=300)
#plt.show()