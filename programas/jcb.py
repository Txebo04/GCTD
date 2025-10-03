import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os

def plot_sistema_estelar():
    # Configuración de figura 3D
    fig = plt.figure(figsize=(14, 10))
    ax = fig.add_subplot(111, projection='3d')
    
    # Paleta de colores Inferno
    cmap = plt.cm.cubehelix
    
    # Constantes de conversión
    KM_TO_AU = 1.496e8  # 1 UA = 1.496e8 km
    RSUN_TO_KM = 695700  # 1 Radio solar en km
    
    # Carga de datos
    script_dir = os.path.dirname(os.path.abspath(__file__))
    dat_path = os.path.join(script_dir, "jcb.dat")
    
    if not os.path.exists(dat_path):
        raise FileNotFoundError(f"Archivo no encontrado: {dat_path}")
    
    # Lectura y procesamiento de datos
    estrellas = []
    radios_km = []  # Almacenar radios en km
    cm = None
    puntos = []
    temperaturas = []
    
    with open(dat_path, 'r') as f:
        for line in f:
            if line.startswith('#POSSTAR'):
                # Convertir posiciones de km a UA
                pos_km = list(map(float, line.split()[1:4]))
                pos_au = [coord / KM_TO_AU for coord in pos_km]
                estrellas.append(pos_au)
            elif line.startswith('#STAR'):
                # Leer radio y temperatura de la estrella
                datos_star = list(map(float, line.split()[1:3]))
                radios_km.append(datos_star[0])  # Radio en km
            elif line.startswith('#CM'):
                # Convertir posición del CM de km a UA
                cm_km = list(map(float, line.split()[1:4]))
                cm = [coord / KM_TO_AU for coord in cm_km]
            elif not line.startswith('#'):
                datos = list(map(float, line.split()))
                # Convertir coordenadas de km a UA
                x_au = datos[0] / KM_TO_AU
                y_au = datos[1] / KM_TO_AU
                puntos.append([x_au, y_au, 0.0])  # Z siempre 0 en plano XY
                temperaturas.append(datos[3])  # Cuarta columna es temperatura
    
    puntos = np.array(puntos)
    temperaturas = np.array(temperaturas)
    estrellas = np.array(estrellas)
    
    # Construcción de la malla estructurada
    x_vals = np.unique(puntos[:,0])
    y_vals = np.unique(puntos[:,1])
    X, Y = np.meshgrid(x_vals, y_vals)
    
    # Crear matriz de temperaturas con forma de malla
    Temp_grid = np.full(X.shape, np.nan)
    for idx, (x, y, _) in enumerate(puntos):
        i = np.where(x_vals == x)[0][0]
        j = np.where(y_vals == y)[0][0]
        Temp_grid[j,i] = temperaturas[idx]
    
    # Normalización para el colormap
    norm = plt.Normalize(vmin=np.nanmin(Temp_grid), vmax=np.nanmax(Temp_grid))
    
    # Crear capa de transparencia
    alpha_layer = np.where(np.isnan(Temp_grid), 0.0, 1.0)
    rgba = cmap(norm(Temp_grid))
    rgba[..., 3] = alpha_layer
    
    # Plot 3D del disco
    ax.plot_surface(X, Y, np.zeros_like(X), 
                   facecolors=rgba, 
                   rstride=1, 
                   cstride=1, 
                   antialiased=False,
                   shade=False)
    
    # Dibujar estrellas a escala real
    estrella_colors = ['gold', 'orange']
    estrella_labels = ['Estrella 1', 'Estrella 2']
    
    for i, estrella in enumerate(estrellas):
        if len(estrella) >= 3:
            x, y, z = estrella[:3]
            # Convertir radio de km a UA y escalar para visualización
            if i < len(radios_km):
                radio_km = radios_km[i]
                radio_au = radio_km / KM_TO_AU
                # Escalar el tamaño para visualización (factor de escala ajustable)
                marker_size = 3000 * radio_au  # Ajusta este factor según sea necesario
            else:
                marker_size = 100
                
            ax.scatter(x, y, z, s=marker_size, c=estrella_colors[i],
                      marker='o', edgecolor='black', alpha=1.0, label=estrella_labels[i])
    
    # Centro de masa
    ax.scatter(cm[0], cm[1], cm[2], c='red', 
              s=20, marker='X', label='CM')
    
    # Configuración estética con unidades en UA
    ax.set_title("Distribución de Temperatura en el Disco", fontsize=14)
    ax.set_xlabel('X (UA)', labelpad=15)
    ax.set_ylabel('Y (UA)', labelpad=15)
    ax.set_zlabel('Z (UA)', labelpad=10)
    
    # Límites y perspectiva
    max_range = np.max(np.abs(puntos[:,:2]))
    ax.set_xlim([-max_range, max_range])
    ax.set_ylim([-max_range, max_range])
    ax.set_zlim([-max_range*0.05, max_range*0.05])
    ax.view_init(elev=40, azim=-45)
    
    # Barra de color para temperatura
    cbar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax,
                       shrink=0.7, pad=0.1)
    cbar.set_label('Temperatura (K)', rotation=270, labelpad=20)
    
    # Leyenda
    ax.legend(loc='upper right')
    
    plt.savefig('pos1.png', dpi=300, bbox_inches='tight')
    # plt.show()

if __name__ == "__main__":
    plot_sistema_estelar()