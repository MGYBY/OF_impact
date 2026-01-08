import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np

# Choose a cmap
cmap_name = 'brg'
cmap = plt.get_cmap(cmap_name)
num_colors = 50 # Number of control points

# Sample the colormap
colors = cmap(np.linspace(0, 1, num_colors))

# Create .map file
with open(f'{cmap_name}.map', 'w') as f:
    f.write("#!MC 1410\n")
    f.write("$! COLORMAP CONTOURCOLORMAP = USERDEF\n")
    f.write("$! COLORMAPCONTROL RESETTOFACTORY\n")
    f.write("$! COLORMAP USERDEFINED {\n")
    f.write(f"  NUMCONTROLPOINTS = {num_colors}\n")

    for i, color in enumerate(colors):
        f.write(f"  CONTROLPOINT {i+1} {{\n")
        f.write(f"    COLORMAPFRACTION = {i/(num_colors-1)}\n")
        f.write(f"    LEADRGB {{ R = {int(color[0]*255)} G = {int(color[1]*255)} B = {int(color[2]*255)} }}\n")
        f.write(f"    TRAILRGB {{ R = {int(color[0]*255)} G = {int(color[1]*255)} B = {int(color[2]*255)} }}\n")
        f.write("  }\n")
    f.write("}\n")
