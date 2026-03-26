import re
import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.colors import Normalize

# ============================================
# READ METADATA (REQUIRED)
# ============================================
def read_run_info(filename="run_info.txt"):
    params = {}
    with open(filename, "r") as f:
        for line in f:
            key, value = line.split()[:2]
            params[key.lower()] = float(value.replace("d", "e").replace("D", "E"))
    return params

params = read_run_info()

dt = params["dt"]
nx = int(params["nx"])
ny = int(params["ny"])


print("Creating gif...")

print(f"Grid: {nx} x {ny}")
print(f"dt = {dt}")

# ============================================
# LOAD FILES
# ============================================
files = sorted(glob.glob("SWE_eta_*.bin"))

if not files:
    raise FileNotFoundError("No SWE_eta_*.bin files found")

print(f"Found {len(files)} frames")

# ============================================
# GRID SETUP
# ============================================
x = np.arange(nx)
y = np.arange(ny)
X, Y = np.meshgrid(x, y, indexing="ij")

stride = 2
Xs = X[::stride, ::stride]
Ys = Y[::stride, ::stride]

# ============================================
# PARSE TIME FROM FILENAMES
# ============================================
def extract_time(filename):
    match = re.search(r"SWE_eta_(\d+)\.bin", filename)
    return int(match.group(1)) * dt if match else None

times = [extract_time(f) for f in files]

# ============================================
# PRE-SCAN FOR CONSISTENT SCALING
# ============================================
zmin, zmax = np.inf, -np.inf

for f in files:
    A = np.fromfile(f, dtype=np.float64).reshape((nx, ny), order="F")
    zmin = min(zmin, A.min())
    zmax = max(zmax, A.max())

norm = Normalize(vmin=zmin, vmax=zmax)

# ============================================
# FIGURE SETUP
# ============================================
fig = plt.figure(figsize=(9, 7), dpi=120)
ax = fig.add_subplot(111, projection="3d")

def update(i):
    ax.clear()

    A = np.fromfile(files[i], dtype=np.float64).reshape((nx, ny), order="F")
    As = A[::stride, ::stride]

    surf = ax.plot_surface(
        Xs, Ys, As,
        cmap="viridis",
        norm=norm,
        linewidth=0,
        antialiased=True
    )

    t = times[i]
    ax.set_title(f"Surface Elevation, t = {t:.3f} s")

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel(r"$\eta$")

    ax.set_xlim(0, nx - 1)
    ax.set_ylim(0, ny - 1)
    ax.set_zlim(zmin, zmax)

    ax.view_init(elev=28, azim=135)
    ax.set_box_aspect((nx, ny, 0.35 * max(nx, ny)))

    return (surf,)

# ============================================
# ANIMATION
# ============================================
anim = FuncAnimation(
    fig,
    update,
    frames=len(files),
    interval=100,
    blit=False
)

anim.save("SWE_eta.gif", writer=PillowWriter(fps=10))
plt.close(fig)

print("Saved SWE_eta.gif")