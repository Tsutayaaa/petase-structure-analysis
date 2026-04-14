from pathlib import Path
import MDAnalysis as mda
from MDAnalysis.analysis.dihedrals import Ramachandran
import matplotlib.pyplot as plt
import numpy as np

# ===== 输入 =====
pdb_path = Path("/Users/shulei/PycharmProjects/Biopython/petase-structure-analysis/ramachandran/5XJU.pdb")

# ===== 读取 =====
u = mda.Universe(str(pdb_path))
protein = u.select_atoms("protein")

rama = Ramachandran(protein).run()
angles = rama.results.angles[0]

phi = angles[:, 0]
psi = angles[:, 1]

mask = np.isfinite(phi) & np.isfinite(psi)
phi = phi[mask]
psi = psi[mask]

# ===== 统一尺度（核心）=====
# ===== 核心参数（针对16 cm海报优化）=====
POINT_SIZE = 42          # 点更大一点（关键）
POINT_RADIUS = np.sqrt(POINT_SIZE)

EDGE_LW = 1           # 点边缘更清晰
SPINE_LW = 2.0          # 坐标轴明显加粗
GRID_LW = 0.8           # 网格略增强（但保持淡）
TICK_LEN = 8
TICK_LW = 1.8

LABEL_SIZE = 10         # 坐标轴字体
TITLE_SIZE = 24
CBAR_SIZE = 16

# ===== matplotlib 参数 =====
plt.rcParams.update({
    "font.size": LABEL_SIZE,
    "axes.linewidth": SPINE_LW,
    "xtick.major.width": TICK_LW,
    "ytick.major.width": TICK_LW,
    "xtick.major.size": TICK_LEN,
    "ytick.major.size": TICK_LEN,
    "svg.fonttype": "none",
})

# ===== 画图 =====
fig = plt.figure(figsize=(7.5, 6.8), dpi=300)

gs = fig.add_gridspec(
    1, 2,
    width_ratios=[1, 0.03],
    wspace=0.08
)

ax = fig.add_subplot(gs[0, 0])
cax = fig.add_subplot(gs[0, 1])

# 背景
rama.plot(ax=ax, ref=True, color="none")

# 点
residue_order = np.arange(len(phi))

sc = ax.scatter(
    phi,
    psi,
    c=residue_order,
    cmap="rainbow",
    s=POINT_SIZE,
    alpha=0.9,
    edgecolors="white",
    linewidths=EDGE_LW,
    zorder=3,
)

# colorbar
cbar = fig.colorbar(sc, cax=cax)
cbar.set_label("Residue order (N → C)", fontsize=CBAR_SIZE)

cbar.outline.set_linewidth(SPINE_LW * 0.8)
cbar.ax.tick_params(
    width=TICK_LW,
    length=TICK_LEN * 0.8,
    labelsize=CBAR_SIZE * 0.8
)

# 坐标
ax.set_xlim(-180, 180)
ax.set_ylim(-180, 180)
ax.set_box_aspect(1)

ax.set_xlabel("φ (phi)", fontsize=LABEL_SIZE)
ax.set_ylabel("ψ (psi)", fontsize=LABEL_SIZE)
ax.set_title("Ramachandran Plot of PETase", fontsize=TITLE_SIZE, pad=12)

ax.set_xticks([-180, -120, -60, 0, 60, 120, 180])
ax.set_yticks([-180, -120, -60, 0, 60, 120, 180])

# 网格
ax.grid(True, linewidth=GRID_LW, alpha=0.18)

for spine in ax.spines.values():
    spine.set_linewidth(SPINE_LW)

# 输出
fig.savefig("petase_ramachandran_scaled.svg", bbox_inches="tight")
fig.savefig("petase_ramachandran_scaled.png", dpi=300, bbox_inches="tight")

plt.show()