# SPT_BH_Radiation
Space Pressure Tensor model for Black Hole Radiation

# SPT_BH_Radiation
Space Pressure Tensor model for unifying Hawking radiation and superradiance in black holes.
- **Paper**: [Final Draft](#) (link to be added)
- **Figures**: 3D visualizations (fig1_schwarzschild.py, fig2_kerr.py, fig3_evaporating_kerr.py)
- **GIF**: Time evolution (gif_evolution.py)
- **Requirements**: numpy, plotly, matplotlib, pillow
- **Author**: Junichi-o
- **Date**: June 15, 2025

- - **コード**: [GitHub: SPT_BH_Radiation](https://github.com/Junichi-o/SPT_BH_Radiation)（2025年6月17日アップ予定）

fig1_schwarzschild.py
import numpy as np
import plotly.graph_objects as go
M0 = 1.0; epsilon, delta = 0.01, 0.01; n = 1.5; spt_alpha, s_osc = 0.1, 1e24
r = np.linspace(2, 10*M0, 100); theta = np.linspace(0, np.pi, 50)
r_grid, theta_grid = np.meshgrid(r, theta); T_H = 1 / (8 * np.pi * M0)
omega = 0.111 / M0 * (1 + spt_alpha * np.cos(2 * np.pi * r_grid / s_osc))
f_r = epsilon * (np.sin(omega * r_grid) / r_grid**n)
rho = (epsilon * ((2 - n) * r_grid**(1 - n) * np.sin(omega * r_grid) + omega * r_grid**(2 - n) * np.cos(omega * r_grid)))**2
x, y, z = r_grid * np.sin(theta_grid), r_grid * np.cos(theta_grid), rho
fig = go.Figure(data=[go.Surface(x=x, y=y, z=z)])
fig.update_layout(title='シュワルツシルトエネルギー密度', scene=dict(xaxis_title='x/M0', yaxis_title='y/M0', zaxis_title='エネルギー密度'))
fig.show()

fig2_kerr.pyimport numpy as np
import plotly.graph_objects as go
M0, alpha_spin0 = 1.0, 0.6; epsilon, delta = 0.01, 0.01; beta, gamma, n = 1.6, 0.06, 1.5; spt_alpha, s_osc = 0.1, 1e24
r = np.linspace(1.8, 10*M0, 100); theta = np.linspace(0, np.pi, 50)
r_grid, theta_grid = np.meshgrid(r, theta); r_plus = M0 + np.sqrt(M0**2 - (alpha_spin0 * M0)**2)
T_H = np.sqrt(M0**2 - (alpha_spin0 * M0)**2) / (4 * np.pi * M0 * (M0 + np.sqrt(M0**2 - (alpha_spin0 * M0)**2)))
omega = 0.111 / M0 * (1 + spt_alpha * np.cos(2 * np.pi * r_grid / s_osc))
f_r = epsilon * (np.sin(omega * r_grid) / r_grid**n) * (1 + beta * np.cos(theta_grid)**2)
div_P_r = epsilon * ((2 - n) * r_grid**(1 - n) * np.sin(omega * r_grid) + omega * r_grid**(2 - n) * np.cos(omega * r_grid)) * (1 + beta * np.cos(theta_grid)**2) - 2 * delta * np.cos(omega * r_grid) * (1 + beta * np.cos(theta_grid)**2) / r_grid**(n + 2) + gamma * (alpha_spin0 * M0 * np.sin(theta_grid)**2 / r_grid**(n + 1)) * np.cos(omega * r_grid)
rho = div_P_r**2; x, y, z = r_grid * np.sin(theta_grid), r_grid * np.cos(theta_grid), rho
fig = go.Figure(data=[go.Surface(x=x, y=y, z=z)])
fig.update_layout(title='カーエネルギー密度', scene=dict(xaxis_title='x/M0', yaxis_title='y/M0', zaxis_title='エネルギー密度'))
fig.show()

fig3_evaporating_kerr.py
import numpy as np
import plotly.graph_objects as go
M0, alpha_spin0 = 1.0, 0.6; kappa, eta = 1e-4, 2; epsilon, delta = 0.01, 0.01; beta, gamma, n = 1.6, 0.06, 1.5; spt_alpha = 0.1; s_osc = 1e24
r = np.linspace(1.8, 10*M0, 100); theta = np.linspace(0, np.pi, 50)
r_grid, theta_grid = np.meshgrid(r, theta); times = [0, 5000, 10000]
fig = go.Figure()
for t in times:
    M = M0 / (1 + 3 * kappa * t / M0**2)**(1/3); alpha_spin = alpha_spin0 * (M / M0)**eta
    r_plus = M + np.sqrt(M**2 - (alpha_spin * M)**2)
    T_H = np.sqrt(M**2 - (alpha_spin * M)**2) / (4 * np.pi * M * (M + np.sqrt(M**2 - (alpha_spin * M)**2)))
    omega = 0.111 / M * (1 + spt_alpha * np.cos(2 * np.pi * r_grid / s_osc))
    f_r = epsilon * (np.sin(omega * r_grid) / r_grid**n) * (1 + beta * np.cos(theta_grid)**2)
    div_P_r = epsilon * ((2 - n) * r_grid**(1 - n) * np.sin(omega * r_grid) + omega * r_grid**(2 - n) * np.cos(omega * r_grid)) * (1 + beta * np.cos(theta_grid)**2) - 2 * delta * np.cos(omega * r_grid) * (1 + beta * np.cos(theta_grid)**2) / r_grid**(n + 2) + gamma * (alpha_spin * M * np.sin(theta_grid)**2 / r_grid**(n + 1)) * np.cos(omega * r_grid)
    rho = div_P_r**2; x, y, z = r_grid * np.sin(theta_grid), r_grid * np.cos(theta_grid), rho
    fig.add_trace(go.Surface(x=x, y=y, z=z, name=f't={t}'))
fig.update_layout(title='蒸発カーBHエネルギー密度', scene=dict(xaxis_title='x/M0', yaxis_title='y/M0', zaxis_title='エネルギー密度'))
fig.show()

gif_evolution.py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
M0, alpha_spin0 = 1.0, 0.6; kappa, eta = 1e-4, 2; epsilon, delta = 0.01, 0.01; beta, gamma, n = 1.6, 0.06, 1.5; spt_alpha = 0.1; s_osc = 1e24
r = np.linspace(1.8, 10*M0, 100); theta = np.linspace(0, np.pi, 50)
r_grid, theta_grid = np.meshgrid(r, theta)
fig, ax = plt.subplots(figsize=(8, 5))
def update(t): ax.clear(); M = M0 / (1 + 3 * kappa * t / M0**2)**(1/3); alpha_spin = alpha_spin0 * (M / M0)**eta
    r_plus = M + np.sqrt(M**2 - (alpha_spin * M)**2)
    T_H = np.sqrt(M**2 - (alpha_spin * M)**2) / (4 * np.pi * M * (M + np.sqrt(M**2 - (alpha_spin * M)**2)))
    omega = 0.111 / M * (1 + spt_alpha * np.cos(2 * np.pi * r_grid / s_osc))
    f_r = epsilon * (np.sin(omega * r_grid) / r_grid**n) * (1 + beta * np.cos(theta_grid)**2)
    div_P_r = epsilon * ((2 - n) * r_grid**(1 - n) * np.sin(omega * r_grid) + omega * r_grid**(2 - n) * np.cos(omega * r_grid)) * (1 + beta * np.cos(theta_grid)**2) - 2 * delta * np.cos(omega * r_grid) * (1 + beta * np.cos(theta_grid)**2) / r_grid**(n + 2) + gamma * (alpha_spin * M * np.sin(theta_grid)**2 / r_grid**(n + 1)) * np.cos(omega * r_grid)
    rho = div_P_r**2; rho_hawking = (T_H)**4 / (r_grid / r_plus)**6
    ax.contourf(r / M0, theta, rho, levels=50, cmap='viridis')
    ax.contour(r / M0, theta, rho_hawking, levels=10, colors='red', linestyles='--')
    ax.set_xlabel('r / M0'); ax.set_ylabel('θ (rad)'); ax.set_title(f't = {int(t)}')
ani = FuncAnimation(fig, update, frames=np.linspace(0, 10000, 50), interval=200)
ani.save('/content/bh_evolution.gif', writer='pillow'); plt.close()
import os; from google.colab import files
if os.path.exists('/content/bh_evolution.gif'): files.download('/content/bh_evolution.gif'); print("ダウンロード開始！")
else: print("ファイルなし。再実行を！")

requirements.txt
numpy
plotly
matplotlib
pillow

