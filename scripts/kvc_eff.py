import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import uproot
import os

plt.rcParams['font.family'] = 'Times New Roman' #全体のフォントを設定
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.size'] = 18
plt.rcParams['axes.linewidth'] = 1.0# 軸の線幅edge linewidth。囲みの太さ
plt.rcParams['axes.grid'] = True
plt.rcParams['axes.axisbelow'] = True
plt.rcParams["xtick.direction"] = "in"               #x軸の目盛線を内向きへ
plt.rcParams["ytick.direction"] = "in"               #y軸の目盛線を内向きへ
plt.rcParams["xtick.minor.visible"] = True           #x軸補助目盛りの追加
plt.rcParams["ytick.minor.visible"] = True           #y軸補助目盛りの追加
plt.rcParams["xtick.major.size"] = 10                #x軸主目盛り線の長さ
plt.rcParams["ytick.major.size"] = 10                #y軸主目盛り線の長さ
plt.rcParams["xtick.minor.size"] = 5                 #x軸補助目盛り線の長さ
plt.rcParams["ytick.minor.size"] = 5                 #y軸補助目盛り線の長さ

# ROOTファイルを開く
script_dir = os.path.dirname(os.path.abspath(__file__))
root_file_path = os.path.join(script_dir, "../results/root/kvc2_eff.root")

file = uproot.open(root_file_path)
tree = file["tree"].arrays(library="np")

# 軸情報とデータ配列を準備
width_vals = np.arange(120.0, 60.0, -120.0/20.0)   # x軸の中心位置
height_vals = np.arange(26.0, 13.0, -26.0/20.0)

ch = 2

efficiency = np.zeros((10, 10))  # shape: (width_index, height_index)
for w_i in range(10):
    for h_i in range(10):
        index = w_i * 10 + h_i
        efficiency[w_i][h_i] = tree["eff_val"][index][ch]  # [2] → そのrunの全ch中のどれか？

# pandas DataFrame に変換（seabornのheatmapに渡しやすい）
df_eff = pd.DataFrame(efficiency.T,  # 転置して高さがindex, 幅がcolumnsになるように
                      index=np.round(height_vals, 1), 
                      columns=np.round(width_vals, 1))

# ヒートマップを描画
fig, ax = plt.subplots(figsize=(8, 8))
sns.heatmap(df_eff,
            cmap='viridis',
            annot=True,
            fmt=".3f",
            annot_kws={"fontsize": 12},
            cbar_kws=dict(label="Efficiency"),
            ax=ax)

ax.set_xlabel("Width")
ax.set_ylabel("Height")
ax.set_title(f"KVC2 Efficiency (Seg. {ch+1})")
plt.savefig(os.path.join(script_dir, f"../results/img/kvc2_eff{ch}.png"), format='png', bbox_inches='tight', dpi=600, transparent=True)

plt.tight_layout()
plt.show()
