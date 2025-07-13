from doctest import debug
import os
import numpy as np
import matplotlib.pyplot as plt
import uproot
import pprint
import sys

import opg_tool

debug_flag = 1

plt.rcParams['font.family'] = 'Times New Roman' #全体のフォントを設定
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.size'] = 28
plt.rcParams['axes.linewidth'] = 1.0# 軸の線幅edge linewidth。囲みの太さ
plt.rcParams['axes.grid'] = False
plt.rcParams['axes.axisbelow'] = True
plt.rcParams["xtick.direction"] = "in"               #x軸の目盛線を内向きへ
plt.rcParams["ytick.direction"] = "in"               #y軸の目盛線を内向きへ
plt.rcParams["xtick.minor.visible"] = True           #x軸補助目盛りの追加
plt.rcParams["ytick.minor.visible"] = True           #y軸補助目盛りの追加
plt.rcParams["xtick.major.size"] = 10                #x軸主目盛り線の長さ
plt.rcParams["ytick.major.size"] = 10                #y軸主目盛り線の長さ
plt.rcParams["xtick.minor.size"] = 5                 #x軸補助目盛り線の長さ
plt.rcParams["ytick.minor.size"] = 5                 #y軸補助目盛り線の長さ

script_dir = os.path.dirname(os.path.abspath(__file__))
root_file_path = os.path.join(script_dir, "../results/root/kvc2_opg.root")

file = uproot.open(root_file_path)
tree = file["tree"].arrays(library="np")

# -- distribute data to dict -----
opg_data = {
    "ch1": dict(),
    "ch2": dict(),
    "ch3": dict(),
    "ch4": dict()    
}
for i in range(len(tree["run_num"])):
    if tree["run_num"][i] > 734:
        key = f"{tree["hv"][i]}-{tree["seg"][i]}-{tree["ch"][i]}-{17-tree["board_b"][i]}"
        if key in opg_data[f"ch{tree["ch"][i]+1}"]:
            opg_data[f"ch{tree["ch"][i]+1}"][key].append([tree["result_val"][i][3], tree["result_err"][i][3]])
        else:
            opg_data[f"ch{tree["ch"][i]+1}"][key] = [[tree["result_val"][i][3], tree["result_err"][i][3]]]

# -- calc weighted mean for each condition -----
opg_mean_data  = np.full_like(np.zeros(4*16*2), np.nan).reshape(4, 16, 2)
for_scale_data = np.full_like(np.zeros(4*4*3*2), np.nan).reshape(4, 4, 3, 2) # seg, board num, HV, values
for ch in range(4):
    for key in opg_data[f"ch{ch+1}"].keys():
        n        = len(opg_data[f"ch{ch+1}"][key])
        mppc_hv  = int(key.split("-")[0])
        seg      = int(key.split("-")[1])
        mppc_num = int(key.split("-")[3]) - 1
        if n > 1:
            opg = np.array(opg_data[f"ch{ch+1}"][key])
            exist_outlier = opg_tool.check_condition(opg[:, 0])
            if exist_outlier:
                print(ch, key)
                pprint.pprint(opg[:, 0])
                print("-----")
            mean, error = opg_tool.weighted_mean_and_error(opg[:, 0], 1.0/opg[:, 1]**2)
        else:
            mean, error = opg_data[f"ch{ch+1}"][key][0]

        if mppc_hv == 58:
            opg_mean_data[ch][mppc_num] = [mean, error]
        if mppc_num in [1, 8, 4]:
            for_scale_data[seg][ch][mppc_hv-56] = [mean, error]

if debug_flag: # for debug
    pprint.pprint(opg_mean_data)
    pprint.pprint(for_scale_data)

# sys.exit()

for HV in [57, 58]:
    fig = plt.figure(figsize=(11, 6))
    ax  = fig.add_subplot(111)
    legend_labels = []
    for ch in range(4):
        label_str = ""
        for seg in range(4):
            if HV == 58:
                scale = 1.0
            else:
                scale = for_scale_data[seg][ch][HV-56][0]/for_scale_data[seg][ch][2][0]

            corrected_data = scale * opg_mean_data[ch][4*seg:4*(seg+1)]
            ax.errorbar(
                np.arange(4*seg+1, 4*(seg+1)+1), corrected_data[:, 0], yerr = corrected_data[:, 1], 
                fmt = "s", capsize = 0, markeredgecolor = "k", ms = 8, ecolor='black',  color=f'C{ch}', markeredgewidth = 0.2, zorder = 3
            )

            try:
                mean, error = opg_tool.weighted_mean_and_error(corrected_data[:, 0], 1.0/corrected_data[:, 1]**2)
                if seg == 0:
                    label_str += f"seg.{seg+1}: mean = {mean:.2f}"
                else:
                    label_str += f"\nseg.{seg+1}: mean = {mean:.2f}"
                print(ch, seg, format(mean, ".3f"), format(error, ".3f"))

                # 重み付き平均の横線を描画
                ax.hlines(mean, 4*seg+0.5, 4*(seg+1)+0.5, color=f'C{ch}', linestyle='--')

                # エラー帯を描画
                ax.fill_between(
                    x=[4*seg+0.5, 4*(seg+1)+0.5],  # x範囲をデータ範囲に拡張
                    y1=mean - error,             # エラー範囲下限
                    y2=mean + error,             # エラー範囲上限
                    color=f'C{ch}',
                    alpha=0.2
                )
            except:
                pass

        legend_labels.append(label_str)

    ax.set_xticks(np.arange(1, 16, 2))
    xtick_labels = [f"{i+1}\n{16-i}" for i in range(0, 16, 2)]
    ax.set_xticklabels(xtick_labels)
    ax.set_xlim(0, 17)
    ax.set_xlabel("MPPC number")
    ax.text(-0.17, -0.057, 'MPPC (a, d)', fontsize = 18, transform=ax.transAxes)
    ax.text(-0.17, -0.143, 'MPPC (b, c)', fontsize = 18, transform=ax.transAxes)

    ax.set_ylabel("One Photon Gain [arb. unit]")
    for i, label in enumerate(legend_labels):
        ax.plot([], [], color=f'C{i}', linestyle='--', label=label)  # ダミー線を作成して凡例に追加
    ax.legend(loc='upper left', fontsize = 18, bbox_to_anchor=(1.0, 1), title = f"KVC2 HV = {HV} V", title_fontsize = 20)
    plt.subplots_adjust(left = 0.13, right = 0.72, top = 0.98, bottom = 0.15)
    plt.savefig(os.path.join(script_dir, f"../results/img/kvc2_opg_{HV:.0f}.png"), format='png', bbox_inches='tight', dpi=600, transparent=True)
    plt.savefig(os.path.join(script_dir, f"../results/img/kvc2_opg_{HV:.0f}.jpg"), format='jpg', bbox_inches='tight', dpi=600, transparent=True)
    plt.show()
    
    