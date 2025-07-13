from doctest import debug
import os
import numpy as np
import matplotlib.pyplot as plt
import uproot
import pprint
import sys
import lmfit as lf
import lmfit.models as lfm

import opg_tool

debug_flag = 1

plt.rcParams['font.family'] = 'Times New Roman' #全体のフォントを設定
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.size'] = 28
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

script_dir = os.path.dirname(os.path.abspath(__file__))
root_file_path = os.path.join(script_dir, "../results/root/kvc2_opg.root")

file = uproot.open(root_file_path)
tree = file["tree"].arrays(library="np")

# -- HV dependence -----
run_num_list = [
    # 378, 381,
    # 678,
    778, 779, 780,
    783, 784,

    # 507, 508,
    808, 809,
    812, 813, 814,

    757, 758, 759,
    762
]
opg_data = {
    "ch1": dict(),
    "ch2": dict(),
    "ch3": dict(),
    "ch4": dict()    
}
for i in range(len(tree["run_num"])):
    if tree["run_num"][i] in run_num_list:
        key = f"{tree["hv"][i]}-{tree["seg"][i]}-{tree["ch"][i]}-{17-tree["board_b"][i]}"
        if key in opg_data[f"ch{tree["ch"][i]+1}"]:
            opg_data[f"ch{tree["ch"][i]+1}"][key].append([tree["result_val"][i][3], tree["result_err"][i][3]])
        else:
            opg_data[f"ch{tree["ch"][i]+1}"][key] = [[tree["result_val"][i][3], tree["result_err"][i][3]]]

# -- calc weighted mean for each condition -----
for_scale_data = np.full_like(np.zeros(4*4*2*2), np.nan).reshape(4, 4, 2, 2) # seg, board num, HV, values
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

        if mppc_num in [1, 8, 4]:
            for_scale_data[seg][ch][mppc_hv-57] = [mean, error]

if debug_flag: # for debug
    pprint.pprint(for_scale_data)


fig = plt.figure(figsize=(10, 6))
ax1 = fig.add_subplot(141)
ax2 = fig.add_subplot(142)
ax3 = fig.add_subplot(143)
ax4 = fig.add_subplot(144)
range_min, range_max = np.inf, -np.inf
markers = ["o", "s", "^"]
labels = ["A", "B", "C", "D"]

linear_fit_result = np.full_like(np.zeros(3*4*2), np.nan).reshape(3, 4, 2) # seg, board num, values
for ch, ax in zip(range(4), [ax1, ax2, ax3, ax4]):
    for seg in range(3):
        ax.errorbar(
            [57, 58], for_scale_data[seg][ch][:, 0], yerr = for_scale_data[seg][ch][:, 1], 
            fmt = markers[seg], capsize = 0, markeredgecolor = "k", ms = 8, ecolor='black',  color=f'C{ch}', markeredgewidth = 0.2, zorder = 3, alpha = 0.9, label = f"seg.{seg+1}"
        )
        ax.set_title(f"Board {labels[ch]}", fontsize = 24)

        model = lfm.LinearModel()
        params = model.guess(x = [57, 58], data = for_scale_data[seg][ch][:, 0])
        result = model.fit(x = [57, 58], data = for_scale_data[seg][ch][:, 0], weights = 1.0/for_scale_data[seg][ch][:, 1], params=params, method='leastsq')
        # print(result.fit_report())

        a = result.result.params["slope"].value
        b = result.result.params["intercept"].value
        linear_fit_result[seg][ch] = [a, b]
        
        fit_x = np.linspace(57, 58)
        fit_y = result.eval_components(x=fit_x)["linear"]
        ax.plot(fit_x, fit_y, "--", color=f'C{ch}')
        
        ymin, ymax = ax.get_ylim()
        if ymin < range_min:
            range_min = ymin
        if ymax > range_max:
            range_max = ymax

for ax in [ax1, ax2, ax3, ax4]:
    ax.set_ylim(range_min, range_max)
    ax.set_xlim(56.5, 58.5)
    if ax != ax1:
        ax.set_yticklabels([])

fig.supxlabel("HV [V]", x = (0.1+0.98)/2.0)
ax1.set_ylabel("One Photon Gain [arb. unit]")
plt.subplots_adjust(left = 0.1, right = 0.98, top = 0.9, bottom = 0.15, wspace=0.01)
plt.savefig(os.path.join(script_dir, "../results/img/kvc2_opg_HV_dependence.png"), format='png', bbox_inches='tight', dpi=600, transparent=True)
plt.show()


# -- temperature dependence -----
run_num_list = [
    # 378, 381,
    678,

    507, 508,
]
opg_data_2 = {
    "ch1": dict(),
    "ch2": dict(),
    "ch3": dict(),
    "ch4": dict()
}
for i in range(len(tree["run_num"])):
    if tree["run_num"][i] in run_num_list:
        key = f"{tree["hv"][i]}-{tree["seg"][i]}-{tree["ch"][i]}-{17-tree["board_b"][i]}"
        if key in opg_data_2[f"ch{tree["ch"][i]+1}"]:
            opg_data_2[f"ch{tree["ch"][i]+1}"][key].append([tree["result_val"][i][3], tree["result_err"][i][3]])
        else:
            opg_data_2[f"ch{tree["ch"][i]+1}"][key] = [[tree["result_val"][i][3], tree["result_err"][i][3]]]

# -- calc weighted mean for each condition -----
for_scale_data_2 = np.full_like(np.zeros(4*4*2*2), np.nan).reshape(4, 4, 2, 2) # seg, board num, HV, values
for ch in range(4):
    for key in opg_data_2[f"ch{ch+1}"].keys():
        n        = len(opg_data_2[f"ch{ch+1}"][key])
        mppc_hv  = int(key.split("-")[0])
        seg      = int(key.split("-")[1])
        mppc_num = int(key.split("-")[3]) - 1
        if n > 1:
            opg = np.array(opg_data_2[f"ch{ch+1}"][key])
            exist_outlier = opg_tool.check_condition(opg[:, 0])
            if exist_outlier:
                print(ch, key)
                pprint.pprint(opg[:, 0])
                print("-----")
            mean, error = opg_tool.weighted_mean_and_error(opg[:, 0], 1.0/opg[:, 1]**2)
        else:
            mean, error = opg_data_2[f"ch{ch+1}"][key][0]

        if mppc_num in [1, 8, 4]:
            for_scale_data_2[seg][ch][mppc_hv-57] = [mean, error]

pprint.pprint(for_scale_data_2)


# -- temperature dependence -----
run_num_list = [
    378, 381,
]
opg_data_3 = {
    "ch1": dict(),
    "ch2": dict(),
    "ch3": dict(),
    "ch4": dict()
}
for i in range(len(tree["run_num"])):
    if tree["run_num"][i] in run_num_list:
        key = f"{tree["hv"][i]}-{tree["seg"][i]}-{tree["ch"][i]}-{17-tree["board_b"][i]}"
        if key in opg_data_3[f"ch{tree["ch"][i]+1}"]:
            opg_data_3[f"ch{tree["ch"][i]+1}"][key].append([tree["result_val"][i][3], tree["result_err"][i][3]])
        else:
            opg_data_3[f"ch{tree["ch"][i]+1}"][key] = [[tree["result_val"][i][3], tree["result_err"][i][3]]]

# -- calc weighted mean for each condition -----
for_scale_data_3 = np.full_like(np.zeros(4*4*2*2), np.nan).reshape(4, 4, 2, 2) # seg, board num, HV, values
for ch in range(4):
    for key in opg_data_3[f"ch{ch+1}"].keys():
        n        = len(opg_data_3[f"ch{ch+1}"][key])
        mppc_hv  = int(key.split("-")[0])
        seg      = int(key.split("-")[1])
        mppc_num = int(key.split("-")[3]) - 1
        if n > 1:
            opg = np.array(opg_data_3[f"ch{ch+1}"][key])
            exist_outlier = opg_tool.check_condition(opg[:, 0])
            if exist_outlier:
                print(ch, key)
                pprint.pprint(opg[:, 0])
                print("-----")
            mean, error = opg_tool.weighted_mean_and_error(opg[:, 0], 1.0/opg[:, 1]**2)
        else:
            mean, error = opg_data_3[f"ch{ch+1}"][key][0]

        if mppc_num in [1, 8, 4]:
            for_scale_data_3[seg][ch][mppc_hv-57] = [mean, error]

pprint.pprint(for_scale_data_3)


fig = plt.figure(figsize=(9, 6))
ax  = fig.add_subplot(111)

seg = 1
for ch in range(4):
    ax.plot(ch, for_scale_data[seg][ch][1][0], "o", color = "C0", ms = 8)
    ax.plot(ch, for_scale_data_2[seg][ch][1][0], "s", color = "C1", ms = 8)
    
    diff = for_scale_data[seg][ch][1][0] - for_scale_data_2[seg][ch][1][0]
    # HV_diff = ( for_scale_data_2[ch][2][0] - linear_fit_result[ch][1] ) / linear_fit_result[ch][0] - 58.0
    # print(HV_diff, HV_diff*1000.0/0.4) 
    HV_diff = diff/linear_fit_result[seg][ch][0]
    print(HV_diff*1000.0, HV_diff*1000.0/4.1) 

    ax.plot([], [], " ", label = f"ch{ch+1}: $\Delta$Gain = {diff:.3f}")

ax.set_xlabel("ch")
ax.set_ylabel("One Photon Gain [arb. unit]")
ax.legend(loc='upper left', fontsize = 18, bbox_to_anchor=(0.98, 1.0), handletextpad=-2.0, title = "KVC2 seg.2", title_fontsize = 20)
plt.subplots_adjust(left = 0.15, right = 0.75, top = 0.98, bottom = 0.15)
plt.savefig(os.path.join(script_dir, "../results/img/kvc2_opg_temp_dependence_seg2.png"), format='png', bbox_inches='tight', dpi=600, transparent=True)
plt.show()


fig = plt.figure(figsize=(9, 6))
ax  = fig.add_subplot(111)

seg = 2
for ch in range(4):
    ax.plot(ch, for_scale_data[seg][ch][1][0], "o", color = "C0", ms = 8)
    ax.plot(ch, for_scale_data_2[seg][ch][1][0], "s", color = "C1", ms = 8)
    ax.plot(ch, for_scale_data_3[seg][ch][1][0], "^", color = "C2", ms = 8)
    
    diff = for_scale_data_2[seg][ch][1][0] - for_scale_data[seg][ch][1][0]
    # HV_diff = ( for_scale_data_2[ch][2][0] - linear_fit_result[ch][1] ) / linear_fit_result[ch][0] - 58.0
    # print(HV_diff, HV_diff*1000.0/0.4) 
    HV_diff = diff/linear_fit_result[seg][ch][0]
    print(HV_diff*1000.0, HV_diff*1000.0/2.2)

    ax.plot([], [], " ", label = f"ch{ch+1}: $\Delta$Gain = {diff:.3f}")

ax.set_xlabel("ch")
ax.set_ylabel("One Photon Gain [arb. unit]")
ax.legend(loc='upper left', fontsize = 18, bbox_to_anchor=(0.98, 1.0), handletextpad=-2.0, title = "KVC2 seg.3", title_fontsize = 20)
plt.subplots_adjust(left = 0.15, right = 0.75, top = 0.98, bottom = 0.15)
plt.savefig(os.path.join(script_dir, "../results/img/kvc2_opg_temp_dependence_seg3.png"), format='png', bbox_inches='tight', dpi=600, transparent=True)
plt.show()
