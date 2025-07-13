import os
import numpy as np
import matplotlib.pyplot as plt
import uproot
import pprint
import lmfit as lf
import lmfit.models as lfm

import opg_tool

debug_flag = 0

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
root_file_path = os.path.join(script_dir, "../results/root/bac_opg.root")

file = uproot.open(root_file_path)
tree = file["tree"].arrays(library="np")

# -- HV dependence -----
opg_data = {
    "ch1": dict(),
    "ch2": dict(),
    "ch3": dict(),
    "ch4": dict()
}
for i in range(len(tree["run_num"])):
    if tree["run_num"][i] >= 484:
        key = f"{tree["hv"][i]}-{tree["board_a"][i]}"
        if key in opg_data[f"ch{tree["ch"][i]+1}"]:
            opg_data[f"ch{tree["ch"][i]+1}"][key].append([tree["result_val"][i][3], tree["result_err"][i][3]])
        else:
            opg_data[f"ch{tree["ch"][i]+1}"][key] = [[tree["result_val"][i][3], tree["result_err"][i][3]]]

pprint.pprint(opg_data)

# -- calc weighted mean for each condition -----
for_scale_data = np.zeros(4*3*2).reshape(4, 3, 2)
for ch in range(4):
    for key in opg_data[f"ch{ch+1}"].keys():
        n        = len(opg_data[f"ch{ch+1}"][key])
        mppc_hv  = int(key.split("-")[0])
        mppc_num = int(key.split("-")[1])-1
        if n > 1:
            opg = np.array(opg_data[f"ch{ch+1}"][key])
            mean, error = opg_tool.weighted_mean_and_error(opg[:, 0], 1.0/opg[:, 1]**2)
        else:
            mean, error = opg_data[f"ch{ch+1}"][key][0]

        for_scale_data[ch][mppc_hv-56] = [mean, error]


linear_fit_result = []
fig = plt.figure(figsize=(8, 6))
ax  = fig.add_subplot(111)
for ch in range(4):
    ax.errorbar(
        [56, 57, 58], for_scale_data[ch][:, 0], yerr = for_scale_data[ch][:, 1], 
        fmt = "s", capsize = 0, markeredgecolor = "k", ms = 8, ecolor='black',  color=f'C{ch}', markeredgewidth = 0.2, zorder = 3, alpha = 0.8
    )

    model = lfm.LinearModel()
    params = model.guess(x = [56, 57, 58], data = for_scale_data[ch][:, 0])
    result = model.fit(x = [56, 57, 58], data = for_scale_data[ch][:, 0], weights = 1.0/for_scale_data[ch][:, 1], params=params, method='leastsq')
    # print(result.fit_report())

    a = result.result.params["slope"].value
    b = result.result.params["intercept"].value
    linear_fit_result.append([a, b])
    
    fit_x = np.linspace(56, 58)
    fit_y = result.eval_components(x=fit_x)["linear"]
    ax.plot(fit_x, fit_y, "--", color=f'C{ch}', label = f"ch{ch+1}")
    
    # corr_vol = 58.0 + 0.4*54.0/1000.0
    # print(for_scale_data[ch][2][0], a*corr_vol + b)


ax.set_xlabel("HV [V]")
ax.set_ylabel("One Photon Gain [arb. unit]")
ax.legend(loc='upper left', fontsize = 18, bbox_to_anchor=(1.0, 1), title = f"Temp: 30.4 deg", title_fontsize = 20)
plt.subplots_adjust(left = 0.13, right = 0.72, top = 0.98, bottom = 0.15)
plt.savefig(os.path.join(script_dir, "../results/img/bac_opg_HV_dependence.png"), format='png', bbox_inches='tight', dpi=600, transparent=True)
plt.show()


# -- temperature dependence -----
opg_data_2 = {
    "ch1": dict(),
    "ch2": dict(),
    "ch3": dict(),
    "ch4": dict()
}
for i in range(len(tree["run_num"])):
    if tree["run_num"][i] < 407:
        key = f"{tree["hv"][i]}-{tree["board_a"][i]}"
        if key in opg_data_2[f"ch{tree["ch"][i]+1}"]:
            opg_data_2[f"ch{tree["ch"][i]+1}"][key].append([tree["result_val"][i][3], tree["result_err"][i][3]])
        else:
            opg_data_2[f"ch{tree["ch"][i]+1}"][key] = [[tree["result_val"][i][3], tree["result_err"][i][3]]]

# -- calc weighted mean for each condition -----
for_scale_data_2 = np.zeros(4*3*2).reshape(4, 3, 2)
for ch in range(4):
    for key in opg_data_2[f"ch{ch+1}"].keys():
        n        = len(opg_data_2[f"ch{ch+1}"][key])
        mppc_hv  = int(key.split("-")[0])
        mppc_num = int(key.split("-")[1])-1
        if n > 1:
            opg = np.array(opg_data_2[f"ch{ch+1}"][key])
            mean, error = opg_tool.weighted_mean_and_error(opg[:, 0], 1.0/opg[:, 1]**2)
        else:
            mean, error = opg_data_2[f"ch{ch+1}"][key][0]

        for_scale_data_2[ch][mppc_hv-56] = [mean, error]

pprint.pprint(for_scale_data_2)

fig = plt.figure(figsize=(9, 6))
ax  = fig.add_subplot(111)
for ch in range(4):
    ax.plot(ch, for_scale_data[ch][2][0], "o", color = "C0", ms = 8)
    ax.plot(ch, for_scale_data_2[ch][2][0], "s", color = "C1", ms = 8)
    diff = for_scale_data_2[ch][2][0] - for_scale_data[ch][2][0]
    # HV_diff = ( for_scale_data_2[ch][2][0] - linear_fit_result[ch][1] ) / linear_fit_result[ch][0] - 58.0
    # print(HV_diff, HV_diff*1000.0/0.4) 
    HV_diff = diff/linear_fit_result[ch][0]
    print(HV_diff, HV_diff*1000.0/0.4) 

    ax.plot([], [], " ", label = f"ch{ch+1}: $\Delta$Gain = {diff:.3f}")

ax.set_xlabel("ch")
ax.set_ylabel("One Photon Gain [arb. unit]")
ax.legend(loc='upper left', fontsize = 18, bbox_to_anchor=(0.98, 1.0), handletextpad=-2.0)
plt.subplots_adjust(left = 0.15, right = 0.75, top = 0.98, bottom = 0.15)
plt.savefig(os.path.join(script_dir, "../results/img/bac_opg_temp_dependence.png"), format='png', bbox_inches='tight', dpi=600, transparent=True)
plt.show()