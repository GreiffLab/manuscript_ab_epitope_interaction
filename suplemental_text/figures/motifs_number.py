#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import math
import matplotlib.pyplot as plt
import seaborn as sns
from palettable.wesanderson import FantasticFox1_5 as fFox


# In[2]:


max_n = 250
C = [[1, 0]]
for n in range(1, max_n):
    C.append([1])
    for k in range(1, n+1):
        C[-1].append(C[n-1][k] + C[n-1][k-1])
    C[-1].append(0)


# In[3]:


def N(L, A):
    res = 0
    if A < L:
        return 1
    if A == L:
        res += 1
    for nx in range(int(L / 2 + 1), L+1):
        ng = L - nx
        res += C[nx - 1][ng] * C[A - nx - 1][ng - 1]
    return res
def NS(L, A):
    return sum([N(L, A1) for A1 in range(L, A+1)])


# In[4]:


def draw_curve(ax, label, xs, ys, color, text_x, text_y, ha="center", va="top"):
    maximum_id = np.argmax(ys)
    maximum_x = xs[maximum_id]
    maximum_y = ys[maximum_id]
    ax.annotate(s="",
                xy=[maximum_x, maximum_y], xycoords="data",
                xytext=[text_x, text_y], textcoords="axes fraction",
                ha="center",
                arrowprops={"width":0.5,
                            "headwidth":0,
                            "color":"grey",
                            "linestyle":"-"})
    ax.plot(xs,
            ys,
            color=color,
            alpha=0.8,
            marker='o',
            markersize=fontsize / 2,
            label=label,
            linewidth=4)
    ax.annotate(s=f"{maximum_y}",
                xy=[maximum_x, maximum_y], xycoords="data",
                xytext=[text_x, text_y], textcoords="axes fraction",
                ha=ha,
                va=va,
                color=color,
                fontsize=fontsize)
def set_ax_params(ax, vmax, title, legend_loc="best"):
    ax.set_ylim(1, vmax)
    ax.set_yscale("log")
    ax.legend(fontsize=fontsize * 0.8, markerscale=1, loc=legend_loc)
    ax.tick_params(axis="x", labelsize=fontsize * 0.8)
    ax.tick_params(axis="y", labelsize=fontsize * 0.8)
    ax.set_title(title, fontsize=fontsize)


# In[ ]:





# In[5]:


sns.set_style("white")

xs = np.arange(1, 11)


L_As = [(23, 24), (11, 17), (15, 16), (7, 12), (32, 40), (9, 15), (13, 13)]
H_As = [(30, 34), (5, 7), (14, 15), (16, 19), (29, 42), (7, 18), (11, 11)]
L_labels = ["LFR1", "CDR-L1", "LFR2", "CDR-L2", "LFR3", "CDR-L3", "LFR4"]
H_labels = ["HFR1", "CDR-H1", "HFR2", "CDR-H2", "HFR3", "CDR-H3", "HFR4"]

n_region_types = len(L_As)

total_max = NS(10, 42)
fontsize=30

all_ys_min = []
all_ys_max = []

ys_dict_min = {}
ys_dict_max = {}

plt.ylim = (0, total_max)

fig = plt.figure(figsize = (56, 14))

common_ax = fig.add_subplot(1, 1, 1)
common_ax.axis("off")

common_ax.text(s="Motif length", x=0.45, y=-0.15, fontsize=fontsize * 3)
common_ax.text(s="Number of motifs", x=-0.04, y=0.06, fontsize=fontsize * 3, rotation=90)

cmap = plt.cm.seismic
colors = [cmap(0.4), cmap(0.1)]

for i in range(n_region_types):
    ax = fig.add_subplot(2, n_region_types, i+1)
    Amin, Amax = L_As[i]
    ys_min = [NS(L, Amin) for L in xs]
    ys_max = [NS(L, Amax) for L in xs]
    all_ys_min.append(ys_min)
    all_ys_max.append(ys_max)
    ys_dict_min[L_labels[i]] = ys_min
    ys_dict_max[L_labels[i]] = ys_max
    title = L_labels[i]


    max_text_x, max_text_y = [0.86, 0.1]
    if title == "CDR-L3":
        max_text_y += 0.35
    draw_curve(ax, label=f"A_max = {Amax}",
               xs=xs, ys=ys_max, color=colors[1],
               text_x=max_text_x, text_y=max_text_y)
    
    min_text_x, min_text_y = [0.72, 0.25]
    draw_curve(ax, label=f"A_min = {Amin}",
               xs=xs, ys=ys_min, color=colors[0],
               text_x=min_text_x, text_y=min_text_y)
    set_ax_params(ax, vmax=total_max, title=title)
    

    
colors = [cmap(0.6), cmap(0.9)]
for i in range(n_region_types):
    ax = fig.add_subplot(2, n_region_types, n_region_types + i+1)
    Amin, Amax = H_As[i]
    ys_min = [NS(L, Amin) for L in xs]
    ys_max = [NS(L, Amax) for L in xs]
    all_ys_min.append(ys_min)
    all_ys_max.append(ys_max)
    ys_dict_min[H_labels[i]] = ys_min
    ys_dict_max[H_labels[i]] = ys_max
    title = H_labels[i]
    
    max_text_x, max_text_y = [0.86, 0.1]
    draw_curve(ax, label=f"A_max = {Amax}",
               xs=xs, ys=ys_max, color=colors[1],
               text_x=max_text_x, text_y=max_text_y)
    
    min_text_x, min_text_y = [0.72, 0.25]
    if title in ("CDR-H1", "CDR-H3"):
        min_text_x = 0.41
        min_text_y = 0.075
    draw_curve(ax, label=f"A_min = {Amin}",
               xs=xs, ys=ys_min, color=colors[0],
               text_x=min_text_x, text_y=min_text_y)
    
    
    set_ax_params(ax, vmax=total_max, title=title)
    
plt.tight_layout()
plt.show()


# In[6]:


ys_min_sum = np.sum(np.array(all_ys_min), axis=0)
ys_max_sum = np.sum(np.array(all_ys_max), axis=0)


# In[7]:


ys_min_frl = np.sum(np.array([value for key, value in ys_dict_min.items() if key.startswith("LFR")]), axis=0)
ys_max_frl = np.sum(np.array([value for key, value in ys_dict_max.items() if key.startswith("LFR")]), axis=0)

ys_min_cdrl = np.sum(np.array([value for key, value in ys_dict_min.items() if key.startswith("CDR-L")]), axis=0)
ys_max_cdrl = np.sum(np.array([value for key, value in ys_dict_max.items() if key.startswith("CDR-L")]), axis=0)

ys_min_frh = np.sum(np.array([value for key, value in ys_dict_min.items() if key.startswith("HFR")]), axis=0)
ys_max_frh = np.sum(np.array([value for key, value in ys_dict_max.items() if key.startswith("HFR")]), axis=0)

ys_min_cdrh = np.sum(np.array([value for key, value in ys_dict_min.items() if key.startswith("CDR-H")]), axis=0)
ys_max_cdrh = np.sum(np.array([value for key, value in ys_dict_max.items() if key.startswith("CDR-H")]), axis=0)


# In[8]:


vmax = ys_max_sum.max() * 10**0.2

markersize=15
linewidth=3
fontsize=30
fig = plt.figure(figsize=(30, 10))

common_ax = fig.add_subplot(1, 1, 1)
common_ax.axis("off")
common_ax.text(s="Motif length", x=0.45, y=-0.15, fontsize=fontsize * 2)
common_ax.text(s="Number of motifs", x=-0.05, y=0.1, fontsize=fontsize * 2, rotation=90)

cmap = plt.cm.seismic
colors = [cmap(0.1), cmap(0.4), cmap(0.9), cmap(0.6)]


# CDR
ax = fig.add_subplot(1, 3, 1)
ax.set_ylim(1, vmax)

max_text_x_1, max_text_y_1 = [0.7, 0.9]
min_text_x_1, min_text_y_1 = [0.6, 0.8]
max_text_x_2, max_text_y_2 = [0.90, 0.48]
min_text_x_2, min_text_y_2 = [0.7, 0.3]

ha_1, va_1 = "right", "center"
ha_2, va_2 = "center", "top"


draw_curve(ax=ax, label=f"CDRH max",
           xs=xs, ys=ys_max_cdrh, color=colors[2],
           text_x=max_text_x_1,
           text_y=max_text_y_1,
           ha=ha_1, va=va_1)
draw_curve(ax=ax, label=f"CDRH min",
           xs=xs, ys=ys_min_cdrh, color=colors[3],
           text_x=min_text_x_1,
           text_y=min_text_y_1,
           ha=ha_1, va=va_1)

draw_curve(ax=ax, label=f"CDRL max",
           xs=xs, ys=ys_max_cdrl, color=colors[0],
           text_x=max_text_x_2,
           text_y=max_text_y_2,
           ha=ha_2, va=va_2)
draw_curve(ax=ax, label=f"CDRL min",
           xs=xs, ys=ys_min_cdrl, color=colors[1],
           text_x=min_text_x_2,
           text_y=min_text_y_2,
           ha=ha_2, va=va_2)
set_ax_params(ax, vmax=vmax, title="(a) Number of motifs across all CDRs (H/L)", legend_loc="upper left")


# FR
ax = fig.add_subplot(1, 3, 2)
ax.set_ylim(1, vmax)

draw_curve(ax=ax, label=f"HFR max",
           xs=xs, ys=ys_max_frh, color=colors[2],
           text_x=max_text_x_1,
           text_y=max_text_y_1,
           ha=ha_1, va=va_1)
draw_curve(ax=ax, label=f"HFR min",
           xs=xs, ys=ys_min_frh, color=colors[3],
           text_x=min_text_x_1,
           text_y=min_text_y_1,
           ha=ha_1, va=va_1)

draw_curve(ax=ax, label=f"LFR max",
           xs=xs, ys=ys_max_frl, color=colors[0],
           text_x=max_text_x_2,
           text_y=max_text_y_2,
           ha=ha_2, va=va_2)
draw_curve(ax=ax, label=f"LFR min",
           xs=xs, ys=ys_min_frl, color=colors[1],
           text_x=min_text_x_2,
           text_y=min_text_y_2,
           ha=ha_2, va=va_2)
set_ax_params(ax, vmax=vmax, title="(b) Number of motifs across all FRs (H/L)", legend_loc="upper left")

# total
ax = fig.add_subplot(1, 3, 3)
ax.set_ylim(1, vmax)
colors = [np.array(t) / 255 for t in fFox.colors]

draw_curve(ax=ax, label=f"total max",
           xs=xs, ys=ys_max_sum, color=colors[1],
           text_x=max_text_x_1,
           text_y=max_text_y_1,
           ha=ha_1, va=va_1)
draw_curve(ax=ax, label=f"total min",
           xs=xs, ys=ys_min_sum, color=colors[0],
           text_x=max_text_x_2,
           text_y=max_text_y_2,
           ha=ha_2, va=va_2)

set_ax_params(ax, vmax=vmax, title="(c) Total number of motifs across all regions", legend_loc="upper left")


plt.tight_layout()
plt.show()


# In[ ]:




