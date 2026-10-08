import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import FancyArrowPatch

fig, ax = plt.subplots(figsize=(11, 6.4), dpi=150)
ax.set_xlim(0, 11)
ax.set_ylim(0, 6.4)
ax.axis('off')

def box(x, y, w, h, text, fc='#eaf2fb', ec='#2b6cb0', fontsize=10.5, weight='normal'):
    rect = patches.FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.06,rounding_size=0.08",
                                   linewidth=1.4, edgecolor=ec, facecolor=fc)
    ax.add_patch(rect)
    ax.text(x + w/2, y + h/2, text, ha='center', va='center', fontsize=fontsize,
             weight=weight)
    return (x, y, w, h)

def straight_arrow(p1, p2, color='#333333', lw=1.3):
    arr = FancyArrowPatch(p1, p2, arrowstyle='-|>', mutation_scale=14,
                           linewidth=lw, color=color, shrinkA=2, shrinkB=2)
    ax.add_patch(arr)

def hbox_arrow(b1, b2):
    x1, y1, w1, h1 = b1
    x2, y2, w2, h2 = b2
    straight_arrow((x1 + w1, y1 + h1/2), (x2, y2 + h2/2))

def vbox_arrow(b1, b2):
    x1, y1, w1, h1 = b1
    x2, y2, w2, h2 = b2
    straight_arrow((x1 + w1/2, y1), (x2 + w2/2, y2 + h2))

# Row 1: transmit chain
b_data   = box(0.3, 5.0, 1.5, 1.0, "Dữ liệu\nký hiệu X", fc='#f5f5f5', ec='#666666')
b_mod    = box(2.15, 5.0, 2.0, 1.0, "Điều chế waveform\nS = M·X\n(OFDM/OTFS/AFDM/RM)", fc='#eaf2fb', ec='#2b6cb0')
b_chan   = box(4.55, 5.0, 1.9, 1.0, "Kênh đa đường/\nDoppler H\n(Phần I Mục 1) + nhiễu", fc='#fdf2e0', ec='#c07a12')
b_demod  = box(6.85, 5.0, 2.0, 1.0, "Giải điều chế\nZ = A·Y\n(Phần I Mục 1)", fc='#eaf2fb', ec='#2b6cb0')

hbox_arrow(b_data, b_mod)
hbox_arrow(b_mod, b_chan)
hbox_arrow(b_chan, b_demod)

# Row 2: two receiver branches, directly below b_demod, side by side
b_eq   = box(6.15, 2.9, 2.9, 1.0, "Cân bằng ZF / LMMSE\nG⁻¹Y  hoặc\n(Heffᴴ Heff+σ²I)⁻¹Heffᴴz\n(Mục 9.2, 9.6)", fc='#eafaf0', ec='#22863a', fontsize=9.6)
b_caf  = box(2.85, 2.9, 2.9, 1.0, "Tương quan CAF /\nmatched filter\n(Mục 9.3, 9.6)", fc='#fdeef3', ec='#c2255c')

demod_bottom = (b_demod[0] + b_demod[2]/2, b_demod[1])
straight_arrow(demod_bottom, (b_eq[0] + b_eq[2]/2, b_eq[1] + b_eq[3]))
straight_arrow(demod_bottom, (b_caf[0] + b_caf[2]/2, b_caf[1] + b_caf[3]))

# Row 3: outputs
b_ber  = box(6.35, 0.9, 2.5, 1.0, "BER/SER\n(Mục 10-11, 13)", fc='#eafaf0', ec='#22863a', fontsize=10.5, weight='bold')
b_sens = box(3.05, 0.9, 2.5, 1.0, "Sai số ước lượng trễ/\nDoppler, cảm biến\n(Mục 12)", fc='#fdeef3', ec='#c2255c', fontsize=10.5, weight='bold')

vbox_arrow(b_eq, b_ber)
vbox_arrow(b_caf, b_sens)

ax.text(5.5, 6.15, "Sơ đồ tổng quan quy trình mô phỏng (Phần II)", ha='center', va='center',
        fontsize=13, weight='bold')

plt.tight_layout()
plt.savefig(r'D:\doubly-dispersive-channel-simulation-main\report_build\system_overview_v1.png', bbox_inches='tight')
print('saved')
