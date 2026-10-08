import sys, math
sys.path.insert(0, "D:/doubly-dispersive-channel-simulation-main/lisat_check")
import torch
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from waveforms import (
    OFDMModulator, OFDMDemodulator,
    OCDMModulator, OCDMDemodulator,
    AFDMModulator, AFDMDemodulator,
    OTFSModulator, OTFSDemodulator,
    qam_mapper, qam_demapper, compute_papr,
)

torch.manual_seed(42)
np.random.seed(42)

N = 64
N_cp = 16
ORDER = 16
BITS_PER_SYM = int(math.log2(ORDER))

# ---------------- Smoke test (reconstruction under no noise) ----------------
print("=== Smoke test (reconstruction, no noise) ===")
batch = 8
bits = torch.randint(0, 2, (batch, N * BITS_PER_SYM))
syms = qam_mapper(bits, order=ORDER)

ofdm_mod, ofdm_dem = OFDMModulator(N=N, N_cp=N_cp), OFDMDemodulator(N=N, N_cp=N_cp)
ocdm_mod, ocdm_dem = OCDMModulator(N=N, N_cp=N_cp), OCDMDemodulator(N=N, N_cp=N_cp)
afdm_mod, afdm_dem = AFDMModulator(N=N, c1=0.25, c2=0.0, N_cp=N_cp), AFDMDemodulator(N=N, c1=0.25, c2=0.0, N_cp=N_cp)

for name, mod, dem in [("OFDM", ofdm_mod, ofdm_dem), ("OCDM", ocdm_mod, ocdm_dem), ("AFDM", afdm_mod, afdm_dem)]:
    sig = mod.modulate(syms)
    rec = dem.demodulate(sig)
    err = (rec - syms).abs().mean().item()
    papr = compute_papr(sig).mean().item()
    print(f"{name}: recon_err={err:.2e}  PAPR={papr:.2f} dB")

M_dd, N_dd = 8, 8  # M_dd*N_dd = 64, matches N; dd_symbols shape is [batch, N, M]
N_cp_otfs = 4  # must be < M_dd (OTFS CP is per delay-bin block, unlike OFDM/AFDM's frame-level CP)
otfs_mod, otfs_dem = OTFSModulator(M=M_dd, N=N_dd, N_cp=N_cp_otfs), OTFSDemodulator(M=M_dd, N=N_dd, N_cp=N_cp_otfs)
dd_syms = syms.reshape(batch, N_dd, M_dd)
sig_t = otfs_mod.modulate(dd_syms)
rec_t = otfs_dem.demodulate(sig_t)
err_t = (rec_t - dd_syms).abs().mean().item()
papr_t = compute_papr(sig_t).mean().item()
print(f"OTFS: recon_err={err_t:.2e}  PAPR={papr_t:.2f} dB")

# ---------------- PAPR comparison (many trials) ----------------
print("\n=== PAPR comparison (2000 random 16-QAM frames) ===")
n_trials = 2000
bits_p = torch.randint(0, 2, (n_trials, N * BITS_PER_SYM))
syms_p = qam_mapper(bits_p, order=ORDER)

papr_vals = {}
papr_vals["OFDM"] = compute_papr(ofdm_mod.modulate(syms_p)).numpy()
papr_vals["OCDM"] = compute_papr(ocdm_mod.modulate(syms_p)).numpy()
papr_vals["AFDM"] = compute_papr(afdm_mod.modulate(syms_p)).numpy()
dd_syms_p = syms_p.reshape(n_trials, N_dd, M_dd)
papr_vals["OTFS"] = compute_papr(otfs_mod.modulate(dd_syms_p)).numpy()

for k, v in papr_vals.items():
    print(f"{k}: mean PAPR={v.mean():.2f} dB, 99th pct={np.percentile(v,99):.2f} dB")

# ---------------- AWGN BER/SER vs SNR ----------------
print("\n=== AWGN BER/SER vs SNR ===")
snr_db_range = np.arange(-5, 20.1, 2.5)
n_mc = 400

def run_ber(mod, dem, is_otfs=False):
    ber = np.zeros(len(snr_db_range))
    ser = np.zeros(len(snr_db_range))
    for i, snr_db in enumerate(snr_db_range):
        bit_errs = 0; sym_errs = 0; tot_bits = 0; tot_syms = 0
        for _ in range(n_mc):
            bits_ = torch.randint(0, 2, (1, N * BITS_PER_SYM))
            syms_ = qam_mapper(bits_, order=ORDER)
            if is_otfs:
                dd_ = syms_.reshape(1, N_dd, M_dd)
                sig_ = mod.modulate(dd_)
            else:
                sig_ = mod.modulate(syms_)
            sig_power = sig_.abs().pow(2).mean()
            noise_var = sig_power / (10 ** (snr_db / 10))
            noise = torch.sqrt(noise_var / 2) * (torch.randn_like(sig_.real) + 1j * torch.randn_like(sig_.imag))
            rx = sig_ + noise
            if is_otfs:
                rec_ = dem.demodulate(rx).reshape(1, N)
                true_ = dd_.reshape(1, N)
            else:
                rec_ = dem.demodulate(rx)
                true_ = syms_
            bits_hat = qam_demapper(rec_, order=ORDER, hard=True)
            bits_true = qam_demapper(true_, order=ORDER, hard=True)
            bit_errs += (bits_hat != bits_true).sum().item()
            tot_bits += bits_true.numel()
            # symbol error: nearest-constellation mismatch via bit-group compare
            bpg = BITS_PER_SYM
            bh = bits_hat.reshape(-1, bpg); bt = bits_true.reshape(-1, bpg)
            sym_errs += (bh != bt).any(dim=-1).sum().item()
            tot_syms += bt.shape[0]
        ber[i] = bit_errs / tot_bits
        ser[i] = sym_errs / tot_syms
    return ber, ser

results = {}
for name, mod, dem, is_otfs in [
    ("OFDM", ofdm_mod, ofdm_dem, False),
    ("OCDM", ocdm_mod, ocdm_dem, False),
    ("AFDM", afdm_mod, afdm_dem, False),
    ("OTFS", otfs_mod, otfs_dem, True),
]:
    ber, ser = run_ber(mod, dem, is_otfs)
    results[name] = (ber, ser)
    print(f"{name}: BER={np.array2string(ber, precision=4)}")

# ---------------- Plots ----------------
colors = {"OFDM": "black", "OCDM": "tab:orange", "AFDM": "tab:green", "OTFS": "tab:blue"}

fig, axes = plt.subplots(1, 2, figsize=(13, 5))
ax = axes[0]
names = list(papr_vals.keys())
means = [papr_vals[k].mean() for k in names]
p99 = [np.percentile(papr_vals[k], 99) for k in names]
x = np.arange(len(names))
ax.bar(x - 0.15, means, width=0.3, label="Mean PAPR", color=[colors[k] for k in names])
ax.bar(x + 0.15, p99, width=0.3, label="99th pct PAPR", color=[colors[k] for k in names], alpha=0.5)
ax.set_xticks(x); ax.set_xticklabels(names)
ax.set_ylabel("PAPR [dB]")
ax.set_title("PAPR comparison (LiSAT waveforms.py, N=64, 16-QAM)")
ax.legend(); ax.grid(True, alpha=0.3)

ax2 = axes[1]
for name in ["OFDM", "OCDM", "AFDM", "OTFS"]:
    ber, ser = results[name]
    ax2.semilogy(snr_db_range, ber, marker='o', label=name, color=colors[name])
ax2.set_xlabel("SNR [dB]"); ax2.set_ylabel("BER")
ax2.set_title("AWGN BER vs SNR (LiSAT waveforms.py, 16-QAM)")
ax2.legend(); ax2.grid(True, which='both', alpha=0.3)

plt.tight_layout()
plt.savefig("D:/doubly-dispersive-channel-simulation-main/lisat_check/lisat_waveform_check.png", dpi=150)
print("\n=== DONE, saved lisat_waveform_check.png ===")
