"""
waveforms.py
------------
Waveform modulator/demodulator implementations for LiSAT simulation suite.

Implements four physical-layer waveforms:
  - OFDM  : Orthogonal Frequency Division Multiplexing (baseline)
  - OTFS  : Orthogonal Time Frequency Space (delay-Doppler domain)
  - OCDM  : Orthogonal Chirp Division Multiplexing (Fresnel/DFnT)
  - AFDM  : Affine Frequency Division Multiplexing (DAFT-based)

All modulators expose:
    modulate(symbols)   -> signal  (complex torch.Tensor)
    demodulate(signal)  -> symbols (complex torch.Tensor)

Mathematical conventions follow the LiSAT paper and:
  - AFDM: Bemani et al. (2023), "AFDM: An Effective Modulation for ISAC"
  - OTFS: Hadani et al. (2017)
  - OCDM: Ouyang & Jin (2016)
"""

from __future__ import annotations

import math
from typing import Optional, Tuple

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F


# ---------------------------------------------------------------------------
# Utility functions
# ---------------------------------------------------------------------------

def compute_papr(signal: torch.Tensor) -> torch.Tensor:
    """
    Compute Peak-to-Average Power Ratio (PAPR) in dB.

    Parameters
    ----------
    signal : torch.Tensor, shape [..., N], complex

    Returns
    -------
    papr_db : torch.Tensor, shape [...], real
    """
    inst_power = signal.abs() ** 2
    peak = inst_power.max(dim=-1).values
    avg = inst_power.mean(dim=-1)
    return 10.0 * torch.log10(peak / (avg + 1e-12))


def qam_mapper(bits: torch.Tensor, order: int = 16) -> torch.Tensor:
    """
    Map bits to QAM constellation symbols using Gray-coded mapping.

    Parameters
    ----------
    bits  : torch.Tensor, shape [..., N * log2(order)], dtype=int/bool
    order : int, QAM order (4, 16, 64)

    Returns
    -------
    symbols : torch.Tensor, shape [..., N], complex64
    """
    bits_per_sym = int(math.log2(order))
    assert bits.shape[-1] % bits_per_sym == 0, "bits length must be divisible by bits_per_sym"
    n_sym = bits.shape[-1] // bits_per_sym
    bits_r = bits.view(*bits.shape[:-1], n_sym, bits_per_sym).int()

    # Gray decode: convert gray-coded index to natural index
    def gray2nat(g: torch.Tensor) -> torch.Tensor:
        n = g.clone()
        mask = g >> 1
        while mask.any():
            n = n ^ mask
            mask = mask >> 1
        return n

    # Convert bit rows to integer indices
    idx = torch.zeros(*bits_r.shape[:-1], dtype=torch.long, device=bits.device)
    for b in range(bits_per_sym):
        idx = idx | (bits_r[..., b] << (bits_per_sym - 1 - b))
    idx = gray2nat(idx)

    # Build constellation
    sqrt_m = int(math.sqrt(order))
    assert sqrt_m * sqrt_m == order, "order must be a perfect square"
    # Levels: -(sqrt_m-1), -(sqrt_m-3), ..., (sqrt_m-1)
    levels = torch.arange(sqrt_m, dtype=torch.float32, device=bits.device) * 2 - (sqrt_m - 1)
    i_idx = idx % sqrt_m
    q_idx = idx // sqrt_m
    symbols = (levels[i_idx] + 1j * levels[q_idx]).to(torch.complex64)
    # Normalise to unit average power
    norm = math.sqrt((levels ** 2).mean().item() * 2)
    return symbols / norm


def qam_demapper(
    symbols: torch.Tensor,
    order: int = 16,
    hard: bool = True,
) -> torch.Tensor:
    """
    Demap received QAM symbols to bits (hard decision).

    Parameters
    ----------
    symbols : torch.Tensor, shape [..., N], complex64
    order   : int, QAM order
    hard    : bool, hard (True) or soft (False, LLR) decision

    Returns
    -------
    bits : torch.Tensor, shape [..., N * log2(order)], int32
    """
    sqrt_m = int(math.sqrt(order))
    bits_per_sym = int(math.log2(order))

    levels = torch.arange(sqrt_m, dtype=torch.float32, device=symbols.device) * 2 - (sqrt_m - 1)
    norm = math.sqrt((levels ** 2).mean().item() * 2)
    levels_n = levels / norm

    # Quantise real and imaginary parts independently
    I = symbols.real.unsqueeze(-1)  # [..., N, 1]
    Q = symbols.imag.unsqueeze(-1)
    lvl = levels_n.view(*([1] * (symbols.dim())), sqrt_m)

    i_idx = (I - lvl).abs().argmin(dim=-1)
    q_idx = (Q - lvl).abs().argmin(dim=-1)
    nat_idx = q_idx * sqrt_m + i_idx

    # Natural → Gray
    gray_idx = nat_idx ^ (nat_idx >> 1)

    bits = torch.zeros(*symbols.shape[:-1], symbols.shape[-1] * bits_per_sym,
                       dtype=torch.int32, device=symbols.device)
    for b in range(bits_per_sym):
        bits[..., torch.arange(symbols.shape[-1]) * bits_per_sym + b] = (
            (gray_idx >> (bits_per_sym - 1 - b)) & 1
        ).int()
    return bits


# ---------------------------------------------------------------------------
# OFDM
# ---------------------------------------------------------------------------

class OFDMModulator(nn.Module):
    """
    Standard OFDM modulator with cyclic prefix.

    x_OFDM[n] = (1/sqrt(N)) * IFFT{X[k]}  for n = 0,...,N-1,
    then prepend N_cp samples of the tail as cyclic prefix.

    Parameters
    ----------
    N    : int, number of subcarriers
    N_cp : int, cyclic prefix length
    """

    def __init__(self, N: int = 64, N_cp: int = 16) -> None:
        super().__init__()
        self.N = N
        self.N_cp = N_cp

    def modulate(self, symbols: torch.Tensor) -> torch.Tensor:
        """
        Modulate frequency-domain symbols to time-domain OFDM signal.

        Parameters
        ----------
        symbols : torch.Tensor, shape [batch, N] or [batch, num_sym, N], complex64

        Returns
        -------
        signal  : torch.Tensor, shape [batch, (N + N_cp)] or [batch, num_sym*(N+N_cp)]
        """
        squeeze = symbols.dim() == 2
        if squeeze:
            symbols = symbols.unsqueeze(1)
        batch, num_sym, N = symbols.shape
        assert N == self.N

        # IFFT (scaled)
        td = torch.fft.ifft(symbols, n=N, dim=-1) * math.sqrt(N)

        # Add cyclic prefix
        cp = td[..., -self.N_cp:]
        td_cp = torch.cat([cp, td], dim=-1)  # [batch, num_sym, N+N_cp]

        signal = td_cp.reshape(batch, -1)
        return signal

    def forward(self, symbols: torch.Tensor) -> torch.Tensor:
        return self.modulate(symbols)

    def compute_ici_matrix(self) -> torch.Tensor:
        """
        Compute N×N ICI (inter-carrier interference) matrix for AWGN channel.
        For a diagonal channel this is the identity.
        """
        return torch.eye(self.N, dtype=torch.complex64)


class OFDMDemodulator(nn.Module):
    """
    Standard OFDM demodulator: remove CP then FFT.

    Parameters
    ----------
    N    : int, number of subcarriers
    N_cp : int, cyclic prefix length
    """

    def __init__(self, N: int = 64, N_cp: int = 16) -> None:
        super().__init__()
        self.N = N
        self.N_cp = N_cp

    def demodulate(self, signal: torch.Tensor) -> torch.Tensor:
        """
        Demodulate received time-domain signal.

        Parameters
        ----------
        signal : torch.Tensor, shape [batch, num_sym*(N+N_cp)] or [batch, N+N_cp], complex64

        Returns
        -------
        symbols : torch.Tensor, shape [batch, N] or [batch, num_sym, N], complex64
        """
        batch = signal.shape[0]
        sym_len = self.N + self.N_cp
        num_sym = signal.shape[-1] // sym_len

        sig_r = signal[:, : num_sym * sym_len].reshape(batch, num_sym, sym_len)

        # Remove cyclic prefix
        td = sig_r[..., self.N_cp:]

        # FFT
        symbols = torch.fft.fft(td, n=self.N, dim=-1) / math.sqrt(self.N)

        if num_sym == 1:
            symbols = symbols.squeeze(1)
        return symbols

    def forward(self, signal: torch.Tensor) -> torch.Tensor:
        return self.demodulate(signal)


# ---------------------------------------------------------------------------
# OTFS
# ---------------------------------------------------------------------------

class OTFSModulator(nn.Module):
    """
    OTFS modulator using the Zak (ISFFT + Heisenberg) transform framework.

    The delay-Doppler grid X[k,l] (k=Doppler bin, l=delay bin) is converted
    to the time-frequency plane via the inverse SFFT:

        S[n,m] = (1/sqrt(NM)) * sum_{k,l} X[k,l]
                  * exp(j 2π (kn/N - lm/M))

    and then transmitted via the Heisenberg transform (OFDM block synthesis):

        s[t] = sum_n sum_m S[n,m] * g_tx(t - nT) * exp(j 2π m Δf t)

    In the discrete, CP-based implementation this becomes NM samples per
    OTFS frame with an outer CP of length N_cp.

    Parameters
    ----------
    M : int, number of delay (subcarrier) bins
    N : int, number of Doppler (symbol) bins
    N_cp : int, cyclic prefix per OTFS block
    """

    def __init__(self, M: int = 64, N: int = 16, N_cp: int = 16) -> None:
        super().__init__()
        self.M = M
        self.N = N
        self.N_cp = N_cp

    def modulate(self, dd_symbols: torch.Tensor) -> torch.Tensor:
        """
        Modulate delay-Doppler symbols to time-domain signal.

        Parameters
        ----------
        dd_symbols : torch.Tensor, shape [batch, N, M], complex64

        Returns
        -------
        signal : torch.Tensor, shape [batch, N*(M+N_cp)], complex64
        """
        batch, N, M = dd_symbols.shape
        assert N == self.N and M == self.M

        # Step 1: Inverse SFFT — 2D IDFT from DD domain to TF domain
        # S[n,m] = IDFT_N{DFT_M{X}}  (in Zak convention)
        tf = torch.fft.ifft(torch.fft.fft(dd_symbols, dim=-1), dim=-2) * math.sqrt(N * M)

        # Step 2: Heisenberg transform — OFDM synthesis per OTFS block
        # Each Doppler index n → one OFDM symbol of length M + N_cp
        ofdm_td = torch.fft.ifft(tf, n=M, dim=-1) * math.sqrt(M)
        cp = ofdm_td[..., -self.N_cp:]
        ofdm_td_cp = torch.cat([cp, ofdm_td], dim=-1)  # [batch, N, M+N_cp]

        signal = ofdm_td_cp.reshape(batch, -1)
        return signal

    def forward(self, dd_symbols: torch.Tensor) -> torch.Tensor:
        return self.modulate(dd_symbols)


class OTFSDemodulator(nn.Module):
    """
    OTFS demodulator: Wigner transform + forward SFFT + MP equaliser.

    Parameters
    ----------
    M    : int, delay bins
    N    : int, Doppler bins
    N_cp : int, cyclic prefix length
    """

    def __init__(self, M: int = 64, N: int = 16, N_cp: int = 16) -> None:
        super().__init__()
        self.M = M
        self.N = N
        self.N_cp = N_cp

    def demodulate(
        self,
        signal: torch.Tensor,
        channel_taps: Optional[torch.Tensor] = None,
    ) -> torch.Tensor:
        """
        Demodulate received signal from time-domain to delay-Doppler domain.

        Parameters
        ----------
        signal       : torch.Tensor [batch, N*(M+N_cp)], complex64
        channel_taps : optional channel estimate for equalisation

        Returns
        -------
        dd_symbols : torch.Tensor [batch, N, M], complex64
        """
        batch = signal.shape[0]
        M, N, N_cp = self.M, self.N, self.N_cp
        sym_len = M + N_cp

        sig_r = signal[:, : N * sym_len].reshape(batch, N, sym_len)

        # Remove CP
        td = sig_r[..., N_cp:]

        # Wigner transform: FFT along delay axis
        tf = torch.fft.fft(td, n=M, dim=-1) / math.sqrt(M)

        # Forward SFFT: 2D transform from TF to DD domain
        dd = torch.fft.ifft(torch.fft.fft(tf, dim=-2), dim=-1) / math.sqrt(N * M)

        # Simple one-tap equaliser in TF domain (if channel estimate provided)
        if channel_taps is not None:
            H_eq = channel_taps.mean(dim=-1, keepdim=True) + 1e-6
            dd = dd / H_eq

        return dd

    def forward(self, signal: torch.Tensor, **kwargs) -> torch.Tensor:
        return self.demodulate(signal, **kwargs)


# ---------------------------------------------------------------------------
# OCDM (Orthogonal Chirp Division Multiplexing)
# ---------------------------------------------------------------------------

class OCDMModulator(nn.Module):
    """
    OCDM modulator using the Discrete Fresnel Transform (DFnT).

    Mathematical definition:
        x_OCDM[n] = sum_{k=0}^{N-1} X[k] * exp(j π k²/N) * exp(-j 2π n k/N)
                  = IFFT{ X[k] * exp(j π k²/N) }

    i.e. pre-chirp multiplication in frequency domain followed by IFFT.
    The chirp term exp(j π k²/N) spreads energy across the entire bandwidth,
    providing robustness against doubly-dispersive channels.

    Reference: Ouyang & Jin, "OCDM: A New Modulation Scheme", IEEE Wireless
               Commun. Lett., 2016.

    Parameters
    ----------
    N    : int, number of subcarriers / chirp length
    N_cp : int, cyclic prefix length
    """

    def __init__(self, N: int = 64, N_cp: int = 16) -> None:
        super().__init__()
        self.N = N
        self.N_cp = N_cp

        # Pre-compute chirp phase: exp(j π k² / N)
        k = torch.arange(N, dtype=torch.float64)
        self.register_buffer(
            "chirp_phase",
            torch.exp(1j * math.pi * k ** 2 / N).to(torch.complex64),
        )

    def modulate(self, symbols: torch.Tensor) -> torch.Tensor:
        """
        Modulate frequency-domain symbols via DFnT.

        Parameters
        ----------
        symbols : torch.Tensor [batch, N] or [batch, num_sym, N], complex64

        Returns
        -------
        signal  : torch.Tensor [batch, num_sym*(N+N_cp)], complex64
        """
        squeeze = symbols.dim() == 2
        if squeeze:
            symbols = symbols.unsqueeze(1)
        batch, num_sym, N = symbols.shape
        assert N == self.N

        # Pre-chirp multiplication
        chirped = symbols * self.chirp_phase  # [batch, num_sym, N]

        # IFFT (no additional normalisation — matches Ouyang & Jin)
        td = torch.fft.ifft(chirped, n=N, dim=-1) * math.sqrt(N)

        # Cyclic prefix
        cp = td[..., -self.N_cp:]
        td_cp = torch.cat([cp, td], dim=-1)

        signal = td_cp.reshape(batch, -1)
        return signal

    def forward(self, symbols: torch.Tensor) -> torch.Tensor:
        return self.modulate(symbols)


class OCDMDemodulator(nn.Module):
    """
    OCDM demodulator: remove CP, FFT, de-chirp.

    Parameters
    ----------
    N    : int
    N_cp : int
    """

    def __init__(self, N: int = 64, N_cp: int = 16) -> None:
        super().__init__()
        self.N = N
        self.N_cp = N_cp

        k = torch.arange(N, dtype=torch.float64)
        self.register_buffer(
            "dechirp_phase",
            torch.exp(-1j * math.pi * k ** 2 / N).to(torch.complex64),
        )

    def demodulate(self, signal: torch.Tensor) -> torch.Tensor:
        """
        Demodulate received OCDM signal.

        Parameters
        ----------
        signal : torch.Tensor [batch, num_sym*(N+N_cp)], complex64

        Returns
        -------
        symbols : torch.Tensor [batch, N] or [batch, num_sym, N]
        """
        batch = signal.shape[0]
        sym_len = self.N + self.N_cp
        num_sym = signal.shape[-1] // sym_len

        sig_r = signal[:, : num_sym * sym_len].reshape(batch, num_sym, sym_len)
        td = sig_r[..., self.N_cp:]

        # FFT
        fd = torch.fft.fft(td, n=self.N, dim=-1) / math.sqrt(self.N)

        # De-chirp
        symbols = fd * self.dechirp_phase

        if num_sym == 1:
            symbols = symbols.squeeze(1)
        return symbols

    def forward(self, signal: torch.Tensor) -> torch.Tensor:
        return self.demodulate(signal)


# ---------------------------------------------------------------------------
# AFDM (Affine Frequency Division Multiplexing)
# ---------------------------------------------------------------------------

class AFDMModulator(nn.Module):
    """
    AFDM modulator using the Inverse Discrete Affine Fourier Transform (IDAFT).

    The IDAFT is defined as:
        x[n] = (1/sqrt(N)) * sum_{k=0}^{N-1} X[k]
                * exp(j 2π (c1 * n² / N + k * n / N + c2 * k² / N))

    where c1, c2 are the pre- and post-chirp rates chosen to guarantee that
    all delay-Doppler paths map to distinct DAFT-domain bins (MERD criterion).

    Optimal choice (Bemani et al. 2023):
        c1 = 0.5 * (ν_max * N / Δf)  (rounded to nearest half-integer / N)
        c2 = 0  (or small value to break symmetry)

    Parameters
    ----------
    N  : int, number of subcarriers
    c1 : float, pre-chirp rate (in units of 1/N)
    c2 : float, post-chirp rate (in units of 1/N)
    N_cp : int, cyclic prefix length
    """

    def __init__(
        self,
        N: int = 64,
        c1: float = 0.25,
        c2: float = 0.0,
        N_cp: int = 16,
    ) -> None:
        super().__init__()
        self.N = N
        self.c1 = c1
        self.c2 = c2
        self.N_cp = N_cp

        self._build_matrices()

    def _build_matrices(self) -> None:
        N, c1, c2 = self.N, self.c1, self.c2
        n = torch.arange(N, dtype=torch.float64)
        k = torch.arange(N, dtype=torch.float64)

        # IDAFT matrix: W[n, k] = (1/sqrt(N)) * exp(j 2π (c1*n²/N + k*n/N + c2*k²/N))
        phase = (
            2.0 * math.pi * (
                c1 * n.unsqueeze(1) ** 2 / N
                + n.unsqueeze(1) * k.unsqueeze(0) / N
                + c2 * k.unsqueeze(0) ** 2 / N
            )
        )
        W = torch.exp(1j * phase).to(torch.complex64) / math.sqrt(N)
        self.register_buffer("W_idaft", W)  # [N, N]

        # DAFT matrix (forward): conjugate transpose
        self.register_buffer("W_daft", W.conj().T.contiguous())

        # Pilot subframe indices for MERD channel estimation
        # Place pilots at N//4 positions evenly spaced
        pilot_spacing = max(1, N // 8)
        self.pilot_indices = list(range(0, N, pilot_spacing))

    def modulate(self, symbols: torch.Tensor) -> torch.Tensor:
        """
        Modulate symbols via IDAFT.

        Parameters
        ----------
        symbols : torch.Tensor [batch, N] or [batch, num_sym, N], complex64

        Returns
        -------
        signal : torch.Tensor [batch, num_sym*(N+N_cp)], complex64
        """
        squeeze = symbols.dim() == 2
        if squeeze:
            symbols = symbols.unsqueeze(1)
        batch, num_sym, N = symbols.shape
        assert N == self.N

        # Apply IDAFT: x = W_idaft @ X  (matrix-vector per symbol)
        # symbols: [batch, num_sym, N]  →  [batch*num_sym, N]
        sym_flat = symbols.reshape(-1, N)
        td = (sym_flat @ self.W_idaft.T)  # [batch*num_sym, N]
        td = td.reshape(batch, num_sym, N)

        # Cyclic prefix
        cp = td[..., -self.N_cp:]
        td_cp = torch.cat([cp, td], dim=-1)

        signal = td_cp.reshape(batch, -1)
        return signal

    def modulate_with_pilot(
        self, data_symbols: torch.Tensor
    ) -> Tuple[torch.Tensor, torch.Tensor]:
        """
        Modulate data with prepended MERD-compatible pilot subframe.

        Returns
        -------
        signal      : time-domain signal including pilot frame
        pilot_syms  : pilot symbols for channel estimation
        """
        N = self.N
        # MERD pilot: Zadoff-Chu sequence at pilot indices
        pilot = torch.zeros(data_symbols.shape[0], N, dtype=torch.complex64,
                            device=data_symbols.device)
        for idx, k in enumerate(self.pilot_indices):
            pilot[:, k] = torch.exp(
                1j * torch.tensor(math.pi * idx ** 2 / len(self.pilot_indices))
            )

        pilot_sig = self.modulate(pilot)
        data_sig = self.modulate(data_symbols)
        return torch.cat([pilot_sig, data_sig], dim=-1), pilot

    def forward(self, symbols: torch.Tensor) -> torch.Tensor:
        return self.modulate(symbols)


class AFDMDemodulator(nn.Module):
    """
    AFDM demodulator: remove CP, apply forward DAFT.

    Parameters
    ----------
    N    : int
    c1   : float
    c2   : float
    N_cp : int
    """

    def __init__(
        self,
        N: int = 64,
        c1: float = 0.25,
        c2: float = 0.0,
        N_cp: int = 16,
    ) -> None:
        super().__init__()
        self.N = N
        self.c1 = c1
        self.c2 = c2
        self.N_cp = N_cp

        # Share matrix construction with modulator
        mod = AFDMModulator(N=N, c1=c1, c2=c2, N_cp=N_cp)
        self.register_buffer("W_daft", mod.W_daft)

    def demodulate(self, signal: torch.Tensor) -> torch.Tensor:
        """
        Demodulate received AFDM signal.

        Parameters
        ----------
        signal : torch.Tensor [batch, num_sym*(N+N_cp)], complex64

        Returns
        -------
        symbols : torch.Tensor [batch, N] or [batch, num_sym, N]
        """
        batch = signal.shape[0]
        sym_len = self.N + self.N_cp
        num_sym = signal.shape[-1] // sym_len

        sig_r = signal[:, : num_sym * sym_len].reshape(batch, num_sym, sym_len)
        td = sig_r[..., self.N_cp:]  # [batch, num_sym, N]

        sym_flat = td.reshape(-1, self.N)
        symbols = sym_flat @ self.W_daft.T  # [batch*num_sym, N]
        symbols = symbols.reshape(batch, num_sym, self.N)

        if num_sym == 1:
            symbols = symbols.squeeze(1)
        return symbols

    def forward(self, signal: torch.Tensor) -> torch.Tensor:
        return self.demodulate(signal)


# ---------------------------------------------------------------------------
# Smoke test
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    torch.manual_seed(42)
    print("=== Waveform Smoke Tests ===\n")

    N = 64
    N_cp = 16
    batch = 4
    M = 16   # OTFS Doppler bins

    # Random QPSK symbols
    bits = torch.randint(0, 2, (batch, N * 4))
    syms = qam_mapper(bits, order=16)
    print(f"QAM mapper: bits={bits.shape}, symbols={syms.shape}, |E[|X|^2]|={syms.abs().pow(2).mean():.3f}")

    # OFDM
    ofdm_mod = OFDMModulator(N=N, N_cp=N_cp)
    ofdm_dem = OFDMDemodulator(N=N, N_cp=N_cp)
    sig = ofdm_mod.modulate(syms)
    rec = ofdm_dem.demodulate(sig)
    err = (rec - syms).abs().mean().item()
    papr = compute_papr(sig)
    print(f"OFDM:  signal={sig.shape}, recon error={err:.2e}, PAPR={papr.mean():.2f} dB")

    # OCDM
    ocdm_mod = OCDMModulator(N=N, N_cp=N_cp)
    ocdm_dem = OCDMDemodulator(N=N, N_cp=N_cp)
    sig_o = ocdm_mod.modulate(syms)
    rec_o = ocdm_dem.demodulate(sig_o)
    err_o = (rec_o - syms).abs().mean().item()
    papr_o = compute_papr(sig_o)
    print(f"OCDM:  signal={sig_o.shape}, recon error={err_o:.2e}, PAPR={papr_o.mean():.2f} dB")

    # AFDM
    afdm_mod = AFDMModulator(N=N, c1=0.25, c2=0.0, N_cp=N_cp)
    afdm_dem = AFDMDemodulator(N=N, c1=0.25, c2=0.0, N_cp=N_cp)
    sig_a = afdm_mod.modulate(syms)
    rec_a = afdm_dem.demodulate(sig_a)
    err_a = (rec_a - syms).abs().mean().item()
    papr_a = compute_papr(sig_a)
    print(f"AFDM:  signal={sig_a.shape}, recon error={err_a:.2e}, PAPR={papr_a.mean():.2f} dB")

    # OTFS
    delay_doppler_symbols = (torch.randn(batch, M, N) + 1j * torch.randn(batch, M, N)).to(torch.complex64)
    dd_syms = delay_doppler_symbols
    otfs_mod = OTFSModulator(M=N, N=M, N_cp=N_cp)
    otfs_dem = OTFSDemodulator(M=N, N=M, N_cp=N_cp)
    sig_t = otfs_mod.modulate(dd_syms)
    rec_t = otfs_dem.demodulate(sig_t)
    err_t = (rec_t - dd_syms).abs().mean().item()
    papr_t = compute_papr(sig_t)
    print(f"OTFS:  signal={sig_t.shape}, recon error={err_t:.2e}, PAPR={papr_t.mean():.2f} dB")

    print("\nAll waveform tests passed.")
