#!/usr/bin/env python3
"""
Extrai duas linhas de visada dos tiles HI4PI CAR e gera a figura dos
espectros do artigo:

  figuras/espectros_hi4pi.pdf   — área hachurada = ∫ T_B dv
                                  caixa com N_HI (comparação visual + numérica)
  figuras/hi4pi_numeros.tex     — macros da Tabela de N_HI

Uso (a partir da raiz do repositório):
    python3 scripts/extract_hi4pi_sightlines.py
    python3 scripts/extract_hi4pi_sightlines.py --data-dir /caminho/para/cubos

Requer: numpy, matplotlib. Cubos FITS: CAR_D14.fits e CAR_D09.fits
(HI4PI, CDS J/A+A/594/A116, CUBES/EQ2000). Pasta: ./data ou HI4PI_DATA_DIR.
"""
from __future__ import annotations

import argparse
import math
import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

_trapz = getattr(np, "trapezoid", None) or np.trapz
C_COL = 1.823e18  # cm^-2 (K km/s)^-1  — Eq. (constante-C) do artigo
REPO_ROOT = Path(__file__).resolve().parents[1]
CDS_CUBES = (
    "https://cdsarc.cds.unistra.fr/viz-bin/qcat?J/A+A/594/A116"
    " (CUBES/EQ2000)"
)

SIGHTLINES = [
    ("CAR_D14.fits", 276.42, -12.33, "(a) plano galáctico", "plano"),
    ("CAR_D09.fits", 170.00, -15.00, "(b) latitude intermediária", "halo"),
]


def read_header(path: Path) -> dict:
    hdr = {}
    with open(path, "rb") as fh:
        n = 0
        while True:
            raw = fh.read(80)
            if len(raw) < 80:
                break
            n += 80
            card = raw.decode("ascii", "replace")
            key = card[:8].strip()
            if key == "END":
                break
            if key and "=" in card[:10]:
                val = card[10:].split("/")[0].strip().strip("'")
                hdr[key] = val
        pad = (2880 - (n % 2880)) % 2880
        hdr["_offset"] = n + pad
    return hdr


def fnum(hdr, key) -> float:
    return float(hdr[key])


def inum(hdr, key) -> int:
    return int(float(hdr[key]))


def load_cube(path: Path):
    hdr = read_header(path)
    nx, ny, nv = inum(hdr, "NAXIS1"), inum(hdr, "NAXIS2"), inum(hdr, "NAXIS3")
    with open(path, "rb") as fh:
        fh.seek(hdr["_offset"])
        buf = fh.read(nv * ny * nx * 4)
    data = np.frombuffer(buf, dtype=">f4").reshape((nv, ny, nx)).astype(np.float64)
    data[~np.isfinite(data)] = 0.0
    data[np.abs(data) > 1e10] = 0.0
    return hdr, data


def world_to_pix(hdr, ra, dec):
    i = (ra - fnum(hdr, "CRVAL1")) / fnum(hdr, "CDELT1") + fnum(hdr, "CRPIX1") - 1.0
    j = (dec - fnum(hdr, "CRVAL2")) / fnum(hdr, "CDELT2") + fnum(hdr, "CRPIX2") - 1.0
    return int(round(i)), int(round(j))


def velocity_kms(hdr, nv):
    """Eixo espectral -> v [km/s] na convenção rádio do HI4PI.

    Cubos CAR entregam CTYPE3 = VRAD em m/s (convenção rádio).
    O /1e3 é o fator entre m/s e km/s. Se o eixo viesse em
    frequência, usar v = c (nu0 - nu)/nu0 (não a convenção óptica).
    """
    k = np.arange(nv, dtype=np.float64)
    v_ms = fnum(hdr, "CRVAL3") + (k + 1.0 - fnum(hdr, "CRPIX3")) * fnum(hdr, "CDELT3")
    ctype = str(hdr.get("CTYPE3", hdr.get("CTYPE1", ""))).upper()
    if "FREQ" in ctype:
        nu = v_ms  # ja em Hz se CDELT3 estiver em Hz
        return 299792458.0 * (1420.405751768e6 - nu) / 1420.405751768e6 / 1e3
    return v_ms / 1e3


def eq_to_gal(ra_deg, dec_deg):
    a_gp = math.radians(192.85948)
    d_gp = math.radians(27.12825)
    l_ncp = math.radians(122.93192)
    a = math.radians(ra_deg)
    d = math.radians(dec_deg)
    sin_b = math.sin(d) * math.sin(d_gp) + math.cos(d) * math.cos(d_gp) * math.cos(a - a_gp)
    b = math.asin(max(-1.0, min(1.0, sin_b)))
    y = math.cos(d) * math.sin(a - a_gp)
    x = math.sin(d) * math.cos(d_gp) - math.cos(d) * math.sin(d_gp) * math.cos(a - a_gp)
    l = (math.degrees(l_ncp - math.atan2(y, x))) % 360.0
    return l, math.degrees(b)


def coluna(v_kms, tb_k):
    """N_HI = C * ∫ T_B dv  (v em km/s, T_B em K)."""
    ordem = np.argsort(v_kms)
    integ = float(_trapz(tb_k[ordem], v_kms[ordem]))
    return C_COL * integ, integ


def extract(tile, ra, dec, data_dir: Path):
    path = data_dir / tile
    hdr, cube = load_cube(path)
    i, j = world_to_pix(hdr, ra, dec)
    nx, ny, nv = inum(hdr, "NAXIS1"), inum(hdr, "NAXIS2"), inum(hdr, "NAXIS3")
    i = min(max(i, 0), nx - 1)
    j = min(max(j, 0), ny - 1)
    tb = cube[:, j, i]
    v = velocity_kms(hdr, nv)
    n_hi, integ = coluna(v, tb)
    ra_pix = fnum(hdr, "CRVAL1") + (i + 1.0 - fnum(hdr, "CRPIX1")) * fnum(hdr, "CDELT1")
    dec_pix = fnum(hdr, "CRVAL2") + (j + 1.0 - fnum(hdr, "CRPIX2")) * fnum(hdr, "CDELT2")
    l, b = eq_to_gal(ra_pix, dec_pix)
    return {
        "tile": tile,
        "ra": ra_pix,
        "dec": dec_pix,
        "l": l,
        "b": b,
        "v": v,
        "tb": tb,
        "n_hi": n_hi,
        "integral": integ,
        "tb_max": float(np.nanmax(tb)),
    }


def sci_tex(x):
    coef, exp = f"{x:.2e}".split("e")
    return rf"{coef.replace('.', '{,}')}\times10^{{{int(exp)}}}"


def deg_tex(x, nd=2):
    return f"{x:.{nd}f}".replace(".", "{,}")


def fmt_sci_math(x):
    coef, exp = f"{x:.2e}".split("e")
    return rf"{coef}\times 10^{{{int(exp)}}}"


def save_figure(fig, pdf_path: Path):
    """PDF vetorial (Type 42). Requer apenas matplotlib."""
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42
    pdf_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(pdf_path, format="pdf", facecolor="white", edgecolor="none")
    print(f"  pdf {pdf_path}")


def plot_espectros(results, pdf_path: Path):
    """Dois painéis: área hachurada = ∫ T_B dv; caixa com N_HI."""
    plt.rcParams.update({
        "font.size": 9,
        "axes.labelsize": 10,
        "mathtext.fontset": "cm",
    })
    # ~15 cm de largura, cabe em textwidth=16 cm do artigo
    fig, axes = plt.subplots(2, 1, figsize=(5.7, 5.1), sharex=True)
    vmax = 180.0
    hatches = ["///", r"\\\\"]
    faces = ["#d9d9d9", "#cfcfcf"]

    for ax, r, hatch, face in zip(axes, results, hatches, faces):
        v, tb = r["v"], np.clip(r["tb"], 0.0, None)
        m = np.abs(v) <= vmax
        vv, tt = v[m], tb[m]
        ordem = np.argsort(vv)
        vv, tt = vv[ordem], tt[ordem]

        ax.fill_between(
            vv, 0.0, tt,
            hatch=hatch, facecolor=face, edgecolor="0.35",
            linewidth=0.4, alpha=1.0,
            label=r"área $=\int T_B\,\mathrm{d}v$",
        )
        ax.plot(vv, tt, color="black", lw=1.2)

        caixa = (
            r"$\int T_B\,\mathrm{d}v="
            + fmt_sci_math(r["integral"])
            + r"\,\mathrm{K\,km\,s^{-1}}$"
            + "\n"
            + r"$N_{\mathrm{HI}}="
            + fmt_sci_math(r["n_hi"])
            + r"\,\mathrm{cm^{-2}}$"
        )
        ax.text(
            0.99, 0.96, caixa,
            transform=ax.transAxes, ha="right", va="top",
            fontsize=8, linespacing=1.4,
            bbox=dict(boxstyle="round,pad=0.28", facecolor="white",
                      edgecolor="0.25", linewidth=0.7),
            clip_on=True, zorder=5,
        )
        ax.set_ylabel(r"$T_B$ [K]")
        ax.set_xlim(-vmax, vmax)
        ymax = max(1.28 * r["tb_max"], 1.0)
        ax.set_ylim(0.0, ymax)
        ax.text(0.02, 0.96, r["label"], transform=ax.transAxes,
                ha="left", va="top", fontsize=9, fontweight="bold")
        ax.legend(loc="upper left", bbox_to_anchor=(0.0, 0.80),
                  frameon=False, fontsize=7.5, handlelength=1.6)
        ax.axhline(0.0, color="0.5", lw=0.4)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    axes[-1].set_xlabel(r"$v$ [km s$^{-1}$]")
    razao = results[0]["n_hi"] / results[1]["n_hi"]
    fig.suptitle(
        r"$N_{\mathrm{HI}}=\mathcal{C}\int T_B\,\mathrm{d}v$  "
        rf"($\mathcal{{C}}=1.823\times 10^{{18}}$)  ·  "
        rf"razão (a)/(b)$={razao:.0f}$",
        fontsize=9, y=0.995,
    )
    fig.tight_layout(rect=(0.01, 0.01, 0.99, 0.96))
    save_figure(fig, pdf_path)
    plt.close(fig)


def resolve_data_dir(cli: Path | None) -> Path:
    if cli is not None:
        return cli.expanduser().resolve()
    env = os.environ.get("HI4PI_DATA_DIR")
    if env:
        return Path(env).expanduser().resolve()
    return (REPO_ROOT / "data").resolve()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=None,
        help="pasta com CAR_D14.fits e CAR_D09.fits (default: ./data ou HI4PI_DATA_DIR)",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=None,
        help="pasta de saída da figura e das macros (default: ./figuras)",
    )
    args = parser.parse_args()
    data_dir = resolve_data_dir(args.data_dir)
    out_dir = (args.out_dir or (REPO_ROOT / "figuras")).expanduser().resolve()

    missing = [tile for tile, *_ in SIGHTLINES if not (data_dir / tile).exists()]
    if missing:
        sys.exit(
            f"Tiles HI4PI ausentes em {data_dir}: {', '.join(missing)}\n"
            f"Baixe os cubos equatoriais CAR em {CDS_CUBES}\n"
            "e coloque-os nessa pasta, ou passe --data-dir / HI4PI_DATA_DIR."
        )

    out_dir.mkdir(parents=True, exist_ok=True)
    results = []
    for tile, ra, dec, label, _key in SIGHTLINES:
        print(f"lendo {tile} ...")
        r = extract(tile, ra, dec, data_dir)
        r["label"] = label
        results.append(r)
        print(
            f"  RA={r['ra']:.2f} Dec={r['dec']:.2f}"
            f"  l={r['l']:.1f} b={r['b']:.1f}"
            f"  int={r['integral']:.1f} K km/s  N_HI={r['n_hi']:.3e}"
            f"  Tmax={r['tb_max']:.2f} K"
        )

    pdf = out_dir / "espectros_hi4pi.pdf"
    plot_espectros(results, pdf)
    print("escrito", pdf)

    r0, r1 = results[0], results[1]
    macros = [
        r"% gerado por scripts/extract_hi4pi_sightlines.py --- nao editar a mao",
        r"\newcommand{\NhiPlano}{" + sci_tex(r0["n_hi"]) + "}",
        r"\newcommand{\IntPlano}{" + sci_tex(r0["integral"]) + "}",
        r"\newcommand{\TbPlano}{" + f"{r0['tb_max']:.0f}" + "}",
        r"\newcommand{\LonPlano}{" + deg_tex(r0["l"], 1) + "}",
        r"\newcommand{\LatPlano}{" + deg_tex(r0["b"], 1) + "}",
        r"\newcommand{\RaPlano}{" + deg_tex(r0["ra"], 2) + "}",
        r"\newcommand{\DecPlano}{" + deg_tex(r0["dec"], 2) + "}",
        r"\newcommand{\NhiHalo}{" + sci_tex(r1["n_hi"]) + "}",
        r"\newcommand{\IntHalo}{" + f"{r1['integral']:.0f}" + "}",
        r"\newcommand{\TbHalo}{" + deg_tex(r1["tb_max"], 1) + "}",
        r"\newcommand{\LonHalo}{" + deg_tex(r1["l"], 1) + "}",
        r"\newcommand{\LatHalo}{" + deg_tex(r1["b"], 1) + "}",
        r"\newcommand{\RaHalo}{" + deg_tex(r1["ra"], 2) + "}",
        r"\newcommand{\DecHalo}{" + deg_tex(r1["dec"], 2) + "}",
        r"\newcommand{\RazaoNhi}{" + f"{r0['n_hi']/r1['n_hi']:.0f}" + "}",
    ]
    (out_dir / "hi4pi_numeros.tex").write_text("\n".join(macros) + "\n")
    print("escrito", out_dir / "hi4pi_numeros.tex")


if __name__ == "__main__":
    main()
