#!/usr/bin/env python3
"""
atividade_nhi.py — Etapas 1-4 da atividade computacional da Secao 7b.

Da equacao de Dirac ao mapa de N_HI: aplica a constante
C = 32 pi kB nu0^2 / (3 h c^3 A_10) = 1.823e18 cm^-2 (K km/s)^-1
a um espectro de 21 cm.

Uso (a partir da raiz do repositório):
    python3 scripts/atividade_nhi.py --demo
        perfil SINTETICO (verifica o encadeamento de unidades).

    python3 scripts/atividade_nhi.py --fits espectro.fits
        espectro HI4PI real (FITS 1D com eixo espectral no header,
        ou um cubo do qual se extrai o pixel central).

Requer: numpy. O modo --fits requer astropy.
"""
import argparse
import sys
import numpy as np

# numpy >= 2.0 usa trapezoid; versoes anteriores, trapz
_trapz = getattr(np, "trapezoid", None) or np.trapz

# ----------------------------------------------------------------------
# Constantes derivadas no artigo (ver verifica_hfs.py para a deducao)
# ----------------------------------------------------------------------
NU0    = 1420.405751768e6      # Hz   — Eq. (nu0)
C_LUZ  = 299792458.0           # m/s
C_COL  = 1.8225e18             # cm^-2 (K km/s)^-1 — Eq. (constante-C)


def freq_para_velocidade(nu_hz):
    """Etapa 1: frequencia -> velocidade radial (km/s), CONVENCAO RADIO
    (v = c(nu0-nu)/nu0, a adotada pelo HI4PI; a convencao optica seria
    v = c(nu0-nu)/nu)."""
    return C_LUZ * (NU0 - nu_hz) / NU0 / 1e3      # o /1e3 e' o fator critico


def coluna_hi(v_kms, tb_k):
    """Etapa 2: N_HI = C * integral T_B dv, com v em km/s e T_B em K."""
    ordem = np.argsort(v_kms)
    integral = _trapz(tb_k[ordem], v_kms[ordem])   # K km/s
    return C_COL * integral, integral


def checagem_dimensional(integral, n_hi):
    """Etapa 2 (cont.): alerta de ordem de grandeza."""
    print(f"  integral T_B dv = {integral:10.2f} K km/s")
    print(f"  N_HI            = {n_hi:10.3e} cm^-2")
    if not (1e18 < n_hi < 1e23):
        print("  [ALERTA] N_HI fora da faixa tipica 1e19-1e22 cm^-2 do MIS.")
        print("           Verifique o fator 1e3 entre m/s e km/s.")
    else:
        print("  [OK] ordem de grandeza compativel com o meio interestelar.")


def massa_hi(n_hi_cm2, area_sr, dist_kpc):
    """Etapa 4: massa de HI em massas solares."""
    KPC_CM, M_H_G, MSOL_G = 3.0857e21, 1.6735e-24, 1.989e33
    area_cm2 = area_sr * (dist_kpc * KPC_CM) ** 2
    return n_hi_cm2 * area_cm2 * M_H_G / MSOL_G


def espectro_demo():
    """Perfil SINTETICO — nao sao dados reais. So testa o caminho de codigo."""
    v = np.linspace(-150.0, 150.0, 400)                 # km/s
    tb = (60.0 * np.exp(-0.5 * ((v - 5.0) / 7.0) ** 2)  # componente fria
          + 12.0 * np.exp(-0.5 * ((v + 30.0) / 25.0) ** 2))  # componente morna
    return v, tb


def espectro_fits(caminho):
    """Le um espectro HI4PI real a partir de um FITS."""
    try:
        from astropy.io import fits
        from astropy.wcs import WCS
    except ImportError:
        sys.exit("astropy nao instalado: pip install astropy")
    with fits.open(caminho) as hdul:
        dados, hdr = hdul[0].data, hdul[0].header
    wcs = WCS(hdr)
    if dados.ndim == 3:                       # cubo -> pixel central
        ny, nx = dados.shape[1], dados.shape[2]
        tb = dados[:, ny // 2, nx // 2]
        n = dados.shape[0]
        eixo = wcs.pixel_to_world_values(np.full(n, nx // 2),
                                         np.full(n, ny // 2),
                                         np.arange(n))[2]
    else:
        tb = np.asarray(dados).ravel()
        eixo = wcs.pixel_to_world_values(np.arange(tb.size))[0]
    eixo = np.asarray(eixo, dtype=float)
    # HI4PI distribui o eixo em m/s (VELO-LSR) ou em Hz (FREQ).
    ctype = hdr.get("CTYPE3", hdr.get("CTYPE1", "")).upper()
    if "FREQ" in ctype:
        v = freq_para_velocidade(eixo)
    else:
        v = eixo / 1e3 if np.nanmax(np.abs(eixo)) > 1e4 else eixo
    return v, np.asarray(tb, dtype=float)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--fits", help="espectro ou cubo HI4PI (FITS)")
    p.add_argument("--demo", action="store_true", help="perfil sintetico")
    p.add_argument("--area-sr", type=float, default=1e-6)
    p.add_argument("--dist-kpc", type=float, default=1.0)
    a = p.parse_args()

    if a.fits:
        v, tb = espectro_fits(a.fits)
        rotulo = f"HI4PI: {a.fits}"
    elif a.demo:
        v, tb = espectro_demo()
        rotulo = "PERFIL SINTETICO (nao e' um dado real)"
    else:
        p.error("escolha --demo ou --fits ARQUIVO")

    print(f"\n{rotulo}")
    print(f"  canais: {v.size}   faixa: {v.min():.1f} a {v.max():.1f} km/s")
    print(f"  T_B maximo: {np.nanmax(tb):.2f} K\n")

    n_hi, integral = coluna_hi(v, tb)
    checagem_dimensional(integral, n_hi)

    m = massa_hi(n_hi, a.area_sr, a.dist_kpc)
    print(f"\n  Etapa 4: para {a.area_sr:.1e} sr a {a.dist_kpc:.1f} kpc"
          f"  ->  M_HI = {m:.3e} M_sol")
    print("\n  Etapa 3 (validacao): compare com o mapa de N_HI do proprio")
    print("  levantamento na mesma linha de visada.\n")


if __name__ == "__main__":
    main()
