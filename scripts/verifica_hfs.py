#!/usr/bin/env python3
"""Verifica numericamente Delta E_hfs, A_10 e C a partir de constantes CODATA."""
import math

# CODATA 2018 / valores usados no artigo
alpha = 7.2973525693e-3
m_e = 9.1093837015e-31
m_p = 1.67262192369e-27
c = 299792458.0
h = 6.62607015e-34
k_B = 1.380649e-23
e = 1.602176634e-19
mu_0 = 1.25663706212e-6
g_p = 5.5856946893
a_e = 1.15965218128e-3
g_e_dirac = 2.0
g_e = 2.0 * (1.0 + a_e)
m_r = m_e * m_p / (m_e + m_p)

def delta_E_J(g_e_val, use_mr=True):
    mass = m_r if use_mr else m_e
    # (2/3) g_e g_p alpha^4 (m_r^3 /(m_e m_p)) c^2
    return (2.0/3.0) * g_e_val * g_p * alpha**4 * (mass**3) / (m_e * m_p) * c**2

nu_exp = 1420.405751768e6  # Hz, Hellwig et al.

print("=== desdobramento hiperfino ===")
for label, ge, mr in [
    ("Dirac puro (g_e=2, m_e)", g_e_dirac, False),
    ("+ massa reduzida", g_e_dirac, True),
    ("+ momento anomalo", g_e, True),
]:
    dE = delta_E_J(ge, mr)
    nu = dE / h
    ppm = (nu - nu_exp) / nu_exp * 1e6
    print(f"{label:28s}  {dE/e*1e6:8.5f} ueV   {nu/1e6:10.3f} MHz   {ppm:+7.0f} ppm")
print(f"{'experimental':28s}  {h*nu_exp/e*1e6:8.5f} ueV   {nu_exp/1e6:10.3f} MHz")

# A_10 (dipolo magnetico, ordem de grandeza com |mu|~mu_B)
mu_B = e * h / (4 * math.pi * m_e)
# formula do artigo: 16 pi^3 nu0^3 mu0 / (3 h c^3) |mu|^2
# valor padrao citado: 2.869e-15 s^-1
A10 = 2.869e-15
C = 32 * math.pi * k_B * nu_exp**2 / (3 * h * c**3 * A10)
# converter para cm^-2 (K km/s)^-1:  * 1e-4 (m->cm)^2 wait
# N_HI [cm^-2] = C_SI * int T dv_SI
# dv in km/s = dv_SI / 1e3; T in K
# C_paper = C_SI * 1e3  in m^-2 (K km/s)^-1, then *1e-4 for cm^-2
C_m2 = C * 1e3
C_cm2 = C_m2 * 1e-4
print("\n=== radio ===")
print(f"A_10 = {A10:.3e} s^-1")
print(f"C    = {C_cm2:.4e} cm^-2 (K km/s)^-1")
