#HEM Model, adapted from Project Sunride, United Kingdom

import numpy as np
import matplotlib.pyplot as plt
import CoolProp.CoolProp as CP
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from mpl_toolkits.mplot3d import Axes3D


def A(d):
    return np.pi * (d / 2) ** 2

def mHEM(A, rho2, h1, h2, N, Cd=0.66):
        return N * A * Cd * rho2 * np.sqrt(2*(h1 - h2))


T1 = 15.5 # Upstream temperature (C)
P2 = 31e5 # Downstream pressure (Pa)
N = 1 # Number of orrifaces

import CoolProp.CoolProp as CP


def HEM_CP(T1, P2, subst="NitrousOxide"):
    state = CP.AbstractState("HEOS", subst)

    # 1. Upstream saturated liquid (Q = 0)
    state.update(CP.QT_INPUTS, 0.0, T1 + 273.15)
    h1 = state.hmass()
    s1 = state.smass()

    # 2. Downstream isentropic state (P = P2, s = s1)
    # Remove the next line if you want standard two-phase HEM flashing
    state.specify_phase(CP.iphase_supercritical_liquid)

    state.update(CP.PSmass_INPUTS, P2, s1)
    assert state.phase() == CP.iphase_supercritical_liquid

    h2 = state.hmass()
    rho2 = state.rhomass()

    return h1, h2, rho2


def plotting(d, T, N, Cd=0.66):
    h1, h2, rho2 = HEM_CP(T, P2)

    m_HEM_vectorized = np.vectorize(mHEM)
    return m_HEM_vectorized(A(d), rho2, h1, h2, N, Cd)

def HEMmassflowrate():
    temps = np.linspace(-20, 32, 500)
    d = np.linspace(0.1, 4, 100) / 1000


    plt.figure()
    # set colors to viridis colormap
    colormap = plt.cm.viridis
    color = iter(colormap(np.linspace(0, 1, len(temps))))

    for T in temps:
        m = plotting(d, T, N)
        plt.plot(d * 1000, m, label=f'{T:.1f} C', color=next(color), linewidth=.3, alpha=1)

    sm = ScalarMappable(cmap=colormap, norm=Normalize(vmin=min(temps), vmax=max(temps)))
    sm.set_array([])
    ax = plt.gca()
    cbar = plt.colorbar(sm, ax=ax)
    cbar.set_label('Tank Temperature (C)')
    plt.xlabel('Single Orifice Diameter (mm)')
    plt.ylabel('Mass Flow Rate (kg/s)')
    plt.title('HEM Model, 1 Orifices')
    plt.show()


def NormalisedHEM():
    # Plot Mass Flow Rate vs Upstream temperature for d =0.5mm and N = 1
    temps = np.linspace(-20, 70, 500)
    d = 1
    N = 1

    m = plotting(d, temps, N) / 0.66
    plt.figure()
    plt.plot(temps, m/1e4)
    plt.xlabel('Tank Temperature (C)')
    plt.ylabel(r'Normalised mMss Flow Rate $\left( \frac{kg}{m^2} \right)$ x10$^{-4}$')
    plt.title('Mass Flow Rate vs Tank Temperature')
    plt.tight_layout()
    plt.show()

NormalisedHEM()
HEMmassflowrate()
# # Plot Mass Flow Rate vs diameter(0 to 2) for ambient temperature of 20degC and N = 12
# temps = 20
# d = np.linspace(0.1, 2.5, 100) / 1000
# N = 12
#
# m = plotting(d, temps, N)
# plt.figure()
# plt.plot(d * 1000, m)
# plt.xlabel('Diameter (mm)')
# plt.ylabel('Mass Flow Rate (kg/s)')
# plt.title('Mass Flow Rate vs Diameter')
# plt.tight_layout()
# plt.show()