import numpy as np
import matplotlib.pyplot as plt
import CoolProp

from fullflow import *

"""
Bluebird Aerospace


Two phase NHNE model injector dP sweep helper. Sweeps injector dPs and plots to mass flux to find choke point

The NHNE orifice:
    G_SPI  = Cd * sqrt(2 * rho_1 * dP)
    G_HEM  = Cd * rho_2 * sqrt(2 * (h_1 - h_2)),  isentropic flash to P2
    kappa  = sqrt(dP / (Psat - P2))
    G_NHNE = (kappa * G_SPI + G_HEM) / (1 + kappa)

"""
psia_to_pa = 6894.76
bar_to_pa = 1.0e5

design_chamber_pressure = 31.0 * bar_to_pa
chamber_pressure = State(design_chamber_pressure)

ox_temperature = 288.65

_sat_ref = CoolProp.AbstractState("HEOS", "NitrousOxide")
_sat_ref.update(CoolProp.QT_INPUTS, 0.0, ox_temperature)
ox_source_pressure = _sat_ref.p()

ox_injector_cda = 5.779e-5

# NHNE discharge coefficient for the ox injector. NOTE: if ox_injector_cda is
# already a Cd*A, leave this at 1.0. If it is a geometric area, use the
# cold-flow-validated Cd (e.g. 0.65) instead. Do not double-count.
ox_injector_nhne_cd = 1.0

ox_feed_cda = 1.0e-2


class NHNEOrifice(Component):
    """Two-phase injector orifice using the Dyer NHNE model.

    forward flow blends the single-phase incompressible and homogeneous
    equilibrium mass flow rates w/ Dyer non-equilibrium parameter kappa

    hopefully this works
    """

    def __init__(
        self,
        name,
        network,
        fluid,
        upstream_pressure,
        upstream_temperature,
        downstream_pressure,
        discharge_coefficient,
        cross_sectional_area,
        mass_flow=None,
        kappa=None,
        G_SPI=None,
        G_HEM=None,
        G_NHNE=None,
        upstream_density=None,
    ):
        backend_fluid = fluid
        self._liquid_state = CoolProp.AbstractState("HEOS", backend_fluid)
        self._liquid_state.specify_phase(CoolProp.iphase_liquid)
        self._flash_state = CoolProp.AbstractState("HEOS", backend_fluid)

        self._ref_T = None
        self._ref_psat = 0.0
        self._ref_rho = 0.0
        self._ref_h = 0.0
        self._ref_s = 0.0
        self._ref_drho_dp = 0.0
        self.setup()

    @property
    def ignored_export_attributes(self):
        return super().ignored_export_attributes | {
            "_liquid_state",
            "_flash_state",
            "_ref_T",
            "_ref_psat",
            "_ref_rho",
            "_ref_h",
            "_ref_s",
            "_ref_drho_dp",
        }

    def _refresh_reference(self, T1):
        """need saturated liquid reference values and the isothermal compressibility slope at init temp"""
        self._flash_state.update(CoolProp.QT_INPUTS, 0.0, T1)
        self._ref_psat = self._flash_state.p()
        self._ref_rho = self._flash_state.rhomass()
        self._ref_h = self._flash_state.hmass()
        self._ref_s = self._flash_state.smass()

        dp_fd = 1.0e5
        self._liquid_state.update(
            CoolProp.PT_INPUTS, self._ref_psat + dp_fd, T1
        )
        self._ref_drho_dp = (
            self._liquid_state.rhomass() - self._ref_rho
        ) / dp_fd
        self._ref_T = T1

    def _upstream_liquid(self, P1, T1):
        if self._ref_T != T1:
            self._refresh_reference(T1)

        if P1 > self._ref_psat:
            self._liquid_state.update(CoolProp.PT_INPUTS, P1, T1)
            return (
                self._liquid_state.rhomass(),
                self._liquid_state.hmass(),
                self._liquid_state.smass(),
            )

        dp = P1 - self._ref_psat
        return (
            self._ref_rho + self._ref_drho_dp * dp,
            self._ref_h + dp / self._ref_rho,
            self._ref_s,
        )

    def evaluate_states(self):
        P1 = self.upstream_pressure.value
        P2 = self.downstream_pressure.value
        T1 = self.upstream_temperature.value
        Cd = self.discharge_coefficient.value
        A = self.cross_sectional_area.value

        dP = P1 - P2

        if dP <= 0.0:
            rho1, _, _ = self._upstream_liquid(max(P2, 1.0), T1)
            self.upstream_density.value = rho1
            self.G_SPI.value = -Cd * (2.0 * rho1 * abs(dP)) ** 0.5
            self.G_HEM.value = 0.0
            self.kappa.value = 0.0
            self.G_NHNE.value = self.G_SPI.value
            self.mass_flow.value = self.G_NHNE.value * A
            return

        if self._ref_T != T1:
            self._refresh_reference(T1)
        Psat = self._ref_psat

        rho1, h1, s1 = self._upstream_liquid(P1, T1)
        self.upstream_density.value = rho1

        self.G_SPI.value = Cd * (2.0 * rho1 * dP) ** 0.5

        if P2 >= Psat:
            self.G_HEM.value = 0.0
            self.kappa.value = 0.0
            self.G_NHNE.value = self.G_SPI.value
        else:
            G_hem = 0.0
            for p in np.linspace(max(0.999 * Psat, P2), P2, 60):
                self._flash_state.update(CoolProp.PSmass_INPUTS, p, s1)
                h2 = self._flash_state.hmass()
                rho2 = self._flash_state.rhomass()
                G_hem = max(G_hem, Cd * rho2 * (2.0 * max(h1 - h2, 0.0)) ** 0.5)
            self.G_HEM.value = G_hem

            kappa = (dP / (Psat - P2)) ** 0.5
            self.kappa.value = kappa
            self.G_NHNE.value = (
                kappa * self.G_SPI.value + self.G_HEM.value
            ) / (1.0 + kappa)

        self.mass_flow.value = self.G_NHNE.value * A


Injector = Network("Two Phase Injector")

ox_injector_liquid_mass_flow = State(0.0)
ox_feed_mass_flow = State(0.0)
ox_injector_pressure = State(ox_source_pressure)

ox_feed_density = _sat_ref.rhomass()

OxFeed = DischargeCoefficient(
    "Ox Feed",
    Injector,
    upstream_pressure=ox_source_pressure,
    downstream_pressure=ox_injector_pressure,
    density=ox_feed_density,
    discharge_coefficient=1.0,
    cross_sectional_area=ox_feed_cda,
    mass_flow=ox_feed_mass_flow,
)

OxInjectorLiquidOutlet = NHNEOrifice(
    "Ox Injector Liquid Outlet",
    Injector,
    fluid="NitrousOxide",
    upstream_pressure=ox_injector_pressure,
    upstream_temperature=ox_temperature,
    downstream_pressure=chamber_pressure,
    discharge_coefficient=ox_injector_nhne_cd,
    cross_sectional_area=ox_injector_cda,
    mass_flow=ox_injector_liquid_mass_flow,
)

OxInjectorInlet = Volume(
    "Ox Injector Inlet",
    Injector,
    pressure=ox_injector_pressure,
    mass_flow_in=OxFeed.mass_flow,
    mass_flow_out=OxInjectorLiquidOutlet.mass_flow,
)

chamber_pressures = np.linspace(0.99 * ox_source_pressure, 2.0 * bar_to_pa, 80)

solver = SteadyState(Injector)

injector_dp = []
mass_flux_nhne = []
mass_flux_spi = []
mass_flux_hem = []

for pressure in chamber_pressures:
    chamber_pressure.value = pressure
    solver.solve()
    injector_dp.append((ox_injector_pressure.value - chamber_pressure.value) / bar_to_pa)
    mass_flux_nhne.append(OxInjectorLiquidOutlet.G_NHNE.value)
    mass_flux_spi.append(OxInjectorLiquidOutlet.G_SPI.value)
    mass_flux_hem.append(OxInjectorLiquidOutlet.G_HEM.value)

chamber_pressure.value = design_chamber_pressure
solver.solve()
design_dp = (ox_injector_pressure.value - chamber_pressure.value) / bar_to_pa
design_mass_flux = OxInjectorLiquidOutlet.G_NHNE.value

fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(injector_dp, mass_flux_nhne, label="NHNE")
ax.plot(injector_dp, mass_flux_spi, "--", label="SPI")
ax.plot(injector_dp, mass_flux_hem, "--", label="HEM")
ax.axvline(design_dp, color="k", linestyle=":", label="Pc")
ax.plot(design_dp, design_mass_flux, "ko")
ax.set_xlabel("dP bar")
ax.set_ylabel("Mass flux kg/m2/s")
ax.set_title("Pressure drop over injector to mass flux, T=15.5 C")
ax.grid(True)
ax.legend()
fig.savefig("injector_dp_vs_mass_flux.png", dpi=200, bbox_inches="tight")
plt.show()
