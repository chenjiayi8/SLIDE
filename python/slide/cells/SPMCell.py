import numpy as np
from typing import Tuple
from .Cell import Cell

class SPMCell(Cell):
    def __init__(self,
                 capacity: float,
                 soc: float,
                 Vmax: float,
                 Vmin: float,
                 Ichargemax: float,
                 Idischargemax: float,
                 T: float,
                 deg: float,
                 Rn: float,
                 Rp: float,
                 Dn: float,
                 Dp: float,
                 Rser: float):
        super().__init__(capacity, soc, Vmax, Vmin, Ichargemax, Idischargemax, T, deg)
        self.Rn = Rn  # Negative particle radius [m]
        self.Rp = Rp  # Positive particle radius [m]
        self.Dn = Dn  # Negative diffusion coefficient [m²/s]
        self.Dp = Dp  # Positive diffusion coefficient [m²/s]
        self.Rser = Rser  # Series resistance [Ohm]

    def get_current(self, V: float, I_prev: float, dt: float) -> float:
        # Simplified SPM current calculation with voltage constraints
        OCV = self._calculate_ocv()
        I = (V - OCV) / self.Rser
        
        # Apply current limits
        I = np.clip(I, -self.Ichargemax, self.Idischargemax)
        return I

    def update_voltage(self, I: float, dt: float) -> Tuple[float, float]:
        # Update SOC using coulomb counting
        self.soc -= I * dt / (self.capacity * 3600)  # Convert dt from seconds to hours
        self.soc = np.clip(self.soc, 0, 1)

        # Calculate surface concentrations
        csn_surf = self._calc_surface_concentration(I, self.Rn, self.Dn, 'negative')
        csp_surf = self._calc_surface_concentration(I, self.Rp, self.Dp, 'positive')

        # Calculate voltage components
        OCV = self._calculate_ocv(csn_surf, csp_surf)
        V = OCV - I * self.Rser
        
        return V, OCV

    def thermal_model(self, I: float, T_env: float, dt: float) -> float:
        # Simple thermal model with Joule heating and convection
        Q_joule = I**2 * self.Rser
        self.T += (Q_joule - (self.T - T_env)/0.1) * dt  # 0.1 is thermal resistance
        return self.T

    def _calculate_ocv(self, csn: float = None, csp: float = None) -> float:
        # Simplified OCV curve (should match C++ implementation)
        return 3.7 + 0.5*(self.soc - 0.5) - 0.1*(1 - self.soc)

    def _calc_surface_concentration(self, I: float, R: float, D: float, electrode: str) -> float:
        # Simplified surface concentration calculation
        F = 96485  # Faraday constant [C/mol]
        a = 3 * (1 - 0.5) / R  # Assume porosity of 0.5
        j = I / (a * self.capacity)  # Current density
        
        if electrode == 'negative':
            return 0.5 - j * R**2 / (5 * D * F)
        else:
            return 0.5 + j * R**2 / (5 * D * F)