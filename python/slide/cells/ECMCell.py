from slide.types import (
    current_t, voltage_t, soc_t,
    temperature_t, State, Measurement
)
import numpy as np
from ..system.StorageUnit import StorageUnit

class ECMCell(StorageUnit):
    """
    Equivalent Circuit Model battery cell
    Translated from C++ ECMCell.hpp
    """
    
    def __init__(self,
                 capacity: float,
                 r_series: float,
                 r_parallel: float,
                 c_parallel: float,
                 soc: soc_t = soc_t(0.5)):
        super().__init__(capacity, soc)
        self.r_series = r_series  # Ohms
        self.r_parallel = r_parallel  # Ohms
        self.c_parallel = c_parallel  # Farads
        self.v_c = voltage_t(0.0)  # Capacitor voltage
        self._state = State.RESTING
        
    def step(self, current: current_t, dt: float) -> Measurement:
        # SOC update
        delta_soc = (current * dt) / (3600 * self.capacity)
        self.soc = np.clip(self.soc - delta_soc, 0, 1)
        
        # RC network update
        tau = self.r_parallel * self.c_parallel
        self.v_c += (current - self.v_c/self.r_parallel) * dt / tau
        
        # Voltage calculation
        ocv = self.get_ocv()
        self.voltage = ocv - current*self.r_series - self.v_c
        
        # Simple degradation model
        self.degradation += abs(current) * dt / (3600 * self.capacity * 1000)
        
    def get_ocv(self):
        # Simplified OCV curve (should be replaced with real data)
        return 3.7 - 0.5*(1 - self.soc)
        
    def reset(self):
        super().reset()
        self.v_c = 0.0