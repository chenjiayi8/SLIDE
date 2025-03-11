import numpy as np
from abc import ABC, abstractmethod

class StorageUnit(ABC):
    """
    Base class for all energy storage units
    Translated from C++ StorageUnit.hpp
    """
    
    def __init__(self, capacity: float, soc: float = 0.5):
        self.capacity = capacity  # Ah
        self.soc = np.clip(soc, 0, 1)
        self.voltage = 0.0
        self.temperature = 298.15  # K
        self.degradation = 0.0
    
    @abstractmethod
    def step(self, current: float, dt: float) -> None:
        """
        Perform time step update
        :param current: Current in A (positive for discharge)
        :param dt: Time step in seconds
        """
        pass
    
    @abstractmethod
    def get_ocv(self) -> float:
        """Get open circuit voltage"""
        pass
    
    def reset(self) -> None:
        """Reset to initial conditions"""
        self.soc = 0.5
        self.voltage = self.get_ocv()
        self.temperature = 298.15
        self.degradation = 0.0
    
    def get_soc(self) -> float:
        return self.soc
    
    def get_voltage(self) -> float:
        return self.voltage
    
    def get_degradation(self) -> float:
        return self.degradation