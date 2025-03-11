from abc import ABC, abstractmethod
from typing import Tuple

class Cell(ABC):
    def __init__(self, 
                 capacity: float,
                 soc: float,
                 Vmax: float,
                 Vmin: float,
                 Ichargemax: float,
                 Idischargemax: float,
                 T: float,
                 deg: float):
        self.capacity = capacity
        self.soc = soc
        self.Vmax = Vmax
        self.Vmin = Vmin
        self.Ichargemax = Ichargemax
        self.Idischargemax = Idischargemax
        self.T = T
        self.deg = deg

    @abstractmethod
    def get_current(self, V: float, I_prev: float, dt: float) -> float:
        pass

    @abstractmethod
    def update_voltage(self, I: float, dt: float) -> Tuple[float, float]:
        pass

    @abstractmethod
    def thermal_model(self, I: float, T_env: float, dt: float) -> float:
        pass

    @property
    def SOC(self) -> float:
        return self.soc

    @property
    def degradation(self) -> float:
        return self.deg