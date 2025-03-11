from enum import Enum
from typing import NewType, Tuple

# Basic numeric types 
voltage_t = NewType('voltage_t', float)
current_t = NewType('current_t', float)
soc_t = NewType('soc_t', float)
temperature_t = NewType('temperature_t', float)
percent_t = NewType('percent_t', float)

# State vector for battery cells
class CellState:
    def __init__(self, voltage: voltage_t, soc: soc_t, temperature: temperature_t):
        self.voltage = voltage
        self.soc = soc
        self.temperature = temperature

# Battery operational states
class State(Enum):
    RESTING = 0
    CHARGING = 1
    DISCHARGING = 2
    FAULT = 3

# Measurement data structure
class Measurement:
    def __init__(self, timestamp: float, voltage: voltage_t, current: current_t):
        self.timestamp = timestamp
        self.voltage = voltage
        self.current = current