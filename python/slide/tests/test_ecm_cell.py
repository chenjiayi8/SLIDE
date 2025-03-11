import numpy as np
import pytest
import numpy as np
from slide.cells.ECMCell import ECMCell, ECMParameters

@pytest.fixture
def basic_ecm_params():
    return ECMParameters(
        R0=0.01,  # Ohms
        R1=0.02,  # Ohms
        C1=3000,  # Farads
        OCV=np.array([3.0, 3.5, 4.0]),  # OCV curve for 0%, 50%, 100% SOC
        capacity=5.0,  # Ah
        mass=0.045,  # kg
        surface_area=0.005,  # m²
        heat_capacity=900,  # J/(kg·K)
        thermal_conductivity=0.5  # W/(m·K)
    )

def test_initialization(basic_ecm_params):
    cell = ECMCell(basic_ecm_params, 0.5, 25.0)
    
    assert cell.state_of_charge == 0.5
    assert cell.core_temperature == 25.0
    assert cell.voltage == pytest.approx(3.6, abs=0.1)

def test_step_method(basic_ecm_params):
    cell = ECMCell(basic_ecm_params, 0.5, 25.0)
    
    # Test charge and discharge scenarios
    for current, dt in [(-5.0, 1800), (2.5, 7200)]:
        initial_soc = cell.state_of_charge
        initial_v = cell.voltage
        
        cell.step(current=current, dt=dt)
        
        # Validate SOC change
        expected_soc = initial_soc - (current * dt) / (basic_ecm_params.capacity * 3600)
        expected_soc = max(0.0, min(1.0, expected_soc))
        assert cell.state_of_charge == pytest.approx(expected_soc, abs=1e-4)
        
        # Validate voltage calculation
        soc_index = int(expected_soc * (len(basic_ecm_params.OCV)-1))
        expected_ocv = basic_ecm_params.OCV[soc_index]
        tau = basic_ecm_params.R1 * basic_ecm_params.C1
        v_r1 = current * basic_ecm_params.R1 * (1 - np.exp(-dt/tau))
        expected_v = expected_ocv - current*basic_ecm_params.R0 - v_r1
        assert cell.voltage == pytest.approx(expected_v, rel=1e-3)

def test_thermal_behavior(basic_ecm_params):
    cell = ECMCell(basic_ecm_params, 0.5, 25.0)
    ambient = 25.0
    
    # Apply 10A current for 10 seconds
    cell.step(10.0, 10.0)
    cell.update_thermal(ambient, 10.0)
    
    # Calculate expected temperature rise
    i_sq_r = 10**2 * basic_ecm_params.R0
    heat = i_sq_r * 10.0
    delta = heat / (basic_ecm_params.mass * basic_ecm_params.heat_capacity)
    assert cell.core_temperature == pytest.approx(25.0 + delta, abs=0.1)

def test_soc_boundaries(basic_ecm_params):
    cell = ECMCell(basic_ecm_params, 0.05, 25.0)
    
    # Try to over-discharge
    cell.step(5.0, 4000)  # 5A for 4000s (5*4000/18000 = 1.111)
    assert cell.state_of_charge == 0.0
    
    # Try to over-charge
    cell.soc = 1.0
    cell.step(-5.0, 4000)
    assert cell.state_of_charge == 1.0

def test_ocv_interpolation(basic_ecm_params):
    # Test OCV at various SOC points
    cell = ECMCell(basic_ecm_params, 0.0, 25.0)
    assert cell.voltage == pytest.approx(3.0, abs=0.01)
    
    cell.soc = 0.5
    cell._update_ocv()
    assert cell.voltage == pytest.approx(3.6, abs=0.01)
    
    cell.soc = 1.0
    cell._update_ocv()
    assert cell.voltage == pytest.approx(4.2, abs=0.01)