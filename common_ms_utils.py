import requests
import json
import re
import numpy as np
from typing import List, Dict, Any, Tuple, Optional, Union


# Configuration
BASE_URL = 'http://localhost:8191/'


def sanitize_value(value: Any) -> Any:
    if isinstance(value, (np.floating, np.integer)):
        return float(value)
    elif isinstance(value, (list, tuple)):
        return [sanitize_value(v) for v in value]
    elif isinstance(value, dict):
        return {k: sanitize_value(v) for k, v in value.items()}
    elif isinstance(value, str):
        return re.sub(r'[\x00-\x1F\x7F]', '', value)
    elif isinstance(value, bool):
        return value
    elif isinstance(value, int):
        return int(value)
    else:
        try:
            return float(value)
        except:
            return value


def call_service(service_name: str, data: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    url = f'{BASE_URL}{service_name}'
    headers = {'Content-Type': 'application/json'}
    
    if data is not None:
        data = sanitize_value(data)
    
    json_data = json.dumps(data)
    
    response = requests.post(url, headers=headers, data=json_data)
    
    if response.status_code != 200:
        raise Exception(f'Error calling service {service_name}: {response.status_code} {response.text}')
    
    if not response.text.strip():
        raise Exception(f'Service {service_name} returned empty response.')
    
    try:
        result = response.json()
        if isinstance(result, str):
            result = json.loads(result)
    except json.JSONDecodeError as e:
        raise Exception(f'Failed to decode JSON from service {service_name}: {response.text}') from e
    
    return result


def create_scan(spectrum: List[Tuple[float, float]], sample_id: str, 
                centroid_mode: bool = False, ms_level: int = 1, 
                precursor_mz: float = -1.0, rt: float = -1.0) -> Dict[str, Any]:
    points = []
    min_mz = float('inf')
    max_mz = float('-inf')
    
    for mz, intensity in spectrum:
        if mz < min_mz:
            min_mz = mz
        if mz > max_mz:
            max_mz = mz
        points.append({'x': float(mz), 'y': float(intensity)})
    
    # Calculate m/z range with padding
    mz_range = max_mz - min_mz
    mid_mz = (max_mz + min_mz) / 2.0
    if mz_range < 1.0:
        mz_range = 1.0
    mz_range = mz_range * 1.1
    min_mz = mid_mz - mz_range / 2.0
    if min_mz < 1.0:
        min_mz = 1.0
    max_mz = min_mz + mz_range
    
    scan = {
        'id': sample_id,
        'msLevel': ms_level,
        'precursorMz': float(precursor_mz),
        'rt': float(rt),
        'points': points,
        'minMz': float(min_mz),
        'maxMz': float(max_mz),
        'centroidMode': bool(centroid_mode),
    }
    
    return scan


def create_sample() -> str:
    response = call_service('io_create_sample', None)
    return response['id']


def add_scan_to_sample(scan: Dict[str, Any]) -> Dict[str, Any]:
    return call_service('io_add_scan', scan)


def flush_sample(sample_id: str, index: int = 0) -> Dict[str, Any]:
    return call_service('io_flush', {'id': sample_id, 'index': index})


def get_nearest_point(mz: float, spectrum: List[Tuple[float, float]]) -> Optional[Tuple[float, float]]:
    nearest_point = None
    min_diff = float('inf')
    
    for point in spectrum:
        diff = abs(point[0] - mz)
        if diff < min_diff:
            min_diff = diff
            nearest_point = point
    
    return nearest_point


def normalize_formula(ion_formula: str) -> str:
    charge_match = re.search(r'([+=]+)$', ion_formula)
    formula = re.sub(r'[+-]+$', '', ion_formula)
    
    delta_h = 0
    if charge_match:
        charge = charge_match.group(1)
        delta_h = -len(charge) if '+' in charge else len(charge)
    
    elements = re.findall(r'([A-Z][a-z]*)(\d*)', formula)
    formula_dict = {}
    for element, count in elements:
        formula_dict[element] = int(count) if count else 1
    
    formula_dict['H'] = formula_dict.get('H', 0) + delta_h
    if formula_dict['H'] < 0:
        formula_dict['H'] = 0
    
    neutral_formula = ''.join(
        f'{element}{formula_dict[element]}' if formula_dict[element] > 0 else '' 
        for element in sorted(formula_dict)
    )
    return neutral_formula


class MassSpecDataProcessor:
    def __init__(self, base_url: str = BASE_URL):
        self.base_url = base_url
        self.sample_id = None
    
    def create_sample(self) -> str:
        self.sample_id = create_sample()
        return self.sample_id
    
    def add_scan(self, spectrum: List[Tuple[float, float]], 
                 centroid_mode: bool = False, ms_level: int = 1,
                 precursor_mz: float = -1.0, rt: float = -1.0) -> Dict[str, Any]:
        if self.sample_id is None:
            self.create_sample()
        
        scan = create_scan(spectrum, self.sample_id, centroid_mode, ms_level, precursor_mz, rt)
        return add_scan_to_sample(scan)
    
    def flush(self, index: int = 0) -> Dict[str, Any]:
        if self.sample_id is None:
            raise ValueError("No sample created. Call create_sample() first.")
        
        return flush_sample(self.sample_id, index)