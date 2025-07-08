from pyteomics import mzxml

import requests
import json
import sys
import re
import numpy as np

base_url = 'http://localhost:8191/'


def sanitize_value(value):
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
            # fallback: try convert to float (e.g., for strange objects)
            return float(value)
        except:
            return value


def call_service(service_name, data):
    url = f'{base_url}{service_name}'
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
            result = json.loads(result)  # ← ここが重要
    except json.JSONDecodeError as e:
        raise Exception(f'Failed to decode JSON from service {service_name}: {response.text}') from e

    print(f'Called service {service_name} with response: {result}')
    return result



def read_mzxml(file_path):
    response = call_service('io_create_sample', None)  
    sample_id = response['id']

    with mzxml.MzXML(file_path) as reader:
        for spectrum in reader:
            ms_level = spectrum.get('msLevel')
            rt = spectrum.get('retentionTime')
            centroid_mode = spectrum.get('centroided', False)

            mz_array = spectrum.get('m/z array', [])
            intensity_array = spectrum.get('intensity array', [])

            points = []
            min_mz = float('inf')
            max_mz = float('-inf')
            for mz, intensity in zip(mz_array, intensity_array):
                if mz < min_mz:
                    min_mz = mz
                if mz > max_mz:
                    max_mz = mz
                points.append({'x': float(mz), 'y': float(intensity)})

            mz_range = max_mz - min_mz
            mid_mz = (max_mz + min_mz) / 2.0
            if mz_range < 1.0:
                mz_range = 1.0
            mz_range = mz_range * 1.1
            min_mz = mid_mz - mz_range / 2.0
            if min_mz < 1.0:
                min_mz = 1.0
            max_mz = min_mz + mz_range

            precursor_mz = -1.0
            if ms_level >= 2:
                precursor_list = spectrum.get('precursorMz')
                if precursor_list:
                    if isinstance(precursor_list, list):
                        precursor_mz = precursor_list[0].get('precursorMz')
                    elif isinstance(precursor_list, dict):
                        precursor_mz = precursor_list.get('precursorMz')

            scan = {
                'id': sample_id,
                'msLevel': ms_level,
                'precursorMz': precursor_mz,
                'rt': rt,
                'points': points,
                'centroidMode': bool(centroid_mode),
                'minMz': float(min_mz),
                'maxMz': float(max_mz),
            }
            scan = call_service('io_add_scan', scan)

    call_service('io_flush', {'id': sample_id, 'index': 0})


if __name__ == "__main__":
    file_path = sys.argv[1] if len(sys.argv) > 1 else 'example.mzXML'
    read_mzxml(file_path)