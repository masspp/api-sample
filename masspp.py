import re
import io
import json
import base64

from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

from PIL import Image


base_url = 'http://localhost:8191/'
image_size = (100, 100)


def call_service(service_name, data):
    import requests
    import json

    url = f'{base_url}{service_name}'
    headers = {'Content-Type': 'application/json'}
    
    response = requests.post(url, headers=headers, data=json.dumps(data))
    
    if response.status_code != 200:
        raise Exception(f'Error calling service {service_name}: {response.text}')
    
    s = response.json()
    print(f'Called service {service_name} with response: {s}')
    data = json.loads(s)
    
    return data


def create_scan(spectrum, sample_id, centroid_mode):
    points = []
    min_mz = float('inf')
    max_mz = float('-inf')
    for mz, intensity in spectrum:
        if mz < min_mz:
            min_mz = mz
        if mz > max_mz:
            max_mz = mz
        points.append({'x': mz, 'y': intensity})

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
        'msLevel': 1,
        'precursorMz': -1.0,
        'rt': -1.0,
        'points': points,
        'minMz': min_mz,
        'maxMz': max_mz,
        'centroidMode': centroid_mode,
    }

    return scan    


def normalize_formula(ion_formula):
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
        f'{element}{formula_dict[element]}' if formula_dict[element] > 0 else '' for element in sorted(formula_dict)
    )
    return neutral_formula


def get_matched_smiles(formula, smiles):
    matched = []
    name = normalize_formula(formula)
    for smile in smiles:
        mol = Chem.MolFromSmiles(smile)
        current_formula = CalcMolFormula(mol)
        if name == normalize_formula(current_formula):
            matched.append(smile)

    return matched[0:1]


def get_base64(smiles):
    images = []

    for smile in smiles:
        mol = Chem.MolFromSmiles(smile)
        img = Draw.MolToImage(mol, size=image_size, kekulize=True, wedgeBonds=True)
        images.append(img)

    combined_height = sum(img.height for img in images)
    combined_width = max(img.width for img in images)
    height = 0

    combined_image = Image.new('RGB', (combined_width, combined_height), color=(255, 255, 255))
    for img in images:
        combined_image.paste(img, (0, height))
        height += img.height

    buffered = io.BytesIO()
    combined_image.save(buffered, format='PNG')
    img_base64 = base64.b64encode(buffered.getvalue()).decode('utf-8')

    return img_base64 


def get_nearest_point(mz, spectrum):
    nearest_point = None
    min_diff = float('inf')

    for point in spectrum:
        diff = abs(point[0] - mz)
        if diff < min_diff:
            min_diff = diff
            nearest_point = point

    return nearest_point


def get_peaks(annotations, id, spectrum, smiles):
    peaks = []
    for annotation in annotations:
        name = annotation['name']
        value = annotation['value']
        point = get_nearest_point(value, spectrum)
        if point is not None:
            matched_smiles = get_matched_smiles(name, smiles)
            base64 = ''
            if matched_smiles is not None and len(matched_smiles) > 0:
                base64 = get_base64(matched_smiles)

            peak = {
                'name': name,
                'mass': value,
                'intensity': point[1],
                'image': base64,
                'id': id,
            }
            peaks.append(peak)

    return peaks


def post_to_masspp(spectrum, annotations, smiles, centroid_mode=False):
    response = call_service('io_create_sample', None)  
    sample_id = response['id']

    scan = create_scan(spectrum, sample_id, centroid_mode)
    call_service('io_add_scan', scan)

    peaks = get_peaks(annotations, sample_id, spectrum, smiles)
    call_service('io_add_annotation', peaks)

    call_service('io_flush', {'id': sample_id, 'index': 0})
