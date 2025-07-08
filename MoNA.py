import requests
import re
import io
import json
import base64
from typing import List, Dict, Any, Tuple

from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from PIL import Image

from common_ms_utils import (
    MassSpecDataProcessor, 
    call_service, 
    get_nearest_point, 
    normalize_formula
)


# Configuration
IMAGE_SIZE = (100, 100)


def get_matched_smiles(formula: str, smiles: List[str]) -> List[str]:
    matched = []
    name = normalize_formula(formula)
    for smile in smiles:
        mol = Chem.MolFromSmiles(smile)
        if mol is not None:
            current_formula = CalcMolFormula(mol)
            if name == normalize_formula(current_formula):
                matched.append(smile)
    
    return matched[0:1]


def get_base64(smiles: List[str]) -> str:
    images = []
    
    for smile in smiles:
        mol = Chem.MolFromSmiles(smile)
        if mol is not None:
            img = Draw.MolToImage(mol, size=IMAGE_SIZE, kekulize=True, wedgeBonds=True)
            images.append(img)
    
    if not images:
        return ''
    
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


def get_peaks(annotations: List[Dict[str, Any]], sample_id: str, 
              spectrum: List[Tuple[float, float]], smiles: List[str]) -> List[Dict[str, Any]]:
    peaks = []
    for annotation in annotations:
        name = annotation['name']
        value = annotation['value']
        point = get_nearest_point(value, spectrum)
        if point is not None:
            matched_smiles = get_matched_smiles(name, smiles)
            base64_img = ''
            if matched_smiles:
                base64_img = get_base64(matched_smiles)
            
            peak = {
                'name': name,
                'mass': value,
                'intensity': point[1],
                'image': base64_img,
                'id': sample_id,
            }
            peaks.append(peak)
    
    return peaks


def get_smiles(data: Dict[str, Any]) -> List[str]:
    smiles = []
    
    compounds = data.get('compound', [])
    for compound in compounds:
        metadata = compound.get('metaData', [])
        for meta in metadata:
            if meta.get('name') == 'SMILES':
                smiles.append(meta['value'])
    
    return smiles


def get_annotations(data: Dict[str, Any]) -> List[Dict[str, Any]]:
    annotations = []
    
    for annotation in data.get('annotations', []):
        annotations.append({
            'name': annotation['name'],
            'value': float(annotation['value'])
        })
    
    return annotations


def get_spectrum(data: Dict[str, Any]) -> List[Tuple[float, float]]:
    points = []
    spectrum = data.get('spectrum', '')
    pairs = spectrum.split(' ')
    for pair in pairs:
        coordinates = pair.strip().split(':')
        if len(coordinates) == 2:
            try:
                mz = float(coordinates[0])
                intensity = float(coordinates[1])
                points.append((mz, intensity))
            except ValueError:
                continue
    
    return points


def process_mona_spectrum(spectrum_id: str, centroid_mode: bool = True):
    url = f'https://mona.fiehnlab.ucdavis.edu/rest/spectra/{spectrum_id}'
    response = requests.get(url)
    
    if response.status_code != 200:
        raise Exception(f'Failed to fetch spectrum {spectrum_id}: {response.status_code}')
    
    data = response.json()
    
    # Extract data components
    smiles = get_smiles(data)
    annotations = get_annotations(data)
    spectrum = get_spectrum(data)
    
    if not spectrum:
        raise Exception(f'No spectrum data found for {spectrum_id}')
    
    # Process using common utilities
    processor = MassSpecDataProcessor()
    sample_id = processor.create_sample()
    
    # Add spectrum scan
    processor.add_scan(spectrum, centroid_mode=centroid_mode)
    
    # Add peak annotations
    peaks = get_peaks(annotations, sample_id, spectrum, smiles)
    if peaks:
        call_service('io_add_annotation', peaks)
    
    # Flush data
    processor.flush()
    
    print(f'Successfully processed MoNA spectrum {spectrum_id}')


if __name__ == '__main__':
    spectrum_id = 'EA299203'
    try:
        process_mona_spectrum(spectrum_id, centroid_mode=True)
    except Exception as e:
        print(f'Error processing spectrum: {e}')