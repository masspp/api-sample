import requests
import re



import masspp

def get_smiles(data):
    smiles = []

    compounds = data['compound']
    for compound in compounds:
        metadata = compound['metaData']
        for meta in metadata:
            if meta['name'] == 'SMILES':
                smiles.append(meta['value'])

    return smiles


def get_annotations(data):
    annotations = []

    for annotation in data['annotations']:
        annotations.append(
            {
                'name': annotation['name'],
                'value': float(annotation['value'])
            }
        )

    return annotations


def get_spectrum(data):
    points = []
    spectrum = data['spectrum']
    pairs = spectrum.split(' ')
    for pair in pairs:
        coordinates = pair.strip().split(':')
        if len(coordinates) == 2:
            mz = float(coordinates[0])
            intensity = float(coordinates[1])
            points.append((mz, intensity))

    return points


if __name__ == '__main__':
    spectrum_id = 'EA299203'
    url = f'https://mona.fiehnlab.ucdavis.edu/rest/spectra/{spectrum_id}'
    response = requests.get(url)

    if response.status_code == 200:
        data = response.json()
        smiles = get_smiles(data)
        annotations = get_annotations(data)
        spectrum  = get_spectrum(data)

        masspp.post_to_masspp(spectrum, annotations, smiles, True)

    else:
        print('The specified spectrum ID does not exist or is not accessible. ')