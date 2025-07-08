from common_ms_utils import (
    MassSpecDataProcessor, 
    call_service, 
    get_nearest_point, 
    normalize_formula
)


def extract_from_masspp():
    result = call_service('io_get_spectra_count')
    count = int(result.get('count', 0))
    print(f'Total spectra count: {count}')

    result = call_service('io_get_current_index')
    index = int(result.get('index', 0))
    print(f'Current index: {index}')

    for i in range(count):
        request = {'index': str(i)}
        spectrum = call_service('io_get_spectrum', request)

        id = spectrum.get('id', '')
        msLevel = int(spectrum.get('msLevel', 1))
        rt = float(spectrum.get('rt', -1.0))
        precursorMz = float(spectrum.get('precursorMz', -1.0))
        minMz = float(spectrum.get('minMz', 0.0))
        maxMz = float(spectrum.get('maxMz', 0.0))
        points = spectrum.get('points', [])

        print(f'Spectrum ID: {id}, MS Level: {msLevel}, RT: {rt}, Precursor MZ: {precursorMz}, Min MZ: {minMz}, Max MZ: {maxMz}')
        if not points:
            print('    No points in spectrum')
            continue
        print(f'    Number of points: {len(points)}')



if __name__ == "__main__":
    extract_from_masspp()