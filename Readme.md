## Document for Mass Spectrometry Data Processing API

## 📋 Contents
1. [Quickstart](#Quickstart)
2. [API Overview](#api-overview)
3. [Setup Guide](#setup-guide)
4. [API Specifications](#api-specifications)
5. [Implementation examples by language](#implementation-examples-by-language)
6. [Troubleshooting](#troubleshooting)

---

## 🚀 Quickstart

### Basic usage flow

```
1. Sample creation        -> POST /io_create_sample
2. Adding Spectrum    -> POST /io_add_scan
3. Adding notes (optional)    -> POST /io_add_annotation
4. Saving data         -> POST /io_flush
5. Data Retrieval         -> POST /io_get_spectrum
```

### Minimal operation check

```bash
# 1. Sample creation
curl -X POST http://localhost:8191/io_create_sample \
  -H "Content-Type: application/json" \
  -d "null"

# 2. Check the total spectrum number
curl -X POST http://localhost:8191/io_get_spectra_count \
  -H "Content-Type: application/json" \
  -d "null"
```

---

## 📖 API Overview

### Basic Information
- **Base URL**: `http://localhost:8191/`
- **Protocol**: HTTP REST API
- **Method**: POST（All endpoints）
- **Data Format**: JSON
- **Authentication**: Unnecessary
- **Character encoding**: UTF-8

### Key Features
- Mass Spectral Data Submission
- Adding Peak Annotations
- Data Persistence
- Retrieving Spectrum Data


---

## Setup Guide

### Server requirements
- **Mass++ Server**: Running on port 8191
- **OS**: Compatible with Windows, macOS, and Linux

### Connection confirmation
```bash
# Server operation check
curl -X POST http://localhost:8191/io_get_spectra_count \
  -H "Content-Type: application/json" \
  -d "null"
```

Success response:
```json
{"count": 0}
```

---

## API Specifications

### Common Specifications

#### HTTP header
```
POST /endpoint
Content-Type: application/json
```

#### Response Format
- **Success**: HTTP status 200 + JSON data
- **Error**: HTTP status 400-500 + error message

#### Error Handling
```json
{
  "error": "Error message",
  "code": "ERROR_CODE"
}
```

---

### Data Management API

#### 1. Sample creation

```http
POST /io_create_sample
Content-Type: application/json

null
```

**Response**
```json
{
  "id": "sample_12345"
}
```

---

#### 2. Adding scan data

```http
POST /io_add_scan
Content-Type: application/json
```

**Request body**
```json
{
  "id": "sample_12345",
  "msLevel": 1,
  "precursorMz": -1.0,
  "rt": 123.45,
  "points": [
    {"x": 100.0, "y": 1000.0},
    {"x": 101.0, "y": 500.0}
  ],
  "centroidMode": true,
  "minMz": 99.5,
  "maxMz": 101.5
}
```

**Parameter details**

| Field | Type | Required | Description | Example |
|-----------|------|------|------|-----|
| `id` | string | Yes | Sample ID | "sample_12345" |
| `msLevel` | integer | Yes | MS level | 1, 2, 3... |
| `precursorMz` | number | Yes | precursor m/z | 123.45 (-1.0 for MS1) |
| `rt` | number | Yes | Retention time (sec) | 60.5 |
| `points` | array | Yes | Spectrum data | [{"x": m/z, "y": intensity}] |
| `centroidMode` | boolean | Yes | Centroid mode | true/false |
| `minMz` | number | Yes | Minimum m/z for display | 99.0 |
| `maxMz` | number | Yes | Maximum m/z for display | 200.0 |

**Response**
```json
{
  "status": "success"
}
```

---

#### 3. Adding peak annotation

```http
POST /io_add_annotation
Content-Type: application/json
```

**Request body**
```json
[
  {
    "name": "C6H6O+",
    "mass": 94.0419,
    "intensity": 1000.0,
    "image": "base64_encoded_image",
    "id": "sample_12345"
  }
]
```

---

#### 4. Saving data

```http
POST /io_flush
Content-Type: application/json
```

**Request body**
```json
{
  "id": "sample_12345",
  "index": 0
}
```

---

### 📤 Data Retrieval API

#### 5. Obtaining the total number of spectra

```http
POST /io_get_spectra_count
Content-Type: application/json

null
```

**Response**
```json
{
  "count": 150
}
```

---

#### 6. Obtaining the current index

```http
POST /io_get_current_index
Content-Type: application/json

null
```

**Response**
```json
{
  "index": 5
}
```

---

#### 7. Obtaining spectrum data

```http
POST /io_get_spectrum
Content-Type: application/json
```

**Request body**
```json
{
  "index": "0"
}
```

**Response**
```json
{
  "id": "sample_12345",
  "msLevel": 1,
  "rt": 123.45,
  "precursorMz": -1.0,
  "minMz": 99.5,
  "maxMz": 200.8,
  "points": [
    {"x": 100.0, "y": 1000.0},
    {"x": 101.0, "y": 500.0}
  ]
}
```

---

## 💻 Implementation examples by language

### 🐍 Python

```python
import requests
import json

def call_api(endpoint, data=None):
    url = f'http://localhost:8191/{endpoint}'
    headers = {'Content-Type': 'application/json'}
    body = json.dumps(data) if data else 'null'
    response = requests.post(url, headers=headers, data=body)
    return response.json()

# Sample creation
sample = call_api('io_create_sample')
sample_id = sample['id']

# Adding Spectrum
scan_data = {
    'id': sample_id,
    'msLevel': 1,
    'precursorMz': -1.0,
    'rt': 60.0,
    'points': [{'x': 100.0, 'y': 1000.0}],
    'centroidMode': True,
    'minMz': 99.0,
    'maxMz': 101.0
}
call_api('io_add_scan', scan_data)
call_api('io_flush', {'id': sample_id, 'index': 0})
```

### ☕ Java

```java
import java.net.http.HttpClient;
import java.net.http.HttpRequest;
import java.net.http.HttpResponse;
import java.net.URI;
import com.fasterxml.jackson.databind.ObjectMapper;

public class MassSpecAPI {
    private static final String BASE_URL = "http://localhost:8191/";
    private final HttpClient client = HttpClient.newHttpClient();
    private final ObjectMapper mapper = new ObjectMapper();

    public String createSample() throws Exception {
        HttpRequest request = HttpRequest.newBuilder()
                .uri(URI.create(BASE_URL + "io_create_sample"))
                .header("Content-Type", "application/json")
                .POST(HttpRequest.BodyPublishers.ofString("null"))
                .build();

        HttpResponse<String> response = client.send(request, 
                HttpResponse.BodyHandlers.ofString());
        
        // Extracting IDs from JSON responses
        return mapper.readTree(response.body()).get("id").asText();
    }

    public void addScan(String sampleId, double[] mzValues, double[] intensities) throws Exception {
        StringBuilder pointsJson = new StringBuilder("[");
        for (int i = 0; i < mzValues.length; i++) {
            if (i > 0) pointsJson.append(",");
            pointsJson.append(String.format("{\"x\":%.3f,\"y\":%.1f}", 
                    mzValues[i], intensities[i]));
        }
        pointsJson.append("]");

        String json = String.format(
            "{\"id\":\"%s\",\"msLevel\":1,\"precursorMz\":-1.0," +
            "\"rt\":60.0,\"points\":%s,\"centroidMode\":true," +
            "\"minMz\":%.1f,\"maxMz\":%.1f}",
            sampleId, pointsJson.toString(), 
            mzValues[0] - 1, mzValues[mzValues.length - 1] + 1
        );

        HttpRequest request = HttpRequest.newBuilder()
                .uri(URI.create(BASE_URL + "io_add_scan"))
                .header("Content-Type", "application/json")
                .POST(HttpRequest.BodyPublishers.ofString(json))
                .build();

        client.send(request, HttpResponse.BodyHandlers.ofString());
    }
}
```

### 🐘 PHP

```php
<?php
class MassSpecAPI {
    private $baseUrl = 'http://localhost:8191/';

    private function callAPI($endpoint, $data = null) {
        $url = $this->baseUrl . $endpoint;
        $json = $data ? json_encode($data) : 'null';
        
        $ch = curl_init();
        curl_setopt($ch, CURLOPT_URL, $url);
        curl_setopt($ch, CURLOPT_POST, true);
        curl_setopt($ch, CURLOPT_POSTFIELDS, $json);
        curl_setopt($ch, CURLOPT_HTTPHEADER, ['Content-Type: application/json']);
        curl_setopt($ch, CURLOPT_RETURNTRANSFER, true);
        
        $response = curl_exec($ch);
        curl_close($ch);
        
        return json_decode($response, true);
    }

    public function createSample() {
        $result = $this->callAPI('io_create_sample');
        return $result['id'];
    }

    public function addScan($sampleId, $mzValues, $intensities) {
        $points = [];
        for ($i = 0; $i < count($mzValues); $i++) {
            $points[] = ['x' => $mzValues[$i], 'y' => $intensities[$i]];
        }

        $scanData = [
            'id' => $sampleId,
            'msLevel' => 1,
            'precursorMz' => -1.0,
            'rt' => 60.0,
            'points' => $points,
            'centroidMode' => true,
            'minMz' => min($mzValues) - 1,
            'maxMz' => max($mzValues) + 1
        ];

        return $this->callAPI('io_add_scan', $scanData);
    }

    public function flushData($sampleId) {
        return $this->callAPI('io_flush', ['id' => $sampleId, 'index' => 0]);
    }
}

// Usage example
$api = new MassSpecAPI();
$sampleId = $api->createSample();
$api->addScan($sampleId, [100.0, 101.0, 102.0], [1000, 500, 750]);
$api->flushData($sampleId);
?>
```

### JavaScript (Node.js)

```javascript
const axios = require('axios');

class MassSpecAPI {
    constructor(baseUrl = 'http://localhost:8191/') {
        this.baseUrl = baseUrl;
        this.headers = { 'Content-Type': 'application/json' };
    }

    async callAPI(endpoint, data = null) {
        const url = `${this.baseUrl}${endpoint}`;
        const body = data ? JSON.stringify(data) : 'null';
        const response = await axios.post(url, body, { headers: this.headers });
        return response.data;
    }

    async createSample() {
        const result = await this.callAPI('io_create_sample');
        return result.id;
    }

    async addScan(sampleId, mzValues, intensities) {
        const points = mzValues.map((mz, i) => ({ x: mz, y: intensities[i] }));
        
        const scanData = {
            id: sampleId,
            msLevel: 1,
            precursorMz: -1.0,
            rt: 60.0,
            points: points,
            centroidMode: true,
            minMz: Math.min(...mzValues) - 1,
            maxMz: Math.max(...mzValues) + 1
        };

        return await this.callAPI('io_add_scan', scanData);
    }
}

// Usage example
(async () => {
    const api = new MassSpecAPI();
    const sampleId = await api.createSample();
    await api.addScan(sampleId, [100.0, 101.0, 102.0], [1000, 500, 750]);
    await api.callAPI('io_flush', { id: sampleId, index: 0 });
})();
```

---

## Typical workflow

### Data registration flow
```
1. Sample creation (/io_create_sample)
   ->
2. Adding Spectrum (/io_add_scan) * N times
   ->
3. Adding annotation (/io_add_annotation) # Optional
   ->
4. Saving data (/io_flush)
```

### Data Extraction Flow
```
1. Check the total number (/io_get_spectra_count)
   ->
2. Obtaining spectrum (/io_get_spectrum) * N times
   ->
3. Data analysis and export
```

---

## Troubleshooting

### Common errors

#### Connection Error
```
Error: Connection refused
```
**Solution**: Check that the Mass++ server is running and that port 8191 is running

#### JSON format error
```
HTTP 400: Invalid JSON
```
**Solution**: Check the format of the request body and the Content-Type header

#### Empty response
```
HTTP 200: Empty response
```
**Solution**: Check API endpoints and validate request parameters

### How to debug

#### 1. Basic Connection Test
```bash
curl -v -X POST http://localhost:8191/io_get_spectra_count \
  -H "Content-Type: application/json" \
  -d "null"
```

#### 2. Response confirmation
```bash
# Detailed HTTP response display
curl -i -X POST http://localhost:8191/io_create_sample \
  -H "Content-Type: application/json" \
  -d "null"
```

---

## 📞 Support

### API Specifications
- **Version**: 1.0
- **Last update**: July, 2025

### Notes
- All endpoints use the POST method.
- Send `null` even if no request body is required.
- Be careful when converting numbers to strings (such as index parameters).
- Base64 encoded image data can be large in size.

---

*This API is designed for efficient processing of mass spectrometry data. Please refer to the implementation examples in each language to develop a client suitable for your environment.**
