[English language version](Readme.md)
## 質量分析データ処理API ドキュメント

## 📋 目次
1. [クイックスタート](#クイックスタート)
2. [API概要](#api概要)
3. [セットアップガイド](#セットアップガイド)
4. [API仕様](#api仕様)
5. [言語別実装例](#言語別実装例)
6. [トラブルシューティング](#トラブルシューティング)

---

## 🚀 クイックスタート

### 基本的な使用フロー

```
1. サンプル作成        → POST /io_create_sample
2. スペクトル追加    → POST /io_add_scan
3. 注釈追加（任意）    → POST /io_add_annotation
4. データ保存         → POST /io_flush
5. データ抽出         → POST /io_get_spectrum
```

### 最小限の動作確認

```bash
# 1. サンプル作成
curl -X POST http://localhost:8191/io_create_sample \
  -H "Content-Type: application/json" \
  -d "null"

# 2. スペクトル総数確認
curl -X POST http://localhost:8191/io_get_spectra_count \
  -H "Content-Type: application/json" \
  -d "null"
```

---

## 📖 API概要

### 基本情報
- **ベースURL**: `http://localhost:8191/`
- **プロトコル**: HTTP REST API
- **メソッド**: POST（全エンドポイント）
- **データ形式**: JSON
- **認証**: 不要
- **文字エンコーディング**: UTF-8

### 主な機能
- 質量分析スペクトルデータの登録
- ピーク注釈の追加
- データの永続化
- スペクトルデータの抽出


---

## セットアップガイド

### サーバー要件
- **Mass++サーバー**: ポート8191で稼働
- **OS**: Windows, macOS, Linux対応

### 接続確認
```bash
# サーバー稼働確認
curl -X POST http://localhost:8191/io_get_spectra_count \
  -H "Content-Type: application/json" \
  -d "null"
```

成功時のレスポンス:
```json
{"count": 0}
```

---

## API仕様

### 共通仕様

#### HTTPヘッダー
```
POST /endpoint
Content-Type: application/json
```

#### レスポンス形式
- **成功**: HTTPステータス200 + JSONデータ
- **エラー**: HTTPステータス400-500 + エラーメッセージ

#### エラーハンドリング
```json
{
  "error": "Error message",
  "code": "ERROR_CODE"
}
```

---

### データ管理API

#### 1. サンプル作成

```http
POST /io_create_sample
Content-Type: application/json

null
```

**レスポンス**
```json
{
  "id": "sample_12345"
}
```

---

#### 2. スキャンデータ追加

```http
POST /io_add_scan
Content-Type: application/json
```

**リクエストボディ**
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

**パラメータ詳細**

| フィールド | 型 | 必須 | 説明 | 例 |
|-----------|------|------|------|-----|
| `id` | string | ○ | サンプルID | "sample_12345" |
| `msLevel` | integer | ○ | MSレベル | 1, 2, 3... |
| `precursorMz` | number | ○ | プリカーサーm/z | 123.45 (MS1は-1.0) |
| `rt` | number | ○ | 保持時間（秒） | 60.5 |
| `points` | array | ○ | スペクトルデータ | [{"x": m/z, "y": intensity}] |
| `centroidMode` | boolean | ○ | セントロイドモード | true/false |
| `minMz` | number | ○ | 表示用最小m/z | 99.0 |
| `maxMz` | number | ○ | 表示用最大m/z | 200.0 |

**レスポンス**
```json
{
  "status": "success"
}
```

---

#### 3. ピーク注釈追加

```http
POST /io_add_annotation
Content-Type: application/json
```

**リクエストボディ**
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

#### 4. データ保存

```http
POST /io_flush
Content-Type: application/json
```

**リクエストボディ**
```json
{
  "id": "sample_12345",
  "index": 0
}
```

---

### 📤 データ抽出API

#### 5. スペクトル総数取得

```http
POST /io_get_spectra_count
Content-Type: application/json

null
```

**レスポンス**
```json
{
  "count": 150
}
```

---

#### 6. 現在のインデックス取得

```http
POST /io_get_current_index
Content-Type: application/json

null
```

**レスポンス**
```json
{
  "index": 5
}
```

---

#### 7. スペクトルデータ取得

```http
POST /io_get_spectrum
Content-Type: application/json
```

**リクエストボディ**
```json
{
  "index": "0"
}
```

**レスポンス**
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

## 💻 言語別実装例

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

# サンプル作成
sample = call_api('io_create_sample')
sample_id = sample['id']

# スペクトル追加
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
        
        // JSONレスポンスからIDを抽出
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

// 使用例
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

// 使用例
(async () => {
    const api = new MassSpecAPI();
    const sampleId = await api.createSample();
    await api.addScan(sampleId, [100.0, 101.0, 102.0], [1000, 500, 750]);
    await api.callAPI('io_flush', { id: sampleId, index: 0 });
})();
```

---

## 一般的なワークフロー

### データ登録フロー
```
1. サンプル作成 (/io_create_sample)
   ↓
2. スペクトル追加 (/io_add_scan) × N回
   ↓
3. 注釈追加 (/io_add_annotation) ※任意
   ↓
4. データ保存 (/io_flush)
```

### データ抽出フロー
```
1. 総数確認 (/io_get_spectra_count)
   ↓
2. スペクトル取得 (/io_get_spectrum) × N回
   ↓
3. データ分析・エクスポート
```

---

## トラブルシューティング

### よくあるエラー

#### 接続エラー
```
Error: Connection refused
```
**解決法**: Mass++サーバーの起動確認、ポート8191の確認

#### JSON形式エラー
```
HTTP 400: Invalid JSON
```
**解決法**: リクエストボディの形式確認、Content-Typeヘッダーの確認

#### 空のレスポンス
```
HTTP 200: Empty response
```
**解決法**: APIエンドポイントの確認、リクエストパラメータの検証

### デバッグ方法

#### 1. 基本接続テスト
```bash
curl -v -X POST http://localhost:8191/io_get_spectra_count \
  -H "Content-Type: application/json" \
  -d "null"
```

#### 2. レスポンス確認
```bash
# 詳細なHTTPレスポンス表示
curl -i -X POST http://localhost:8191/io_create_sample \
  -H "Content-Type: application/json" \
  -d "null"
```

---

## 📞 サポート

### API仕様
- **Version**: 1.0
- **最終更新**: 2025年7月

### 注意事項
- 全エンドポイントはPOSTメソッドを使用
- リクエストボディが不要な場合も`null`を送信
- 数値の文字列変換に注意（indexパラメータなど）
- Base64エンコード画像データは大きなサイズになる可能性あり

---

*このAPIは質量分析データの効率的な処理を目的として設計されています。各言語での実装例を参考に、お使いの環境に適したクライアントを開発してください。*
