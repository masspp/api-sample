# How to Set Up and Run the Sample Program

## 1. Create the Virtual Environment
```bash
python -m venv venv
```

## 2. Activate the Virtual Environment
```bash
source venv/bin/activate
```

## 3. Install Required libraries
```bash
pip install --upgrade pip
pip install -r requirements.txt
```

## 4. Run the sample program
```bash
python MoNa.py`
```

# 質量分析データ処理API ドキュメント

## 概要
このAPIは質量分析データ（mzXMLファイル）の処理と管理を行うためのRESTful APIです。

**ベースURL:** `http://localhost:8191/`

## 認証
現在のところ認証は不要です。

## 共通仕様

### リクエストヘッダー
```
Content-Type: application/json
```

### レスポンス形式
- 成功時: JSON形式のデータ
- エラー時: HTTPステータスコード + エラーメッセージ

### エラーハンドリング
- HTTPステータスコードが200以外の場合はエラー
- レスポンスボディが空の場合はエラー
- JSON形式でないレスポンスの場合はエラー

## エンドポイント

### 1. サンプル作成

**URL:** `/io_create_sample`  
**メソッド:** POST  
**説明:** 新しいサンプルを作成し、一意のIDを取得します。

#### リクエスト
```json
null
```

#### レスポンス
```json
{
  "id": "sample_id_string_or_number"
}
```

#### 例
```bash
curl -X POST http://localhost:8191/io_create_sample \
  -H "Content-Type: application/json" \
  -d "null"
```

### 2. スキャンデータ追加

**URL:** `/io_add_scan`  
**メソッド:** POST  
**説明:** 指定されたサンプルに質量分析スキャンデータを追加します。

#### リクエスト
```json
{
  "id": "sample_id",
  "msLevel": 1,
  "precursorMz": -1.0,
  "rt": 123.45,
  "points": [
    {
      "x": 100.0,
      "y": 1000.0
    },
    {
      "x": 101.0,
      "y": 500.0
    }
  ],
  "centroidMode": true,
  "minMz": 99.95,
  "maxMz": 201.05
}
```

#### パラメータ
- `id` (string/number): サンプルID
- `msLevel` (integer): MSレベル（1 = MS1, 2 = MS2, etc.）
- `precursorMz` (float): プリカーサーイオンのm/z値（MS1の場合は-1.0）
- `rt` (float): 保持時間（秒）
- `points` (array): スペクトラムデータポイント
  - `x` (float): m/z値
  - `y` (float): 強度値
- `centroidMode` (boolean): セントロイドモードかどうか
- `minMz` (float): 最小m/z値（表示範囲用）
- `maxMz` (float): 最大m/z値（表示範囲用）

#### レスポンス
```json
{
  "status": "success",
  "message": "Scan added successfully"
}
```

### 3. データフラッシュ

**URL:** `/io_flush`  
**メソッド:** POST  
**説明:** 指定されたサンプルのデータをフラッシュ（永続化）します。

#### リクエスト
```json
{
  "id": "sample_id",
  "index": 0
}
```

#### パラメータ
- `id` (string/number): サンプルID
- `index` (integer): インデックス番号（通常は0）

#### レスポンス
```json
{
  "status": "success",
  "message": "Data flushed successfully"
}
```

## データフロー

1. **サンプル作成**: `io_create_sample`でサンプルIDを取得
2. **スキャンデータ追加**: `io_add_scan`で各スペクトラムデータを追加
3. **データフラッシュ**: `io_flush`でデータを永続化

## 使用例

### Python
```python
import requests
import json

base_url = 'http://localhost:8191/'

# 1. サンプル作成
response = requests.post(f'{base_url}io_create_sample', 
                        headers={'Content-Type': 'application/json'}, 
                        data='null')
sample_id = response.json()['id']

# 2. スキャンデータ追加
scan_data = {
    'id': sample_id,
    'msLevel': 1,
    'precursorMz': -1.0,
    'rt': 60.0,
    'points': [
        {'x': 100.0, 'y': 1000.0},
        {'x': 101.0, 'y': 500.0}
    ],
    'centroidMode': True,
    'minMz': 99.0,
    'maxMz': 102.0
}

response = requests.post(f'{base_url}io_add_scan',
                        headers={'Content-Type': 'application/json'},
                        data=json.dumps(scan_data))

# 3. データフラッシュ
flush_data = {'id': sample_id, 'index': 0}
response = requests.post(f'{base_url}io_flush',
                        headers={'Content-Type': 'application/json'},
                        data=json.dumps(flush_data))
```

