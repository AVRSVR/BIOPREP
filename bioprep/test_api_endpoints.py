import requests
import os
import json

BASE_URL = "http://127.0.0.1:5000"

def test_analyze():
    print("Testing /api/analyze...")
    pdb_path = r"c:\Users\anant\Documents\AG PLUGINS\protein prep\bioprep\test_1crn.pdb"
    with open(pdb_path, 'rb') as f:
        files = {'file': f}
        response = requests.post(f"{BASE_URL}/api/analyze", files=files)
    
    if response.status_code == 200:
        data = response.json()
        print("✅ Analyze success")
        return data.get('session_id')
    else:
        print(f"❌ Analyze failed: {response.text}")
        return None

def test_site_analyzer(session_id):
    if not session_id:
        print("Skipping Site Analyzer test: No session_id")
        return

    print(f"Testing /api/analyze-site/{session_id}...")
    response = requests.post(f"{BASE_URL}/api/analyze-site/{session_id}")
    
    if response.status_code == 200:
        data = response.json()
        print(f"✅ Site Analyzer success: Found {len(data.get('sites', []))} sites")
    else:
        print(f"❌ Site Analyzer failed: {response.text}")

if __name__ == "__main__":
    try:
        session_id = test_analyze()
        test_site_analyzer(session_id)
    except Exception as e:
        print(f"Test crashed: {e}")
