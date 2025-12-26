"""
RIR Study - NHANES Data Download Script
Downloads all required XPT files for the Residual Inflammatory Risk study

Available cycles with hs-CRP:
- 2005-2006, 2007-2008, 2009-2010 (CRP_D, CRP_E, CRP_F)
- 2015-2016, 2017-2020 (HSCRP_I, P_HSCRP)
NOTE: 2011-2014 does NOT have CRP data - excluded from analysis
"""

import os
import requests
from pathlib import Path
import urllib3

# Disable SSL warnings
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

# Project paths
PROJECT_ROOT = Path(__file__).parent.parent
DATA_RAW = PROJECT_ROOT / "data" / "raw"

# Headers
HEADERS = {
    'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36',
    'Accept': '*/*',
}

# Base URL
BASE_URL = "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public"

# Drug database for statin identification
DRUG_DB_URL = "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/1988/DataFiles/RXQ_DRUG.xpt"

# =============================================================================
# FILE DEFINITIONS BY CYCLE
# =============================================================================

FILES_2005_2006 = {
    # Demographics
    "DEMO_D.xpt": f"{BASE_URL}/2005/DataFiles/DEMO_D.xpt",
    
    # CRP
    "CRP_D.xpt": f"{BASE_URL}/2005/DataFiles/CRP_D.xpt",
    
    # Lipids
    "TCHOL_D.xpt": f"{BASE_URL}/2005/DataFiles/TCHOL_D.xpt",
    "HDL_D.xpt": f"{BASE_URL}/2005/DataFiles/HDL_D.xpt",
    "TRIGLY_D.xpt": f"{BASE_URL}/2005/DataFiles/TRIGLY_D.xpt",
    
    # Prescription medications
    "RXQ_RX_D.xpt": f"{BASE_URL}/2005/DataFiles/RXQ_RX_D.xpt",
    
    # Fasting subsample weights
    "FASTQX_D.xpt": f"{BASE_URL}/2005/DataFiles/FASTQX_D.xpt",
    
    # Covariates
    "BMX_D.xpt": f"{BASE_URL}/2005/DataFiles/BMX_D.xpt",
    "BPX_D.xpt": f"{BASE_URL}/2005/DataFiles/BPX_D.xpt",
    "BPQ_D.xpt": f"{BASE_URL}/2005/DataFiles/BPQ_D.xpt",
    "DIQ_D.xpt": f"{BASE_URL}/2005/DataFiles/DIQ_D.xpt",
    "GHB_D.xpt": f"{BASE_URL}/2005/DataFiles/GHB_D.xpt",
    "GLU_D.xpt": f"{BASE_URL}/2005/DataFiles/GLU_D.xpt",
    "SMQ_D.xpt": f"{BASE_URL}/2005/DataFiles/SMQ_D.xpt",
    "MCQ_D.xpt": f"{BASE_URL}/2005/DataFiles/MCQ_D.xpt",
}

FILES_2007_2008 = {
    "DEMO_E.xpt": f"{BASE_URL}/2007/DataFiles/DEMO_E.xpt",
    "CRP_E.xpt": f"{BASE_URL}/2007/DataFiles/CRP_E.xpt",
    "TCHOL_E.xpt": f"{BASE_URL}/2007/DataFiles/TCHOL_E.xpt",
    "HDL_E.xpt": f"{BASE_URL}/2007/DataFiles/HDL_E.xpt",
    "TRIGLY_E.xpt": f"{BASE_URL}/2007/DataFiles/TRIGLY_E.xpt",
    "RXQ_RX_E.xpt": f"{BASE_URL}/2007/DataFiles/RXQ_RX_E.xpt",
    "FASTQX_E.xpt": f"{BASE_URL}/2007/DataFiles/FASTQX_E.xpt",
    "BMX_E.xpt": f"{BASE_URL}/2007/DataFiles/BMX_E.xpt",
    "BPX_E.xpt": f"{BASE_URL}/2007/DataFiles/BPX_E.xpt",
    "BPQ_E.xpt": f"{BASE_URL}/2007/DataFiles/BPQ_E.xpt",
    "DIQ_E.xpt": f"{BASE_URL}/2007/DataFiles/DIQ_E.xpt",
    "GHB_E.xpt": f"{BASE_URL}/2007/DataFiles/GHB_E.xpt",
    "GLU_E.xpt": f"{BASE_URL}/2007/DataFiles/GLU_E.xpt",
    "SMQ_E.xpt": f"{BASE_URL}/2007/DataFiles/SMQ_E.xpt",
    "MCQ_E.xpt": f"{BASE_URL}/2007/DataFiles/MCQ_E.xpt",
}

FILES_2009_2010 = {
    "DEMO_F.xpt": f"{BASE_URL}/2009/DataFiles/DEMO_F.xpt",
    "CRP_F.xpt": f"{BASE_URL}/2009/DataFiles/CRP_F.xpt",
    "TCHOL_F.xpt": f"{BASE_URL}/2009/DataFiles/TCHOL_F.xpt",
    "HDL_F.xpt": f"{BASE_URL}/2009/DataFiles/HDL_F.xpt",
    "TRIGLY_F.xpt": f"{BASE_URL}/2009/DataFiles/TRIGLY_F.xpt",
    "RXQ_RX_F.xpt": f"{BASE_URL}/2009/DataFiles/RXQ_RX_F.xpt",
    "FASTQX_F.xpt": f"{BASE_URL}/2009/DataFiles/FASTQX_F.xpt",
    "BMX_F.xpt": f"{BASE_URL}/2009/DataFiles/BMX_F.xpt",
    "BPX_F.xpt": f"{BASE_URL}/2009/DataFiles/BPX_F.xpt",
    "BPQ_F.xpt": f"{BASE_URL}/2009/DataFiles/BPQ_F.xpt",
    "DIQ_F.xpt": f"{BASE_URL}/2009/DataFiles/DIQ_F.xpt",
    "GHB_F.xpt": f"{BASE_URL}/2009/DataFiles/GHB_F.xpt",
    "GLU_F.xpt": f"{BASE_URL}/2009/DataFiles/GLU_F.xpt",
    "SMQ_F.xpt": f"{BASE_URL}/2009/DataFiles/SMQ_F.xpt",
    "MCQ_F.xpt": f"{BASE_URL}/2009/DataFiles/MCQ_F.xpt",
}

# Skip 2011-2014: NO CRP DATA

FILES_2015_2016 = {
    "DEMO_I.xpt": f"{BASE_URL}/2015/DataFiles/DEMO_I.xpt",
    "HSCRP_I.xpt": f"{BASE_URL}/2015/DataFiles/HSCRP_I.xpt",
    "TCHOL_I.xpt": f"{BASE_URL}/2015/DataFiles/TCHOL_I.xpt",
    "HDL_I.xpt": f"{BASE_URL}/2015/DataFiles/HDL_I.xpt",
    "TRIGLY_I.xpt": f"{BASE_URL}/2015/DataFiles/TRIGLY_I.xpt",
    "RXQ_RX_I.xpt": f"{BASE_URL}/2015/DataFiles/RXQ_RX_I.xpt",
    "FASTQX_I.xpt": f"{BASE_URL}/2015/DataFiles/FASTQX_I.xpt",
    "BMX_I.xpt": f"{BASE_URL}/2015/DataFiles/BMX_I.xpt",
    "BPX_I.xpt": f"{BASE_URL}/2015/DataFiles/BPX_I.xpt",
    "BPQ_I.xpt": f"{BASE_URL}/2015/DataFiles/BPQ_I.xpt",
    "DIQ_I.xpt": f"{BASE_URL}/2015/DataFiles/DIQ_I.xpt",
    "GHB_I.xpt": f"{BASE_URL}/2015/DataFiles/GHB_I.xpt",
    "GLU_I.xpt": f"{BASE_URL}/2015/DataFiles/GLU_I.xpt",
    "SMQ_I.xpt": f"{BASE_URL}/2015/DataFiles/SMQ_I.xpt",
    "MCQ_I.xpt": f"{BASE_URL}/2015/DataFiles/MCQ_I.xpt",
}

FILES_2017_2020 = {
    "P_DEMO.xpt": f"{BASE_URL}/2017/DataFiles/P_DEMO.xpt",
    "P_HSCRP.xpt": f"{BASE_URL}/2017/DataFiles/P_HSCRP.xpt",
    "P_TCHOL.xpt": f"{BASE_URL}/2017/DataFiles/P_TCHOL.xpt",
    "P_HDL.xpt": f"{BASE_URL}/2017/DataFiles/P_HDL.xpt",
    "P_TRIGLY.xpt": f"{BASE_URL}/2017/DataFiles/P_TRIGLY.xpt",
    "RXQ_RX_J.xpt": f"{BASE_URL}/2017/DataFiles/RXQ_RX_J.xpt",
    "P_FASTQX.xpt": f"{BASE_URL}/2017/DataFiles/P_FASTQX.xpt",
    "P_BMX.xpt": f"{BASE_URL}/2017/DataFiles/P_BMX.xpt",
    "P_BPXO.xpt": f"{BASE_URL}/2017/DataFiles/P_BPXO.xpt",
    "P_BPQ.xpt": f"{BASE_URL}/2017/DataFiles/P_BPQ.xpt",
    "P_DIQ.xpt": f"{BASE_URL}/2017/DataFiles/P_DIQ.xpt",
    "P_GHB.xpt": f"{BASE_URL}/2017/DataFiles/P_GHB.xpt",
    "P_GLU.xpt": f"{BASE_URL}/2017/DataFiles/P_GLU.xpt",
    "P_SMQ.xpt": f"{BASE_URL}/2017/DataFiles/P_SMQ.xpt",
    "P_MCQ.xpt": f"{BASE_URL}/2017/DataFiles/P_MCQ.xpt",
}


def download_file(url: str, dest_path: Path, min_size: int = 1000) -> bool:
    """Download a file with proper headers."""
    try:
        if dest_path.exists() and dest_path.stat().st_size > min_size:
            print(f"  [EXISTS] {dest_path.name}")
            return True
            
        print(f"  [DOWNLOADING] {dest_path.name}...")
        response = requests.get(url, headers=HEADERS, verify=False, timeout=120)
        
        if response.status_code == 200 and len(response.content) > min_size:
            with open(dest_path, 'wb') as f:
                f.write(response.content)
            print(f"  [SUCCESS] {dest_path.name} ({len(response.content):,} bytes)")
            return True
        else:
            print(f"  [FAILED] {dest_path.name} - Status {response.status_code}")
            return False
            
    except Exception as e:
        print(f"  [ERROR] {dest_path.name}: {e}")
        return False


def download_cycle(cycle_name: str, files: dict):
    """Download all files for a cycle."""
    print(f"\n{'='*60}")
    print(f"CYCLE: {cycle_name}")
    print('='*60)
    
    cycle_dir = DATA_RAW / cycle_name
    cycle_dir.mkdir(parents=True, exist_ok=True)
    
    success_count = 0
    for filename, url in files.items():
        dest_path = cycle_dir / filename
        if download_file(url, dest_path):
            success_count += 1
    
    print(f"\n  Downloaded {success_count}/{len(files)} files")


def download_drug_database():
    """Download the NHANES drug database for statin identification."""
    print(f"\n{'='*60}")
    print("DRUG DATABASE (for statin identification)")
    print('='*60)
    
    db_dir = DATA_RAW / "drug_database"
    db_dir.mkdir(parents=True, exist_ok=True)
    
    dest_path = db_dir / "RXQ_DRUG.xpt"
    download_file(DRUG_DB_URL, dest_path)


def main():
    print("\n" + "="*60)
    print("RIR STUDY - NHANES DATA DOWNLOAD")
    print("Cycles: 2005-2010, 2015-2020 (skipping 2011-2014: no CRP)")
    print("="*60)
    
    DATA_RAW.mkdir(parents=True, exist_ok=True)
    
    # Download each cycle
    download_cycle("2005-2006", FILES_2005_2006)
    download_cycle("2007-2008", FILES_2007_2008)
    download_cycle("2009-2010", FILES_2009_2010)
    # Note: 2011-2014 skipped - no CRP data
    download_cycle("2015-2016", FILES_2015_2016)
    download_cycle("2017-2020", FILES_2017_2020)
    
    # Download drug database
    download_drug_database()
    
    print("\n" + "="*60)
    print("DOWNLOAD COMPLETE")
    print("="*60)
    print(f"\nData saved to: {DATA_RAW}")
    print("\nNote: 2011-2014 cycles excluded due to missing hs-CRP data")


if __name__ == "__main__":
    main()
