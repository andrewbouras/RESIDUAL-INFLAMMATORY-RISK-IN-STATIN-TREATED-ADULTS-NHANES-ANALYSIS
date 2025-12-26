"""
RIR Study - Cohort Construction
Builds the analytic cohort for Residual Inflammatory Risk analysis

Cohort definition:
- Adults ≥20 years
- Fasting subsample with lipid measurements
- Current statin use (from RXQ_RX prescription files)
- hs-CRP measured
- Exclude hs-CRP >10 mg/L

Primary analytic cohort: Statin users with LDL-C <70 mg/dL

RIR Definition: statin use + LDL-C <70 + hs-CRP ≥2 mg/L
"""

import pandas as pd
import numpy as np
from pathlib import Path
import warnings
import re

warnings.filterwarnings('ignore')

# Project paths
PROJECT_ROOT = Path(__file__).parent.parent
DATA_RAW = PROJECT_ROOT / "data" / "raw"
DATA_PROCESSED = PROJECT_ROOT / "data" / "processed"

# =============================================================================
# STATIN DRUG IDENTIFICATION
# =============================================================================

# Known statin generic names for pattern matching
STATIN_PATTERNS = [
    r'atorvastatin',
    r'fluvastatin', 
    r'lovastatin',
    r'pitavastatin',
    r'pravastatin',
    r'rosuvastatin',
    r'simvastatin',
    # Combination products
    r'vytorin',        # ezetimibe + simvastatin
    r'advicor',        # niacin + lovastatin
    r'caduet',         # amlodipine + atorvastatin
    r'liptruzet',      # ezetimibe + atorvastatin
]

# Therapeutic class codes for statins (if available in drug database)
STATIN_CLASS_CODES = ['CV350', 'CV351', 'CV352', 'CV359']  # HMG-CoA reductase inhibitors


def load_xpt(filepath: Path) -> pd.DataFrame:
    """Load XPT file into pandas DataFrame."""
    try:
        return pd.read_sas(filepath, format='xport', encoding='latin1')
    except Exception as e:
        print(f"  Warning: Could not load {filepath.name}: {e}")
        return pd.DataFrame()


def find_file(cycle_dir: Path, pattern: str) -> Path:
    """Find a file matching pattern in cycle directory."""
    matches = list(cycle_dir.glob(f"*{pattern}*"))
    if matches:
        return matches[0]
    return None


def identify_statin_users(rx_df: pd.DataFrame, drug_db: pd.DataFrame = None) -> pd.Series:
    """
    Identify statin users from prescription medication file.
    
    RXQ_RX files contain:
    - SEQN: participant ID
    - RXDDRUG: drug name (generic)
    - RXDDRGID: drug ID linking to drug database
    - RXDDAYS: days taking drug (30-day recall)
    """
    if rx_df.empty:
        return pd.Series(dtype=int)
    
    # Get drug name column (varies by cycle)
    drug_col = None
    for col in ['RXDDRUG', 'RXD240B', 'RXDDRGID']:
        if col in rx_df.columns:
            drug_col = col
            break
    
    if drug_col is None:
        print("  Warning: No drug name column found")
        return pd.Series(dtype=int)
    
    # Convert drug names to lowercase for matching
    rx_df = rx_df.copy()
    rx_df['drug_lower'] = rx_df[drug_col].astype(str).str.lower()
    
    # Match statin patterns
    statin_pattern = '|'.join(STATIN_PATTERNS)
    rx_df['is_statin'] = rx_df['drug_lower'].str.contains(statin_pattern, regex=True, na=False)
    
    # Aggregate to participant level
    statin_users = rx_df.groupby('SEQN')['is_statin'].any().astype(int)
    
    return statin_users


def compute_ldl_friedewald(tchol: pd.Series, hdl: pd.Series, trigly: pd.Series) -> pd.Series:
    """
    Calculate LDL using Friedewald equation: LDL = TC - HDL - (TG/5)
    Only valid when TG < 400 mg/dL.
    """
    ldl = tchol - hdl - (trigly / 5)
    ldl[trigly >= 400] = np.nan
    ldl[ldl < 0] = np.nan  # Invalid negative values
    return ldl


def define_hypertension(df: pd.DataFrame) -> pd.Series:
    """Define hypertension: SBP ≥130 OR DBP ≥80 OR on BP medication."""
    sbp_cols = [c for c in df.columns if 'BPXSY' in c or 'BPXOSY' in c]
    dbp_cols = [c for c in df.columns if 'BPXDI' in c or 'BPXODI' in c]
    
    mean_sbp = df[sbp_cols].mean(axis=1) if sbp_cols else pd.Series([np.nan] * len(df))
    mean_dbp = df[dbp_cols].mean(axis=1) if dbp_cols else pd.Series([np.nan] * len(df))
    
    bp_med = df.get('BPQ040A', pd.Series([np.nan] * len(df)))
    
    htn = ((mean_sbp >= 130) | (mean_dbp >= 80) | (bp_med == 1)).astype(float)
    return htn


def define_diabetes(df: pd.DataFrame) -> pd.Series:
    """Define diabetes: HbA1c ≥6.5% OR FPG ≥126 OR told by doctor OR on meds."""
    hba1c = df.get('LBXGH', pd.Series([np.nan] * len(df)))
    glucose = df.get('LBXGLU', pd.Series([np.nan] * len(df)))
    told = df.get('DIQ010', pd.Series([np.nan] * len(df)))
    
    dm = ((hba1c >= 6.5) | (glucose >= 126) | (told == 1)).astype(float)
    return dm


def define_smoking_status(df: pd.DataFrame) -> pd.Series:
    """Define smoking: 0=Never, 1=Former, 2=Current."""
    smoke_100 = df.get('SMQ020', pd.Series([np.nan] * len(df)))
    smoke_now = df.get('SMQ040', pd.Series([np.nan] * len(df)))
    
    status = pd.Series([np.nan] * len(df), index=df.index)
    status[smoke_100 == 2] = 0  # Never
    status[(smoke_100 == 1) & smoke_now.isin([1, 2])] = 2  # Current
    status[(smoke_100 == 1) & (smoke_now == 3)] = 1  # Former
    
    return status


def define_cvd_history(df: pd.DataFrame) -> pd.Series:
    """Define self-reported CVD history."""
    chd = df.get('MCQ160C', pd.Series([2] * len(df)))
    mi = df.get('MCQ160E', pd.Series([2] * len(df)))
    stroke = df.get('MCQ160F', pd.Series([2] * len(df)))
    angina = df.get('MCQ160D', pd.Series([2] * len(df)))
    chf = df.get('MCQ160B', pd.Series([2] * len(df)))
    
    cvd = ((chd == 1) | (mi == 1) | (stroke == 1) | (angina == 1) | (chf == 1)).astype(float)
    return cvd


def process_cycle(cycle: str, cycle_dir: Path) -> pd.DataFrame:
    """Process a single NHANES cycle."""
    print(f"\n{'='*60}")
    print(f"Processing {cycle}")
    print('='*60)
    
    if not cycle_dir.exists():
        print(f"  Directory not found: {cycle_dir}")
        return pd.DataFrame()
    
    # Determine if pre-pandemic 2017-2020
    is_prepandemic = "2017" in cycle
    
    # Load demographics
    demo_file = find_file(cycle_dir, "DEMO")
    if not demo_file:
        print("  No DEMO file found")
        return pd.DataFrame()
    
    df = load_xpt(demo_file)
    print(f"  Loaded {len(df):,} participants from DEMO")
    
    # Merge component files
    components = [
        ("CRP", "CRP/hs-CRP"),
        ("HSCRP", "hs-CRP"),
        ("TCHOL", "Total cholesterol"),
        ("HDL", "HDL"),
        ("TRIGLY", "Triglycerides"),
        ("FASTQX", "Fasting weights"),
        ("BMX", "Body measures"),
        ("BPX", "Blood pressure"),
        ("BPXO", "Blood pressure (oscillometric)"),
        ("BPQ", "BP questionnaire"),
        ("DIQ", "Diabetes"),
        ("GHB", "Glycohemoglobin"),
        ("GLU", "Glucose"),
        ("SMQ", "Smoking"),
        ("MCQ", "Medical conditions"),
    ]
    
    for pattern, description in components:
        file_path = find_file(cycle_dir, pattern)
        if file_path:
            comp_df = load_xpt(file_path)
            if 'SEQN' in comp_df.columns and len(comp_df) > 0:
                df = df.merge(comp_df, on='SEQN', how='left', suffixes=('', f'_{pattern}'))
                print(f"  Merged {description}: {len(comp_df):,} records")
    
    # Load prescription data and identify statin users
    rx_file = find_file(cycle_dir, "RXQ_RX")
    if rx_file:
        rx_df = load_xpt(rx_file)
        print(f"  Loaded prescription data: {len(rx_df):,} medication records")
        statin_users = identify_statin_users(rx_df)
        df['statin_use'] = df['SEQN'].map(statin_users).fillna(0).astype(int)
        print(f"  Identified {df['statin_use'].sum():,} statin users")
    else:
        print("  Warning: No prescription file found")
        df['statin_use'] = 0
    
    # Add cycle info
    df['cycle'] = cycle
    df['is_prepandemic'] = is_prepandemic
    
    # Harmonize variable names across cycles
    # CRP: LBXCRP (2005-2010) is in mg/dL; LBXHSCRP (2015+) is in mg/L
    # Convert 2005-2010 CRP from mg/dL to mg/L by multiplying by 10
    if 'LBXCRP' in df.columns and 'LBXHSCRP' not in df.columns:
        df['LBXHSCRP'] = df['LBXCRP'] * 10  # Convert mg/dL to mg/L
        print(f"  Converted LBXCRP (mg/dL) to mg/L for {cycle}")
    
    # Fasting weights: WTSAF2YR (2005-2016) vs WTSAFPRP (2017-2020)
    if 'WTSAFPRP' in df.columns and 'WTSAF2YR' not in df.columns:
        df['WTSAF2YR'] = df['WTSAFPRP']
    
    # MEC weights: WTMEC2YR (2005-2016) vs WTMECPRP (2017-2020)
    if 'WTMECPRP' in df.columns and 'WTMEC2YR' not in df.columns:
        df['WTMEC2YR'] = df['WTMECPRP']
    
    # HDL: LBDHDD vs LBXHDD
    if 'LBXHDD' in df.columns and 'LBDHDD' not in df.columns:
        df['LBDHDD'] = df['LBXHDD']
    
    return df


def build_analytic_cohort(df: pd.DataFrame) -> pd.DataFrame:
    """Build the RIR analytic cohort with all derived variables."""
    print("\n" + "="*60)
    print("Building analytic cohort")
    print("="*60)
    
    # Demographics
    df['age'] = df['RIDAGEYR']
    df['sex'] = df['RIAGENDR']  # 1=Male, 2=Female
    df['female'] = (df['sex'] == 2).astype(int)
    
    # Race/ethnicity (use RIDRETH1 for 2005-2010, RIDRETH3 for 2015+)
    if 'RIDRETH3' in df.columns:
        df['race_eth'] = df['RIDRETH3']
    else:
        df['race_eth'] = df.get('RIDRETH1', pd.Series([np.nan] * len(df)))
    
    # hs-CRP (harmonized to LBXHSCRP in process_cycle)
    df['hscrp'] = df.get('LBXHSCRP', pd.Series([np.nan] * len(df), index=df.index))
    print(f"  hs-CRP available for {df['hscrp'].notna().sum():,} participants")
    
    # Lipids (LBDHDD harmonized in process_cycle)
    df['tchol'] = df.get('LBXTC', pd.Series([np.nan] * len(df), index=df.index))
    df['hdl'] = df.get('LBDHDD', pd.Series([np.nan] * len(df), index=df.index))
    df['trigly'] = df.get('LBXTR', pd.Series([np.nan] * len(df), index=df.index))
    
    # Calculate LDL (Friedewald)
    df['ldl'] = compute_ldl_friedewald(df['tchol'], df['hdl'], df['trigly'])
    print(f"  LDL-C available for {df['ldl'].notna().sum():,} participants")
    
    # Fasting status and weights (harmonized to WTSAF2YR in process_cycle)
    df['fasting_weight'] = df.get('WTSAF2YR', pd.Series([np.nan] * len(df), index=df.index))
    df['in_fasting_subsample'] = df['fasting_weight'].notna().astype(int)
    print(f"  Fasting subsample: {df['in_fasting_subsample'].sum():,} participants")
    
    # Survey design variables
    df['mec_weight'] = df.get('WTMEC2YR', df.get('WTMECPRP', pd.Series([np.nan] * len(df))))
    df['psu'] = df['SDMVPSU']
    df['strata'] = df['SDMVSTRA']
    
    # BMI
    df['bmi'] = df.get('BMXBMI', pd.Series([np.nan] * len(df)))
    df['bmi_cat'] = pd.cut(df['bmi'], 
                           bins=[0, 18.5, 25, 30, 100],
                           labels=['Underweight', 'Normal', 'Overweight', 'Obese'])
    df['obesity'] = (df['bmi'] >= 30).astype(int)
    
    # Clinical conditions
    df['hypertension'] = define_hypertension(df)
    df['diabetes'] = define_diabetes(df)
    df['hba1c'] = df.get('LBXGH', pd.Series([np.nan] * len(df)))
    df['fasting_glucose'] = df.get('LBXGLU', pd.Series([np.nan] * len(df)))
    
    # Smoking
    df['smoking_status'] = define_smoking_status(df)
    df['current_smoker'] = (df['smoking_status'] == 2).astype(int)
    
    # CVD history
    df['cvd_history'] = define_cvd_history(df)
    
    # ==========================================================================
    # ELIGIBILITY AND RIR DEFINITION
    # ==========================================================================
    
    print("\n" + "-"*40)
    print("Applying eligibility criteria")
    print("-"*40)
    
    n_total = len(df)
    print(f"  Total participants: {n_total:,}")
    
    # Age ≥20
    df['age_eligible'] = (df['age'] >= 20).astype(int)
    n_age = df['age_eligible'].sum()
    print(f"  Age ≥20: {n_age:,}")
    
    # Fasting subsample
    df['fasting_eligible'] = df['in_fasting_subsample']
    n_fasting = (df['age_eligible'] & df['fasting_eligible']).sum()
    print(f"  In fasting subsample: {n_fasting:,}")
    
    # Non-missing hs-CRP
    df['crp_available'] = df['hscrp'].notna().astype(int)
    n_crp = (df['age_eligible'] & df['fasting_eligible'] & df['crp_available']).sum()
    print(f"  hs-CRP available: {n_crp:,}")
    
    # Non-missing LDL
    df['ldl_available'] = df['ldl'].notna().astype(int)
    n_ldl = (df['age_eligible'] & df['fasting_eligible'] & df['crp_available'] & df['ldl_available']).sum()
    print(f"  LDL-C available: {n_ldl:,}")
    
    # Exclude hs-CRP >10 (acute inflammation)
    df['crp_valid'] = (df['hscrp'] <= 10).astype(int)
    
    # Overall eligibility
    df['eligible'] = (
        df['age_eligible'] & 
        df['fasting_eligible'] & 
        df['crp_available'] & 
        df['ldl_available'] &
        df['crp_valid']
    ).astype(int)
    n_eligible = df['eligible'].sum()
    print(f"\n  Overall eligible: {n_eligible:,}")
    
    # Statin users in eligible population
    df['statin_eligible'] = (df['eligible'] & df['statin_use']).astype(int)
    n_statin = df['statin_eligible'].sum()
    print(f"  Statin users (eligible): {n_statin:,}")
    
    # ==========================================================================
    # RIR DEFINITION
    # ==========================================================================
    
    print("\n" + "-"*40)
    print("Defining RIR")
    print("-"*40)
    
    # Primary analytic cohort: statin users with LDL <70
    df['ldl_lt70'] = (df['ldl'] < 70).astype(int)
    df['primary_cohort'] = (df['statin_eligible'] & df['ldl_lt70']).astype(int)
    n_primary = df['primary_cohort'].sum()
    print(f"  Statin + LDL <70 (primary cohort): {n_primary:,}")
    
    # RIR definition: hs-CRP ≥2 mg/L
    df['crp_elevated'] = (df['hscrp'] >= 2).astype(int)
    df['rir'] = (df['primary_cohort'] & df['crp_elevated']).astype(int)
    n_rir = df['rir'].sum()
    print(f"  RIR (hs-CRP ≥2): {n_rir:,}")
    
    if n_primary > 0:
        rir_pct = 100 * n_rir / n_primary
        print(f"  RIR prevalence in primary cohort: {rir_pct:.1f}%")
    
    # Sensitivity definitions
    df['rir_crp3'] = (df['primary_cohort'] & (df['hscrp'] >= 3)).astype(int)
    
    df['ldl_lt55'] = (df['ldl'] < 55).astype(int)
    df['cohort_ldl55'] = (df['statin_eligible'] & df['ldl_lt55']).astype(int)
    df['rir_ldl55'] = (df['cohort_ldl55'] & df['crp_elevated']).astype(int)
    
    print(f"  Sensitivity: RIR (hs-CRP ≥3): {df['rir_crp3'].sum():,}")
    print(f"  Sensitivity: LDL <55 cohort: {df['cohort_ldl55'].sum():,}")
    
    return df


def main():
    """Main cohort construction function."""
    print("\n" + "="*60)
    print("RIR STUDY - COHORT CONSTRUCTION")
    print("="*60)
    
    # Process each cycle
    cycles = [
        ("2005-2006", DATA_RAW / "2005-2006"),
        ("2007-2008", DATA_RAW / "2007-2008"),
        ("2009-2010", DATA_RAW / "2009-2010"),
        # 2011-2014 skipped - no CRP
        ("2015-2016", DATA_RAW / "2015-2016"),
        ("2017-2020", DATA_RAW / "2017-2020"),
    ]
    
    all_data = []
    for cycle_name, cycle_dir in cycles:
        df = process_cycle(cycle_name, cycle_dir)
        if len(df) > 0:
            all_data.append(df)
    
    if not all_data:
        print("\nNo data loaded. Run 01_download_data.py first.")
        return
    
    # Combine cycles
    combined = pd.concat(all_data, ignore_index=True)
    print(f"\n{'='*60}")
    print(f"Combined: {len(combined):,} participants across all cycles")
    print('='*60)
    
    # Build analytic cohort
    combined = build_analytic_cohort(combined)
    
    # Select key variables for export
    export_vars = [
        # Identifiers
        'SEQN', 'cycle', 'is_prepandemic',
        
        # Survey design
        'mec_weight', 'fasting_weight', 'psu', 'strata',
        
        # Demographics
        'age', 'sex', 'female', 'race_eth',
        
        # Exposures
        'statin_use', 'ldl', 'ldl_lt70', 'ldl_lt55',
        'hscrp', 'crp_elevated',
        
        # Lipids
        'tchol', 'hdl', 'trigly',
        
        # Clinical
        'bmi', 'bmi_cat', 'obesity',
        'hypertension', 'diabetes', 'hba1c', 'fasting_glucose',
        'smoking_status', 'current_smoker',
        'cvd_history',
        
        # Eligibility flags
        'eligible', 'in_fasting_subsample', 'crp_valid',
        'statin_eligible', 'primary_cohort',
        
        # Outcomes
        'rir', 'rir_crp3', 'cohort_ldl55', 'rir_ldl55',
    ]
    
    export_vars = [v for v in export_vars if v in combined.columns]
    
    # Save
    DATA_PROCESSED.mkdir(parents=True, exist_ok=True)
    
    output_path = DATA_PROCESSED / "rir_cohort.parquet"
    combined[export_vars].to_parquet(output_path)
    print(f"\nSaved: {output_path}")
    
    csv_path = DATA_PROCESSED / "rir_cohort.csv"
    combined[export_vars].to_csv(csv_path, index=False)
    print(f"Saved: {csv_path}")
    
    # Summary
    print("\n" + "="*60)
    print("COHORT SUMMARY")
    print("="*60)
    
    print(f"\nTotal participants: {len(combined):,}")
    print(f"Eligible for analysis: {combined['eligible'].sum():,}")
    print(f"Statin users (eligible): {combined['statin_eligible'].sum():,}")
    print(f"Primary cohort (statin + LDL<70): {combined['primary_cohort'].sum():,}")
    print(f"RIR (hs-CRP ≥2): {combined['rir'].sum():,}")
    
    print("\nBy cycle:")
    cycle_summary = combined.groupby('cycle').agg({
        'eligible': 'sum',
        'statin_eligible': 'sum',
        'primary_cohort': 'sum',
        'rir': 'sum'
    })
    print(cycle_summary)


if __name__ == "__main__":
    main()

