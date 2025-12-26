"""
RIR Study - Data Processing and Cohort Creation
1. Load and merge NHANES files (2005-2020)
2. Identify statin users from prescription medication files
3. Calculate LDL-C (Friedewald, with direct LDL fallback)
4. Create analytic cohort with all variables
5. Apply inclusion/exclusion criteria
6. Export clean dataset for R analysis
"""

import pandas as pd
import numpy as np
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Project paths
PROJECT_ROOT = Path(__file__).parent.parent
DATA_RAW = PROJECT_ROOT / "data" / "raw"
DATA_PROCESSED = PROJECT_ROOT / "data" / "processed"
OUTPUT_DIR = PROJECT_ROOT / "output"

# NHANES cycles with CRP data available
# Note: 2011-2014 HSCRP files are not publicly available via direct download
CYCLES = ["2005-2006", "2007-2008", "2009-2010", 
          "2015-2016", "2017-2020"]

# ============================================================================
# STATIN IDENTIFICATION
# ============================================================================

# Statin drug names (generic names used in NHANES RXQ_RX)
STATIN_NAMES = [
    'atorvastatin',
    'simvastatin', 
    'rosuvastatin',
    'pravastatin',
    'lovastatin',
    'fluvastatin',
    'pitavastatin',
    # Combination products
    'atorvastatin/amlodipine',
    'simvastatin/ezetimibe',
    'simvastatin/niacin',
    'atorvastatin/ezetimibe',
    'rosuvastatin/ezetimibe',
    'lovastatin/niacin',
    # Brand names that may appear
    'lipitor',
    'zocor',
    'crestor',
    'pravachol',
    'mevacor',
    'lescol',
    'livalo',
    'vytorin',
    'caduet',
    'advicor',
]

# Therapeutic class codes for statins (HMG-CoA reductase inhibitors)
# RXDDCI1A/B/C/D contain drug class codes
STATIN_THERAPEUTIC_CLASSES = [
    '358',  # HMG-CoA reductase inhibitors
]


def load_xpt(filepath: Path) -> pd.DataFrame:
    """Load XPT file into pandas DataFrame."""
    try:
        df = pd.read_sas(filepath, format='xport', encoding='latin1')
        return df
    except Exception as e:
        print(f"  Warning: Could not load {filepath.name}: {e}")
        return pd.DataFrame()


def find_file(cycle_dir: Path, pattern: str) -> Path | None:
    """Find a file matching pattern in cycle directory (case-insensitive)."""
    for ext in ['.XPT', '.xpt']:
        matches = list(cycle_dir.glob(f"*{pattern}*{ext}"))
        if matches:
            return matches[0]
    return None


def identify_statin_users(rxq_df: pd.DataFrame) -> pd.DataFrame:
    """
    Identify statin users from RXQ_RX prescription medication data.
    Returns DataFrame with SEQN and statin_user flag.
    """
    if rxq_df.empty or 'SEQN' not in rxq_df.columns:
        return pd.DataFrame(columns=['SEQN', 'statin_user'])
    
    # Standardize column names to uppercase
    rxq_df.columns = [c.upper() for c in rxq_df.columns]
    
    statin_users = set()
    
    # Method 1: Check drug name (RXDDRUG)
    drug_name_cols = ['RXDDRUG', 'RXDRUG']
    for col in drug_name_cols:
        if col in rxq_df.columns:
            for idx, row in rxq_df.iterrows():
                drug_name = str(row.get(col, '')).lower()
                if any(statin in drug_name for statin in STATIN_NAMES):
                    statin_users.add(row['SEQN'])
    
    # Method 2: Check therapeutic class codes
    class_cols = ['RXDDCI1A', 'RXDDCI1B', 'RXDDCI1C', 'RXDDCI2A', 'RXDDCI2B', 'RXDDCI2C']
    for col in class_cols:
        if col in rxq_df.columns:
            for idx, row in rxq_df.iterrows():
                class_code = str(row.get(col, ''))
                if class_code in STATIN_THERAPEUTIC_CLASSES:
                    statin_users.add(row['SEQN'])
    
    # Method 3: Check generic drug code (RXDDCN1A, RXDDCN1B for drug category)
    # Category 019 = cardiovascular agents, but too broad - rely on names/classes
    
    # Create output dataframe with all unique SEQNs
    all_seqn = rxq_df['SEQN'].unique()
    result = pd.DataFrame({'SEQN': all_seqn})
    result['statin_user'] = result['SEQN'].isin(statin_users).astype(int)
    
    return result


def calculate_ldl(tc: pd.Series, hdl: pd.Series, tg: pd.Series, 
                  direct_ldl: pd.Series = None) -> pd.Series:
    """
    Calculate LDL using Friedewald equation: LDL = TC - HDL - (TG/5)
    Uses direct LDL as fallback when Friedewald is invalid or missing.
    Friedewald only valid when TG < 400 mg/dL.
    """
    # Friedewald calculation
    ldl_calc = tc - hdl - (tg / 5)
    
    # Set to NaN if TG >= 400 (Friedewald not valid)
    ldl_calc = ldl_calc.where(tg < 400, np.nan)
    
    # Use direct LDL as fallback if available
    if direct_ldl is not None:
        ldl_calc = ldl_calc.fillna(direct_ldl)
    
    return ldl_calc


def process_cycle(cycle: str) -> pd.DataFrame:
    """Process a single NHANES cycle and return merged data."""
    print(f"\n[Processing {cycle}]")
    
    cycle_dir = DATA_RAW / cycle
    if not cycle_dir.exists():
        print(f"  Directory not found: {cycle_dir}")
        return pd.DataFrame()
    
    # -------------------------------------------------------------------------
    # Load Demographics
    # -------------------------------------------------------------------------
    demo_file = find_file(cycle_dir, "DEMO")
    if not demo_file:
        print(f"  No DEMO file found")
        return pd.DataFrame()
    
    df = load_xpt(demo_file)
    if df.empty:
        return pd.DataFrame()
    
    df.columns = [c.upper() for c in df.columns]
    print(f"  Loaded {len(df):,} participants from {demo_file.name}")
    
    # Keep relevant demographic columns
    demo_cols = ['SEQN', 'RIDAGEYR', 'RIAGENDR', 'RIDRETH1', 'RIDRETH3',
                 'DMDEDUC2', 'INDFMPIR', 'WTMEC2YR', 'WTINT2YR', 
                 'SDMVPSU', 'SDMVSTRA']
    
    # Handle 2017-2020 pre-pandemic weight naming
    if 'WTMECPRP' in df.columns:
        df['WTMEC2YR'] = df['WTMECPRP']
    if 'WTINTPRP' in df.columns:
        df['WTINT2YR'] = df['WTINTPRP']
    
    df = df[[c for c in demo_cols if c in df.columns]].copy()
    
    # -------------------------------------------------------------------------
    # Load and merge laboratory data
    # -------------------------------------------------------------------------
    lab_files = {
        'TCHOL': ['LBXTC'],                          # Total cholesterol
        'HDL': ['LBDHDD', 'LBXHDD'],                 # HDL cholesterol
        'TRIGLY': ['LBXTR', 'WTSAF2YR', 'WTSAFPRP'], # Triglycerides + fasting weight
        'GHB': ['LBXGH'],                             # HbA1c
        'GLU': ['LBXGLU'],                            # Fasting glucose
    }
    
    for file_base, var_names in lab_files.items():
        lab_file = find_file(cycle_dir, file_base)
        if lab_file:
            lab_df = load_xpt(lab_file)
            if not lab_df.empty:
                lab_df.columns = [c.upper() for c in lab_df.columns]
                # Keep only SEQN and target variables
                keep_cols = ['SEQN'] + [v for v in var_names if v in lab_df.columns]
                if len(keep_cols) > 1:
                    df = df.merge(lab_df[keep_cols], on='SEQN', how='left')
    
    # Load CRP/HSCRP separately to handle different file naming
    # Try HSCRP first (2011+), then CRP (2005-2010)
    crp_loaded = False
    for pattern in ['HSCRP', 'CRP']:
        crp_file = find_file(cycle_dir, pattern)
        if crp_file and not crp_loaded:
            crp_df = load_xpt(crp_file)
            if not crp_df.empty:
                crp_df.columns = [c.upper() for c in crp_df.columns]
                # CRP variable may be named LBXCRP or LBXHSCRP
                crp_vars = [v for v in ['LBXHSCRP', 'LBXCRP'] if v in crp_df.columns]
                if crp_vars:
                    keep_cols = ['SEQN'] + crp_vars
                    df = df.merge(crp_df[keep_cols], on='SEQN', how='left')
                    crp_loaded = True
    
    # -------------------------------------------------------------------------
    # Load and merge examination data (BMI, Blood Pressure)
    # -------------------------------------------------------------------------
    # Body measures
    bmx_file = find_file(cycle_dir, "BMX")
    if bmx_file:
        bmx_df = load_xpt(bmx_file)
        if not bmx_df.empty:
            bmx_df.columns = [c.upper() for c in bmx_df.columns]
            bmx_cols = ['SEQN', 'BMXBMI', 'BMXWT', 'BMXHT']
            keep = [c for c in bmx_cols if c in bmx_df.columns]
            df = df.merge(bmx_df[keep], on='SEQN', how='left')
    
    # Blood pressure
    bp_file = find_file(cycle_dir, "BPX")
    if not bp_file:
        bp_file = find_file(cycle_dir, "BPXO")  # 2017-2020 uses BPXO
    if bp_file:
        bp_df = load_xpt(bp_file)
        if not bp_df.empty:
            bp_df.columns = [c.upper() for c in bp_df.columns]
            bp_cols = ['SEQN', 'BPXSY1', 'BPXSY2', 'BPXSY3', 
                      'BPXDI1', 'BPXDI2', 'BPXDI3',
                      'BPXOSY1', 'BPXOSY2', 'BPXOSY3',  # 2017-2020 oscillometric
                      'BPXODI1', 'BPXODI2', 'BPXODI3']
            keep = [c for c in bp_cols if c in bp_df.columns]
            df = df.merge(bp_df[keep], on='SEQN', how='left')
    
    # -------------------------------------------------------------------------
    # Load and merge questionnaire data
    # -------------------------------------------------------------------------
    # Diabetes questionnaire
    diq_file = find_file(cycle_dir, "DIQ")
    if diq_file:
        diq_df = load_xpt(diq_file)
        if not diq_df.empty:
            diq_df.columns = [c.upper() for c in diq_df.columns]
            diq_cols = ['SEQN', 'DIQ010', 'DIQ050', 'DIQ070']
            keep = [c for c in diq_cols if c in diq_df.columns]
            df = df.merge(diq_df[keep], on='SEQN', how='left')
    
    # Blood pressure questionnaire
    bpq_file = find_file(cycle_dir, "BPQ")
    if bpq_file:
        bpq_df = load_xpt(bpq_file)
        if not bpq_df.empty:
            bpq_df.columns = [c.upper() for c in bpq_df.columns]
            bpq_cols = ['SEQN', 'BPQ020', 'BPQ040A']  # Ever told high BP, taking meds
            keep = [c for c in bpq_cols if c in bpq_df.columns]
            df = df.merge(bpq_df[keep], on='SEQN', how='left')
    
    # Smoking questionnaire
    smq_file = find_file(cycle_dir, "SMQ")
    if smq_file:
        smq_df = load_xpt(smq_file)
        if not smq_df.empty:
            smq_df.columns = [c.upper() for c in smq_df.columns]
            smq_cols = ['SEQN', 'SMQ020', 'SMQ040']
            keep = [c for c in smq_cols if c in smq_df.columns]
            df = df.merge(smq_df[keep], on='SEQN', how='left')
    
    # -------------------------------------------------------------------------
    # Load prescription medications and identify statin users
    # -------------------------------------------------------------------------
    rxq_file = find_file(cycle_dir, "RXQ_RX")
    if rxq_file:
        rxq_df = load_xpt(rxq_file)
        if not rxq_df.empty:
            statin_df = identify_statin_users(rxq_df)
            if not statin_df.empty:
                df = df.merge(statin_df, on='SEQN', how='left')
                df['statin_user'] = df['statin_user'].fillna(0).astype(int)
                print(f"  Identified {statin_df['statin_user'].sum():,} statin users")
    else:
        df['statin_user'] = 0
        print(f"  Warning: No RXQ_RX file found for {cycle}")
    
    # Add cycle identifier
    df['cycle'] = cycle
    
    return df


def harmonize_variables(df: pd.DataFrame) -> pd.DataFrame:
    """Harmonize variable names and create derived variables."""
    print("\n[Harmonizing variables]")
    
    # -------------------------------------------------------------------------
    # Standardize CRP variable name
    # Coalesce LBXHSCRP (2011+) and LBXCRP (2005-2010) since both may exist
    # after concatenating cycles
    # -------------------------------------------------------------------------
    hscrp_col = df.get('LBXHSCRP', pd.Series([np.nan] * len(df), index=df.index))
    crp_col = df.get('LBXCRP', pd.Series([np.nan] * len(df), index=df.index))
    df['hscrp'] = hscrp_col.fillna(crp_col)
    
    # -------------------------------------------------------------------------
    # Standardize lipid variables
    # -------------------------------------------------------------------------
    df['total_chol'] = df.get('LBXTC', np.nan)
    
    # HDL - may be LBDHDD or LBXHDD depending on cycle
    if 'LBDHDD' in df.columns:
        df['hdl'] = df['LBDHDD']
    elif 'LBXHDD' in df.columns:
        df['hdl'] = df['LBXHDD']
    else:
        df['hdl'] = np.nan
    
    df['trigly'] = df.get('LBXTR', np.nan)
    
    # Calculate LDL (Friedewald)
    df['ldl_calc'] = calculate_ldl(df['total_chol'], df['hdl'], df['trigly'])
    
    # -------------------------------------------------------------------------
    # Standardize fasting subsample weight
    # Coalesce WTSAF2YR (most cycles) and WTSAFPRP (2017-2020)
    # -------------------------------------------------------------------------
    wtsaf2yr = df.get('WTSAF2YR', pd.Series([np.nan] * len(df), index=df.index))
    wtsafprp = df.get('WTSAFPRP', pd.Series([np.nan] * len(df), index=df.index))
    df['fasting_weight'] = wtsaf2yr.fillna(wtsafprp)
    
    # -------------------------------------------------------------------------
    # Standardize demographics
    # -------------------------------------------------------------------------
    df['age'] = df.get('RIDAGEYR', np.nan)
    df['sex'] = df.get('RIAGENDR', np.nan)  # 1=Male, 2=Female
    df['race_eth'] = df.get('RIDRETH1', np.nan)  # 1=Mex Am, 2=Other Hisp, 3=NH White, 4=NH Black, 5=Other
    df['education'] = df.get('DMDEDUC2', np.nan)
    df['pir'] = df.get('INDFMPIR', np.nan)  # Poverty income ratio
    
    # -------------------------------------------------------------------------
    # Standardize blood pressure
    # -------------------------------------------------------------------------
    # Average available SBP readings
    sbp_cols = [c for c in df.columns if c.startswith(('BPXSY', 'BPXOSY'))]
    dbp_cols = [c for c in df.columns if c.startswith(('BPXDI', 'BPXODI'))]
    
    if sbp_cols:
        df['sbp_mean'] = df[sbp_cols].mean(axis=1)
    else:
        df['sbp_mean'] = np.nan
        
    if dbp_cols:
        df['dbp_mean'] = df[dbp_cols].mean(axis=1)
    else:
        df['dbp_mean'] = np.nan
    
    # -------------------------------------------------------------------------
    # Standardize other variables
    # -------------------------------------------------------------------------
    df['bmi'] = df.get('BMXBMI', np.nan)
    df['hba1c'] = df.get('LBXGH', np.nan)
    df['fasting_glucose'] = df.get('LBXGLU', np.nan)
    
    # Survey design variables
    df['psu'] = df.get('SDMVPSU', np.nan)
    df['strata'] = df.get('SDMVSTRA', np.nan)
    df['mec_weight'] = df.get('WTMEC2YR', np.nan)
    
    # -------------------------------------------------------------------------
    # Create derived clinical variables
    # -------------------------------------------------------------------------
    # Diabetes: HbA1c ≥6.5% OR FPG ≥126 OR told by doctor OR on insulin/oral meds
    df['diabetes'] = (
        (df['hba1c'] >= 6.5) |
        (df['fasting_glucose'] >= 126) |
        (df.get('DIQ010', 2) == 1) |  # Doctor told diabetes
        (df.get('DIQ050', 2) == 1) |  # Taking insulin
        (df.get('DIQ070', 2) == 1)    # Taking oral diabetes meds
    ).astype(float)
    df.loc[df[['hba1c', 'fasting_glucose', 'DIQ010']].isna().all(axis=1), 'diabetes'] = np.nan
    
    # Prediabetes: HbA1c 5.7-6.4% OR FPG 100-125 mg/dL (excluding those with diabetes)
    df['prediabetes'] = (
        ((df['hba1c'] >= 5.7) & (df['hba1c'] < 6.5)) |
        ((df['fasting_glucose'] >= 100) & (df['fasting_glucose'] < 126))
    ).astype(float)
    # Exclude those already classified as diabetic
    df.loc[df['diabetes'] == 1, 'prediabetes'] = 0
    df.loc[df[['hba1c', 'fasting_glucose']].isna().all(axis=1), 'prediabetes'] = np.nan
    
    # Combined glycemic status: 0=Normal, 1=Prediabetes, 2=Diabetes
    df['glycemic_status'] = 0
    df.loc[df['prediabetes'] == 1, 'glycemic_status'] = 1
    df.loc[df['diabetes'] == 1, 'glycemic_status'] = 2
    df.loc[df[['diabetes', 'prediabetes']].isna().all(axis=1), 'glycemic_status'] = np.nan
    
    # Hypertension: SBP ≥130 OR DBP ≥80 OR taking BP medication
    df['hypertension'] = (
        (df['sbp_mean'] >= 130) |
        (df['dbp_mean'] >= 80) |
        (df.get('BPQ040A', 2) == 1)  # Taking BP medication
    ).astype(float)
    df.loc[df[['sbp_mean', 'dbp_mean', 'BPQ040A']].isna().all(axis=1), 'hypertension'] = np.nan
    
    # Smoking status: 1=Current, 2=Former, 3=Never
    smoke_100 = df.get('SMQ020', np.nan)  # Smoked 100+ cigarettes
    smoke_now = df.get('SMQ040', np.nan)  # Currently smoke
    
    df['smoking_status'] = np.nan
    df.loc[smoke_100 == 2, 'smoking_status'] = 3  # Never (didn't smoke 100+)
    df.loc[(smoke_100 == 1) & (smoke_now.isin([1, 2])), 'smoking_status'] = 1  # Current
    df.loc[(smoke_100 == 1) & (smoke_now == 3), 'smoking_status'] = 2  # Former
    
    df['current_smoker'] = (df['smoking_status'] == 1).astype(float)
    df.loc[df['smoking_status'].isna(), 'current_smoker'] = np.nan
    
    # Obesity
    df['obese'] = (df['bmi'] >= 30).astype(float)
    df.loc[df['bmi'].isna(), 'obese'] = np.nan
    
    return df


def create_rir_cohort(df: pd.DataFrame) -> pd.DataFrame:
    """
    Create the RIR study cohort with inclusion/exclusion criteria.
    
    Inclusion:
    - Adults ≥20 years
    - Statin user
    - Has fasting LDL-C measured
    - Has hs-CRP measured
    
    Exclusion:
    - hs-CRP >10 mg/L (acute inflammation)
    - Missing key variables
    
    Primary analytic cohort: Statin users with LDL-C <70 mg/dL
    """
    print("\n[Creating RIR cohort]")
    print(f"  Starting N: {len(df):,}")
    
    # Adults ≥20 years
    df = df[df['age'] >= 20].copy()
    print(f"  After age ≥20: {len(df):,}")
    
    # Has hs-CRP measured
    df = df[df['hscrp'].notna()].copy()
    print(f"  After requiring hs-CRP: {len(df):,}")
    
    # Exclude hs-CRP >10 (acute inflammation)
    df = df[df['hscrp'] <= 10].copy()
    print(f"  After excluding hs-CRP >10: {len(df):,}")
    
    # Has statin use information
    df = df[df['statin_user'].notna()].copy()
    print(f"  After requiring statin info: {len(df):,}")
    
    # Has LDL-C calculated (requires TC, HDL, TG with TG<400)
    df = df[df['ldl_calc'].notna()].copy()
    print(f"  After requiring LDL-C: {len(df):,}")
    
    # Has fasting subsample weight
    df = df[df['fasting_weight'].notna()].copy()
    print(f"  After requiring fasting weight: {len(df):,}")
    
    # Has survey design variables
    df = df[(df['psu'].notna()) & (df['strata'].notna())].copy()
    print(f"  After requiring survey design vars: {len(df):,}")
    
    # -------------------------------------------------------------------------
    # Create RIR-specific variables
    # -------------------------------------------------------------------------
    # LDL <70 mg/dL flag
    df['ldl_under_70'] = (df['ldl_calc'] < 70).astype(int)
    
    # LDL <55 mg/dL flag (sensitivity)
    df['ldl_under_55'] = (df['ldl_calc'] < 55).astype(int)
    
    # High hs-CRP ≥2 mg/L
    df['hscrp_high'] = (df['hscrp'] >= 2).astype(int)
    
    # Very high hs-CRP ≥3 mg/L (sensitivity)
    df['hscrp_very_high'] = (df['hscrp'] >= 3).astype(int)
    
    # Primary RIR definition: statin + LDL<70 + hs-CRP≥2
    df['rir'] = ((df['statin_user'] == 1) & 
                 (df['ldl_under_70'] == 1) & 
                 (df['hscrp_high'] == 1)).astype(int)
    
    # Alternative RIR with hs-CRP ≥3
    df['rir_strict'] = ((df['statin_user'] == 1) & 
                        (df['ldl_under_70'] == 1) & 
                        (df['hscrp_very_high'] == 1)).astype(int)
    
    # RIR with LDL <55 threshold
    df['rir_ldl55'] = ((df['statin_user'] == 1) & 
                       (df['ldl_under_55'] == 1) & 
                       (df['hscrp_high'] == 1)).astype(int)
    
    # -------------------------------------------------------------------------
    # PI-requested 6-group LDL × Statin stratification
    # -------------------------------------------------------------------------
    # LDL-C categories: ≤70, 70-130, >130
    def ldl_category(ldl):
        if pd.isna(ldl):
            return np.nan
        elif ldl <= 70:
            return 1  # ≤70
        elif ldl <= 130:
            return 2  # 70-130
        else:
            return 3  # >130
    
    df['ldl_cat'] = df['ldl_calc'].apply(ldl_category)
    
    # 6-group variable: LDL category × Statin status
    # 1 = LDL≤70, statin user
    # 2 = LDL≤70, no statin
    # 3 = LDL 70-130, statin user
    # 4 = LDL 70-130, no statin
    # 5 = LDL >130, statin user
    # 6 = LDL >130, no statin
    def ldl_statin_group(row):
        if pd.isna(row['ldl_cat']) or pd.isna(row['statin_user']):
            return np.nan
        ldl_c = int(row['ldl_cat'])
        statin = int(row['statin_user'])
        if ldl_c == 1:  # ≤70
            return 1 if statin == 1 else 2
        elif ldl_c == 2:  # 70-130
            return 3 if statin == 1 else 4
        else:  # >130
            return 5 if statin == 1 else 6
    
    df['ldl_statin_group'] = df.apply(ldl_statin_group, axis=1)
    
    # -------------------------------------------------------------------------
    # Summary statistics
    # -------------------------------------------------------------------------
    print(f"\n  Final cohort: {len(df):,}")
    print(f"  Statin users: {df['statin_user'].sum():,} ({100*df['statin_user'].mean():.1f}%)")
    print(f"  Statin users with LDL<70: {(df['statin_user'] & df['ldl_under_70']).sum():,}")
    print(f"  RIR cases (statin + LDL<70 + hs-CRP≥2): {df['rir'].sum():,}")
    
    return df


def main():
    """Main processing pipeline."""
    print("="*70)
    print("RIR STUDY - DATA PROCESSING")
    print("="*70)
    
    # Process all cycles
    all_data = []
    for cycle in CYCLES:
        cycle_data = process_cycle(cycle)
        if not cycle_data.empty:
            all_data.append(cycle_data)
    
    if not all_data:
        print("\nERROR: No data loaded. Run 01_download_data.py first.")
        return
    
    # Combine all cycles
    print("\n[Combining cycles]")
    df = pd.concat(all_data, ignore_index=True)
    print(f"  Combined N: {len(df):,}")
    
    # Harmonize variables
    df = harmonize_variables(df)
    
    # Create RIR cohort
    df = create_rir_cohort(df)
    
    # -------------------------------------------------------------------------
    # Select final analysis variables
    # -------------------------------------------------------------------------
    analysis_vars = [
        # Identifiers
        'SEQN', 'cycle',
        
        # Survey design
        'psu', 'strata', 'mec_weight', 'fasting_weight',
        
        # Demographics
        'age', 'sex', 'race_eth', 'education', 'pir',
        
        # Biomarkers
        'hscrp', 'ldl_calc', 'total_chol', 'hdl', 'trigly',
        'hba1c', 'fasting_glucose', 'sbp_mean', 'dbp_mean', 'bmi',
        
        # Clinical conditions
        'statin_user', 'diabetes', 'prediabetes', 'glycemic_status',
        'hypertension', 'current_smoker', 'smoking_status', 'obese',
        
        # RIR variables
        'ldl_under_70', 'ldl_under_55', 'hscrp_high', 'hscrp_very_high',
        'rir', 'rir_strict', 'rir_ldl55',
        
        # LDL × Statin stratification (PI request)
        'ldl_cat', 'ldl_statin_group',
    ]
    
    df_final = df[[c for c in analysis_vars if c in df.columns]].copy()
    
    # -------------------------------------------------------------------------
    # Save processed data
    # -------------------------------------------------------------------------
    DATA_PROCESSED.mkdir(parents=True, exist_ok=True)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    (OUTPUT_DIR / "tables").mkdir(exist_ok=True)
    (OUTPUT_DIR / "figures").mkdir(exist_ok=True)
    
    # Save as parquet for efficiency
    output_path = DATA_PROCESSED / "rir_analytic_cohort.parquet"
    df_final.to_parquet(output_path)
    print(f"\n[Saved] {output_path}")
    
    # Also save as CSV for R
    csv_path = DATA_PROCESSED / "rir_analytic_cohort.csv"
    df_final.to_csv(csv_path, index=False)
    print(f"[Saved] {csv_path}")
    
    # -------------------------------------------------------------------------
    # Summary report
    # -------------------------------------------------------------------------
    print("\n" + "="*70)
    print("PROCESSING COMPLETE - SUMMARY")
    print("="*70)
    print(f"\nTotal participants: {len(df_final):,}")
    print(f"\nBy cycle:")
    print(df_final.groupby('cycle').size().to_string())
    print(f"\nStatin users: {df_final['statin_user'].sum():,}")
    print(f"Statin users with LDL<70: {(df_final['statin_user'] & df_final['ldl_under_70']).sum():,}")
    print(f"\nPrimary RIR cases: {df_final['rir'].sum():,}")
    print(f"  (statin + LDL<70 + hs-CRP≥2)")
    
    # Create statin + LDL<70 subcohort for primary analysis
    statin_ldl70 = df_final[(df_final['statin_user'] == 1) & (df_final['ldl_under_70'] == 1)]
    if len(statin_ldl70) > 0:
        rir_prev = statin_ldl70['rir'].mean() * 100
        print(f"\nPrimary analysis cohort (statin + LDL<70):")
        print(f"  N = {len(statin_ldl70):,}")
        print(f"  RIR prevalence: {rir_prev:.1f}% (unweighted)")


if __name__ == "__main__":
    main()

