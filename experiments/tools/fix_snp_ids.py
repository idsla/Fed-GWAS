#!/usr/bin/env python3
"""
Fix SNP IDs in BIM files to rsXXX format

Since the VCF file doesn't contain rsIDs, we generate synthetic rsIDs
based on chromosome and position (sorted order).
"""

import pandas as pd
import sys
from pathlib import Path

def fix_snp_ids_to_rsformat(bim_file):
    """Fix SNP IDs to rsXXX format based on sorted chromosome and position"""
    path = Path(bim_file)
    if not path.exists():
        print(f'Skipping {bim_file} (not found)', file=sys.stderr)
        return False
    
    print(f'Fixing {bim_file}...')
    
    # Read BIM file
    bim = pd.read_csv(bim_file, sep='\t', header=None,
                     names=['CHR', 'SNP', 'CM', 'BP', 'A1', 'A2'])
    
    # Create a mapping from (CHR, BP) to rsID based on sorted order
    # Sort by CHR, then BP to ensure consistent ordering
    sorted_bim = bim.sort_values(['CHR', 'BP']).reset_index(drop=True)
    
    # Create rsID mapping: rs{1-based sequential number padded to 7 digits}
    rsid_map = {}
    for idx, row in sorted_bim.iterrows():
        key = (int(row['CHR']), int(row['BP']))
        rsid_map[key] = f'rs{(idx+1):07d}'
    
    # Apply mapping to original order
    bim['SNP'] = bim.apply(
        lambda row: rsid_map.get((int(row['CHR']), int(row['BP'])), f"{row['CHR']}:{row['BP']}"),
        axis=1
    )
    
    # Write back
    bim.to_csv(bim_file, sep='\t', header=False, index=False)
    print(f'  ✓ Fixed {len(bim)} SNPs')
    print(f'  Sample IDs: {bim["SNP"].head(3).tolist()}')
    return True


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python fix_snp_ids.py <bim_file1> [bim_file2] ...")
        sys.exit(1)
    
    success_count = 0
    for bim_file in sys.argv[1:]:
        if fix_snp_ids_to_rsformat(bim_file):
            success_count += 1
    
    print(f'\nFixed {success_count} BIM file(s)')
    sys.exit(0 if success_count > 0 else 1)
