"""
Author: Jyotirmoy Das
Maintainer: Bioinformatics Core Facility and Clinical Genomics Linkoping
Version: 0.2.3
Purpose: methylation dashboard data injection from gms-nallo result.
"""

import os
import json
import csv
import re
import argparse
import sys
from datetime import datetime
from pathlib import Path

def print_banner():
    current_time = datetime.now().strftime("%Y-%m-%d %H:%M")
    banner = f"""
    =======================================================
     Methylation Dashboard Manager for gms-nallo 
    =======================================================
    Author: Jyotirmoy Das
    Run Date: {current_time}
    =======================================================
    """
    print(banner)

def parse_methbat_tsv(tsv_path):
    """Parses the methbat profile TSV/CSV into a list of dicts for AG-Grid."""
    rows = []
    try:
        # Detect delimiter based on extension
        delimiter = ',' if tsv_path.suffix == '.csv' else '\t'
        with open(tsv_path, 'r', encoding='utf-8-sig') as f:
            reader = csv.reader(f, delimiter=delimiter)
            header = None
            for row in reader:
                if not row or row[0].startswith('#'):
                    continue
                header = [h.strip().lower() for h in row]
                break
            
            if not header:
                return []

            # Map indices
            def get_idx(name):
                for i, h in enumerate(header):
                    if name.lower() in h:
                        return i
                return -1

            idx_chrom = get_idx("chrom")
            idx_start = get_idx("start")
            idx_end = get_idx("end")
            idx_cpg = get_idx("cpg_label")
            idx_exp_cat = get_idx("expected_category")
            idx_sum_lab = get_idx("summary_label")
            idx_rep_sum = get_idx("report_summary")
            idx_qc = get_idx("qc_warnings")
            idx_combined = get_idx("mean_combined_methyl")
            idx_delta = get_idx("mean_meth_delta")
            idx_h1 = get_idx("mean_hap1_methyl")
            idx_h2 = get_idx("mean_hap2_methyl")

            for row in reader:
                if not row or row[0].startswith('#'):
                    continue
                
                try:
                    start_val = row[idx_start] if idx_start != -1 else "0"
                    end_val = row[idx_end] if idx_end != -1 else "0"
                    start = int(start_val)
                    end = int(end_val)
                    interval = end - start
                except (ValueError, TypeError, IndexError):
                    interval = 0

                def safe_get(idx, default=""):
                    try:
                        return row[idx].strip() if idx != -1 else default
                    except IndexError:
                        return default

                # Combine report_summary and qc_warnings
                summary = safe_get(idx_rep_sum) or "PASS"
                warnings = safe_get(idx_qc) or "PASS"
                
                # 0: PASS, 1: One fail, 2: Both fail
                qc_status = 0
                if summary == "PASS" and warnings == "PASS":
                    qc = "PASS"
                    qc_status = 0
                elif summary != "PASS" and warnings == "PASS":
                    qc = summary
                    qc_status = 1
                elif summary == "PASS" and warnings != "PASS":
                    qc = warnings
                    qc_status = 1
                else:
                    qc = f"{summary}; {warnings}"
                    qc_status = 2

                # Status logic: use summary_label as requested
                status = safe_get(idx_sum_lab) or safe_get(idx_exp_cat) or "Uncategorized"
                
                combined_raw = safe_get(idx_combined, "0")
                combined_val = float(combined_raw) if combined_raw else 0.0

                cleaned_row = {
                    "chrom": safe_get(idx_chrom),
                    "start": safe_get(idx_start, "0"),
                    "end": safe_get(idx_end, "0"),
                    "interval": str(interval),
                    "cpg_label": safe_get(idx_cpg, "Unknown"),
                    "expected_category": safe_get(idx_exp_cat, "Uncategorized"),
                    "summary_label": safe_get(idx_sum_lab, "Uncategorized"),
                    "status_label": status,
                    "qc_warnings": qc,
                    "qc_status": qc_status,
                    "mean_combined_methyl": str(round(combined_val, 3)),
                    "meth_score": str(round(combined_val, 3)),
                    "unmeth_score": str(round(1.0 - combined_val, 3)),
                    "mean_meth_delta": str(round(float(safe_get(idx_delta, "0") or 0), 3)),
                    "mean_hap1_methyl": str(round(float(safe_get(idx_h1, "0") or 0), 3)),
                    "mean_hap2_methyl": str(round(float(safe_get(idx_h2, "0") or 0), 3))
                }
                rows.append(cleaned_row)
    except Exception as e:
        print(f" [!] Error parsing {tsv_path.name}: {e}")
    return rows

def calculate_stats(rows, sid):
    """Calculates summary statistics for the dashboard cards."""
    if not rows:
        return {"cov": "N/A", "h1Avg": "0%", "h2Avg": "0%", "passCount": 0}
    
    h1_vals = [float(r['mean_hap1_methyl']) for r in rows if r['mean_hap1_methyl']]
    h2_vals = [float(r['mean_hap2_methyl']) for r in rows if r['mean_hap2_methyl']]
    
    h1_avg = (sum(h1_vals) / len(h1_vals) * 100) if h1_vals else 0
    h2_avg = (sum(h2_vals) / len(h2_vals) * 100) if h2_vals else 0
    
    pass_count = sum(1 for r in rows if r.get('qc_status') == 0)
    
    return {
        "cov": "Nallo-Phased",
        "h1Avg": f"{h1_avg:.1f}%",
        "h2Avg": f"{h2_avg:.1f}%",
        "passCount": pass_count
    }

def main():
    print_banner()
    
    parser = argparse.ArgumentParser(
        description="gms-nallo Methylation Dashboard.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example Usage:
  python3 nallo_methylation_manager.py --results ./results_260409 --template dashboard-template.html --output My_Audit_Report.html
        """
    )
    
    # Auto-generate date for default output
    default_date = datetime.now().strftime("%y%m%d")
    default_name = f"GMS_Nallo_Methylation_Report_{default_date}.html"

    parser.add_argument("--results", help="Directory containing nallo pipeline results (e.g. results)")
    parser.add_argument("--template", help="Path to the HTML template file (dashboard-template.html)")
    parser.add_argument("--output", default=default_name, help=f"Path for final report (default: {default_name})")
    
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
        
    args = parser.parse_args()
    
    results_dir = Path(args.results)
    template_path = Path(args.template)
    
    # Validation
    if not results_dir.exists():
        print(f"\n [!] ERROR: Results directory does not exist: {args.results}")
        sys.exit(1)
    if not template_path.exists():
        print(f"\n [!] ERROR: HTML template not found: {args.template}")
        sys.exit(1)

    sample_data = {}
    sample_stats = {}
    
    print(f" [+] Scanning results: {results_dir.name}")
    
    # Search in multiple locations
    search_dirs = [
        results_dir / "methylation" / "report",
        results_dir / "tables",
        results_dir
    ]
    
    for base_dir in search_dirs:
        if not base_dir.exists():
            continue
            
        # Check subdirectories (original logic)
        for entry in base_dir.iterdir():
            if entry.is_dir():
                sid = entry.name
                if sid in sample_data or sid in ["tables", "figures", "reports", "methylation"]: 
                    continue
                
                target_file = None
                candidates = [
                    entry / f"{sid}_methbat_profile.tsv",
                    entry / f"{sid}.methbat_mod1.csv",
                    entry / f"{sid}_methbat_profile.csv"
                ]
                for c in candidates:
                    if c.exists():
                        target_file = c
                        break
                if not target_file:
                    for f in entry.glob("*.methbat*"):
                        if f.suffix in ['.csv', '.tsv']:
                            target_file = f
                            break
                if target_file:
                    print(f" [+] Found Sample Profile: {target_file.name} in {entry.name}")
                    rows = parse_methbat_tsv(target_file)
                    if rows:
                        sample_data[sid] = rows
                        sample_stats[sid] = calculate_stats(rows, sid)
            
            # Also check files directly in base_dir if they match the pattern
            elif entry.is_file() and (".methbat" in entry.name):
                # Try to extract SID from filename
                sid = entry.name.split('.')[0].split('_')[0]
                if sid not in sample_data:
                    print(f" [+] Found Sample Profile File: {entry.name}")
                    rows = parse_methbat_tsv(entry)
                    if rows:
                        sample_data[sid] = rows
                        sample_stats[sid] = calculate_stats(rows, sid)

    if not sample_data:
        print("\n [!] ERROR: No valid sample profiles were identified.")
        sys.exit(1)

    # Read the template
    try:
        with open(template_path, 'r', encoding='utf-8') as f:
            template_html = f.read()

        # Output directory: if --output is a file, use its directory. 
        # If it's a directory, use it directly.
        output_path = Path(args.output)
        if output_path.suffix:
            out_dir = output_path.parent
        else:
            out_dir = output_path
            out_dir.mkdir(parents=True, exist_ok=True)

        print(f"\n [+] Generating separate reports in: {out_dir}")

        for sid, rows in sample_data.items():
            # Injection logic for single sample
            sample_report_html = template_html.replace('{{SAMPLE_DATA_JSON}}', json.dumps({sid: rows}))
            sample_report_html = sample_report_html.replace('{{SAMPLE_STATS_JSON}}', json.dumps({sid: sample_stats[sid]}, indent=8))

            sample_output = out_dir / f"{sid}_methylation_report.html"
            with open(sample_output, 'w', encoding='utf-8') as f:
                f.write(sample_report_html)
            
            print(f" [SUCCESS] Report for {sid} saved to: {sample_output}")

        print(f"\n [+] Total Samples Processed: {len(sample_data)}")
        print(f" [+] Total Loci Processed: {sum(len(v) for v in sample_data.values())}")
        print("-" * 55)

    except Exception as e:
        print(f"\n [!] Critical Error during report generation: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
