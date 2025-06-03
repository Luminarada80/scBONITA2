import numpy as np
import random
import argparse
import logging
import sys
from tqdm import tqdm

# ────────────────────────────────────────────────────────────────────────────────
def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Evaluate scBONITA rules vs. simulated rules.")
    parser.add_argument(
        "--scbonita_rules",
        type=str,
        required=True,
        help="Path to scBONITA output rules file (one gene per line)."
    )
    parser.add_argument(
        "--simulated_rules",
        type=str,
        required=True,
        help="Path to “truth” rules file (one gene per line)."
    )
    parser.add_argument(
        "--data_file",
        type=str,
        required=True,
        help="Path to CSV of 0/1 expression data (header row with cells, first column gene names)."
    )
    parser.add_argument(
        "--output_file",
        type=str,
        default="output.txt",
        help="File to append summary results to."
    )
    parser.add_argument(
        "--n_reps",
        type=int,
        default=50,
        help="How many random replicates to run (each samples a random subset of cells)."
    )
    parser.add_argument(
        "--cell_frac",
        type=float,
        default=0.8,
        help="Fraction of cells to sample each replicate (0 < cell_frac ≤ 1)."
    )
    return parser.parse_args()


# ────────────────────────────────────────────────────────────────────────────────
def get_rules(path):
    """
    Parse a rules file where each line is:
      GeneX = GeneA AND GeneB OR NOT GeneC

    Returns:
      (ruleset, raw_ruletext)

      ruleset: list of list-of-tokens, e.g. [['0','and','1','or','not','2'], …]
      raw_ruletext: same length, each entry is original “GeneX = …” line stripped.
    """
    ruleset = []
    raw_ruletext = []
    with open(path, 'r') as rules_file:
        for lineno, line in enumerate(rules_file, start=1):
            line = line.strip()
            if not line or line.startswith('#'):
                continue  # skip blank lines or comments

            # Split on the first '='
            if '=' not in line:
                logging.warning(f"Line {lineno}: no '=' found → skipping: {line}")
                continue

            left, right = line.split('=', 1)
            left = left.strip()
            right = right.strip()

            # Warn if “left” is not “Gene…”
            if not left.lower().startswith('gene'):
                logging.warning(f"Line {lineno}: left side doesn’t start with 'Gene': {left}")

            raw_ruletext.append(line)

            # Tokenize the RHS by whitespace (strip parentheses)
            tokens = right.replace('(', ' ( ').replace(')', ' ) ').split()
            processed = []
            for tok in tokens:
                tok_lower = tok.lower()
                if tok_lower in {'and', 'or', 'not'}:
                    processed.append(tok_lower)
                elif tok_lower == '(' or tok_lower == ')':
                    continue
                else:
                    # assume “GeneN” → strip “gene” prefix
                    if tok_lower.startswith('gene'):
                        idx_str = tok_lower[4:]
                    else:
                        idx_str = tok_lower
                    if not idx_str.isdigit():
                        raise ValueError(f"Invalid gene token on line {lineno}: '{tok}'")
                    processed.append(idx_str)

            if len(processed) == 0:
                logging.warning(f"Line {lineno}: no valid tokens parsed, skipping.")
                continue

            ruleset.append(processed)

    if len(ruleset) != len(raw_ruletext):
        logging.debug("Warning: some lines were skipped; ruleset and raw_ruletext lengths may differ")

    return ruleset, raw_ruletext


# ────────────────────────────────────────────────────────────────────────────────
def get_matrix(path):
    """
    Read a CSV where first row is header (cell names),
    subsequent rows: geneName, 0/1, 0/1, …

    Returns:
      data_matrix: np.ndarray shape = (n_genes, n_cells), dtype=int
    """
    with open(path, 'r') as f:
        lines = [line.strip().split(',') for line in f]

    if len(lines) < 2:
        raise ValueError(f"Data file {path} has no data rows")

    header = lines[0][1:]  # skip the first column (gene name); these are cell names
    data_rows = lines[1:]
    n_cells = len(header)

    matrix = []
    for row in data_rows:
        if len(row) < n_cells + 1:
            raise ValueError(f"Row too short: {row}")
        try:
            nums = [int(x) for x in row[1 : 1 + n_cells]]
        except ValueError:
            raise ValueError(f"Non-integer value in row: {row}")
        matrix.append(nums)

    mat = np.array(matrix, dtype=int)
    if mat.ndim != 2:
        raise ValueError("Converted data is not 2D")
    return mat


# ────────────────────────────────────────────────────────────────────────────────
def check_logic(rule, state):
    """
    Evaluate a single rule (list of tokens) against a 1D state array (0/1 for each gene).
    No parentheses or real precedence: tokens are evaluated left→right, with 'not' applying only
    to the very next token and 'and'/'or' combining the accumulated result.

    Returns:
      1 or 0
    """
    result = None
    operator = None

    for token in rule:
        if token in {'and', 'or', 'not'}:
            operator = token
            continue

        idx = int(token)
        if idx < 0 or idx >= len(state):
            raise IndexError(f"Rule references gene index {idx}, but state length = {len(state)}")

        value = bool(state[idx])
        if operator == 'not':
            value = not value
            operator = None

        if result is None:
            result = value
        else:
            if operator == 'and':
                result = result and value
            elif operator == 'or':
                result = result or value
            else:
                # no operator means just override
                result = value

        operator = None

    if result is None:
        # no tokens parsed
        return 0
    return 1 if result else 0


# ────────────────────────────────────────────────────────────────────────────────
def check_rules(dataset, rules, raw_rule_text, cell_indices):
    """
    For the given subset of cells (columns = cell_indices) in dataset,
    evaluate each rule on each chosen cell, count mismatches.

    Args:
      dataset      : np.ndarray of shape (n_genes, n_cells_full)
      rules        : list of rule‐token‐lists (length = n_genes)
      raw_rule_text: list of the original rule strings
      cell_indices : list/array of column indices to use in this replicate

    Returns:
      mismatch_count, total_checks, error_pct, num_rules
    """
    n_genes = dataset.shape[0]
    mismatch_count = 0
    total_checks = 0
    rule_errors = {}

    for cell_idx in cell_indices:
        state = dataset[:, cell_idx]
        for gene_idx, rule in enumerate(rules):
            key = f"Gene{gene_idx}: {raw_rule_text[gene_idx]}"
            rule_errors.setdefault(key, 0)

            predicted = check_logic(rule, state)
            actual = int(state[gene_idx])
            if predicted != actual:
                mismatch_count += 1
                rule_errors[key] += 1
            total_checks += 1

    error_pct = 100.0 * mismatch_count / total_checks if total_checks > 0 else 0.0
    return mismatch_count, total_checks, error_pct, len(rules)


# ────────────────────────────────────────────────────────────────────────────────
def evaluate_scbonita_output():
    args = parse_args()

    # 1) Load both rule files (static; do not change across replicates)
    scbonita_ruleset, scbonita_raw = get_rules(args.scbonita_rules)
    simulated_ruleset, simulated_raw = get_rules(args.simulated_rules)

    # 2) Load the full expression matrix once
    matrix = get_matrix(args.data_file)
    n_genes, n_cells_full = matrix.shape
    print(f"Loaded data matrix: {n_genes} genes × {n_cells_full} cells")

    # 3) Decide how many replicates and what fraction of cells per replicate
    n_reps = args.n_reps
    frac = args.cell_frac
    if not (0 < frac <= 1):
        raise ValueError("cell_frac must be between 0 and 1")

    # 4) Prepare output file header
    with open(args.output_file, 'a') as f:
        f.write("replicate,rule_set,num_rules,mismatches,total_checks,error_pct\n")

    # 5) Run the replicates
    for i in tqdm(range(n_reps), desc="Replicates"):
        # Randomly sample (frac * n_cells_full) distinct cell indices
        k = max(1, int(n_cells_full * frac))
        cell_indices = np.random.choice(n_cells_full, size=k, replace=False)

        # Evaluate scBONITA rules on this subset
        sc_m, sc_t, sc_e, sc_num = check_rules(matrix, scbonita_ruleset, scbonita_raw, cell_indices)
        with open(args.output_file, 'a') as f:
            f.write(f"{i},scBONITA,{sc_num},{sc_m},{sc_t},{sc_e:.2f}\n")

        # Evaluate simulated “ground truth” rules on the same subset
        sim_m, sim_t, sim_e, sim_num = check_rules(matrix, simulated_ruleset, simulated_raw, cell_indices)
        with open(args.output_file, 'a') as f:
            f.write(f"{i},simulated,{sim_num},{sim_m},{sim_t},{sim_e:.2f}\n")

    print(f"\nFinished {n_reps} replicates. Results appended to {args.output_file}")


# ────────────────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    evaluate_scbonita_output()
