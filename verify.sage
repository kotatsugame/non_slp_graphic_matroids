import csv
import os
import sys

# Increase recursion depth for complex polynomial computations
sys.setrecursionlimit(20000)

# Directory configurations
DIR_GRAPHS  = 'graphs'
DIR_POLYS   = 'polynomial'
DIR_HESSIAN = 'H'
DIR_KERNEL  = 'F'
DIR_BASIS   = 'basis'

def load_poly_from_txt(path, R):
    """Load a polynomial from a text file and parse it in the given ring R."""
    with open(path) as f:
        return R(f.read())

def verify_graph_consistency(index, g6_str, binary_str):
    """
    Verify that the graph data from three different sources are identical:
    1. graph6 string in table.csv
    2. Edge list in the graphs/ directory
    3. Binary edge set representation
    """
    g_ref = Graph(g6_str)
    
    # 1. Reconstruct from edge list file
    g_file = Graph(8)
    with open(os.path.join(DIR_GRAPHS, f'graphs_{index}.txt')) as f:
        num_edges = int(f.readline())
        for _ in range(num_edges):
            u, v = f.readline().split()
            g_file.add_edge(int(u), int(v))
            
    # 2. Reconstruct from binary string (subset of all possible edges)
    g_bin = Graph(8)
    edges_all = Combinations(range(8), 2).list()
    for i, is_present in enumerate(binary_str):
        if is_present == '1':
            g_bin.add_edge(edges_all[i])
            
    assert g_ref == g_file == g_bin
    return g_ref

def verify_basis_and_hessian(f_G, basis_indices, H_actual, R):
    """Verify that the stored Hessian matrix matches the derivatives of f_G by the basis."""
    m = len(basis_indices)
    for i in range(m):
        # First derivative by basis element B[i]
        df_i = f_G
        for var_idx in basis_indices[i]:
            df_i = df_i.derivative(R.gen(var_idx))
            
        for j in range(m):
            # Second derivative by basis element B[j]
            entry = df_i
            for var_idx in basis_indices[j]:
                entry = entry.derivative(R.gen(var_idx))
            # Compare with the stored Hessian entry
            assert entry == H_actual[i, j]

def run_verification(entry):
    """Run the full suite of algebraic and combinatorial verifications for a single graph."""
    index = int(entry['index'])
    n_edges = int(entry['num_edges'])
    m_dim = int(entry['dim_a3'])
    R = PolynomialRing(QQ, 'x', n_edges)
    
    print(f"Checking Index {index:03d}: (num_edges={n_edges}, dim_A3={m_dim})")

    # 1. Graph and Matroid consistency check
    G = verify_graph_consistency(index, entry['graph6'], entry['edge_set_binary'])
    M = Matroid(G)
    
    # 2. Basis Generating Function (f_G) verification
    f_G_stored = load_poly_from_txt(os.path.join(DIR_POLYS, f'polynomial_{index}.txt'), R)
    e_to_x = dict(zip(sorted(M.groundset()), R.gens()))
    f_G_computed = sum(prod(e_to_x[e] for e in b) for b in M.bases())
    assert f_G_stored == f_G_computed

    # 3. Hessian Matrix (H) loading
    H = matrix(R, m_dim)
    with open(os.path.join(DIR_HESSIAN, f'H_{index}.csv')) as f:
        for i, row in enumerate(csv.reader(f)):
            for j, cell in enumerate(row):
                H[i, j] = R(cell)

    # 4. Kernel Vector (F) verification
    # Note: Skip verification for cases where degree is 10 (estimated cases)
    kernel_files = []
    if entry['deg_f_1'] != '10':
        kernel_files.append(os.path.join(DIR_KERNEL, f'F_{index}.txt'))
    if entry['deg_f_2']:
        kernel_files.append(os.path.join(DIR_KERNEL, f'F_{index}_2.txt'))

    for fpath in kernel_files:
        with open(fpath) as f:
            # Reconstruct kernel vector F and check H * F == 0
            F = vector(R, [R(line) for line in f])
            assert not F.is_zero(), f"Kernel vector F in {fpath} is zero"
            assert (H * F).is_zero(), f"H * F != 0 for index {index}"

    # 5. Basis (B) consistency with Hessian matrix
    basis_indices = []
    with open(os.path.join(DIR_BASIS, f'basis_{index}.csv')) as f:
        for row in csv.reader(f):
            basis_indices.append([int(x) for x in row])
    
    verify_basis_and_hessian(f_G_stored, basis_indices, H, R)

	# 6. Basis Selection Validity (Checking linear independence via pivot rows)
    comb_3 = Combinations(range(n_edges), 3).list()
    derivs = []
    all_monomials = set()

    for vars_idx in comb_3:
        d = f_G_stored
        for v in vars_idx:
            d = d.derivative(R.gen(v))
        derivs.append(d)
        # Collect all monomials present in these derivatives
        all_monomials.update(d.monomials())

    # Sort monomials to ensure a consistent matrix representation
    sorted_monomials = sorted(list(all_monomials))
    mono_to_idx = {m: i for i, m in enumerate(sorted_monomials)}

    # Construct the coefficient matrix
    matrix_rows = []
    for d in derivs:
        row = [0] * len(sorted_monomials)
        for coeff, mono in d:
            row[mono_to_idx[mono]] = coeff
        matrix_rows.append(row)

    mat_derivs = matrix(QQ, matrix_rows)
    pivots = mat_derivs.pivot_rows()
    assert basis_indices == [comb_3[p] for p in pivots], f"Basis mismatch for index {index}"

# --- Main Execution ---
if __name__ == "__main__":
    try:
        with open('table.csv') as f:
            for entry in csv.DictReader(f):
                run_verification(entry)
        print("\n" + "="*40)
        print("SUCCESS: All verifications passed.")
        print("="*40)
    except Exception as e:
        print(f"\nVerification failed: {e}")
        sys.exit(1)
