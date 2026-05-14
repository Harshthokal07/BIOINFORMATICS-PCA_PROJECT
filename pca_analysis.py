# PCA Analysis - Breast Cancer Gene Expression (GSE5325 Dataset)
# Reproducing Figure 1 from Nature Biotechnology PCA Primer
# Using XBP1 and GATA3 gene expression across 105 patients

import matplotlib
matplotlib.use('Agg')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# --- Loading the dataset ---

# Read patient class labels (ER+ = 1, ER- = 0)
patient_labels = pd.read_csv('data/class.tsv', header=None, names=['class'])
print(f"Number of patients: {len(patient_labels)}")
print(f"ER+ patients: {(patient_labels['class'] == 1).sum()}")
print(f"ER- patients: {(patient_labels['class'] == 0).sum()}")

# Read gene expression data (patients as rows, genes as columns)
expr_data = pd.read_csv('data/filtered.tsv.gz', sep='\t', compression='gzip')

# Strip any whitespace from column headers
expr_data.columns = expr_data.columns.str.strip()
print(f"\nExpression matrix shape: {expr_data.shape}")
print(f"(rows=patients, columns=genes)")


# --- Extracting the two genes of interest ---
# Gene mapping (from columns.tsv.gz):
#   ID 4404 -> XBP1 (X-box binding protein 1)
#   ID 4359 -> GATA3 (GATA binding protein 3)

xbp1_expr = expr_data['4404'].values
gata3_expr = expr_data['4359'].values
sample_class = patient_labels['class'].values

print(f"\nXBP1 expression range: [{xbp1_expr.min():.4f}, {xbp1_expr.max():.4f}]")
print(f"GATA3 expression range: [{gata3_expr.min():.4f}, {gata3_expr.max():.4f}]")


# --- Identifying patient groups ---
mask_positive = sample_class == 1  # ER+ patients
mask_negative = sample_class == 0  # ER- patients


# --- Figure 1a: Scatter plot (GATA3 vs XBP1) ---

f1, a1 = plt.subplots(figsize=(6, 5))

# ER- shown as black squares, ER+ as red squares (matching the reference paper)
a1.scatter(gata3_expr[mask_negative], xbp1_expr[mask_negative],
           c='black', s=30, marker='s', label='ER- (class 0)', zorder=2)
a1.scatter(gata3_expr[mask_positive], xbp1_expr[mask_positive],
           c='red', s=30, marker='s', label='ER+ (class 1)', zorder=3)

a1.set_xlabel('GATA3', fontsize=13, fontstyle='italic')
a1.set_ylabel('XBP1', fontsize=13, fontstyle='italic')
a1.set_title('Figure 1a: XBP1 vs GATA3 Expression', fontsize=14, fontweight='bold')
a1.legend(fontsize=10)
a1.grid(False)
a1.spines['top'].set_visible(False)
a1.spines['right'].set_visible(False)

f1.tight_layout()
f1.savefig('figure_1a.png', dpi=150, bbox_inches='tight')
print("\nFigure 1a saved as 'figure_1a.png'")


# --- Performing PCA manually on the 2-gene data ---

# Stack GATA3 (x) and XBP1 (y) into a 105x2 matrix
combined = np.column_stack([gata3_expr, xbp1_expr])

# Step 1: Center the data by removing the mean
avg = combined.mean(axis=0)
centered = combined - avg

# Step 2: Covariance matrix of centered data
cov_mat = np.cov(centered.T)
print(f"\nCovariance matrix:\n{cov_mat}")

# Step 3: Eigendecomposition
eig_vals, eig_vecs = np.linalg.eigh(cov_mat)

# Step 4: Sort eigenpairs (largest eigenvalue first)
sort_order = np.argsort(eig_vals)[::-1]
eig_vals = eig_vals[sort_order]
eig_vecs = eig_vecs[:, sort_order]

# PC1 and PC2 directions
pc1_vec = eig_vecs[:, 0]
pc2_vec = eig_vecs[:, 1]

# Flip PC1 direction so ER+ projects to positive side (eigenvectors have arbitrary sign)
pc1_vec = -pc1_vec

print(f"\nPC1 direction: {pc1_vec}")
print(f"PC2 direction: {pc2_vec}")
print(f"Eigenvalues: {eig_vals}")
print(f"Variance explained by PC1: {eig_vals[0]/eig_vals.sum()*100:.1f}%")


# --- Figure 1b: PCA direction arrows overlaid on scatter ---

f1b, a1b = plt.subplots(figsize=(6, 5))

# Replot the scatter points
a1b.scatter(gata3_expr[mask_negative], xbp1_expr[mask_negative],
            c='black', s=30, marker='s', zorder=2)
a1b.scatter(gata3_expr[mask_positive], xbp1_expr[mask_positive],
            c='red', s=30, marker='s', zorder=3)

# Arrow scale factor for visualization
arrow_len = 3.5

# Draw PC1 arrow through the data mean
a1b.annotate('', xy=(avg[0] + arrow_len * pc1_vec[0],
                      avg[1] + arrow_len * pc1_vec[1]),
             xytext=(avg[0] - arrow_len * pc1_vec[0],
                     avg[1] - arrow_len * pc1_vec[1]),
             arrowprops=dict(arrowstyle='->', lw=1.5, color='black'))
a1b.annotate('PC1', xy=(avg[0] + arrow_len * pc1_vec[0] + 0.15,
                         avg[1] + arrow_len * pc1_vec[1] + 0.15),
             fontsize=12, fontweight='bold')

# Draw PC2 arrow through the data mean
a1b.annotate('', xy=(avg[0] + arrow_len * pc2_vec[0],
                      avg[1] + arrow_len * pc2_vec[1]),
             xytext=(avg[0] - arrow_len * pc2_vec[0],
                     avg[1] - arrow_len * pc2_vec[1]),
             arrowprops=dict(arrowstyle='->', lw=1.5, color='black'))
a1b.annotate('PC2', xy=(avg[0] + arrow_len * pc2_vec[0] + 0.15,
                         avg[1] + arrow_len * pc2_vec[1] + 0.15),
             fontsize=12, fontweight='bold')

a1b.set_xlabel('GATA3', fontsize=13, fontstyle='italic')
a1b.set_ylabel('XBP1', fontsize=13, fontstyle='italic')
a1b.set_title('Figure 1b: PCA Directions on Scatter Plot', fontsize=14, fontweight='bold')
a1b.grid(False)
a1b.spines['top'].set_visible(False)
a1b.spines['right'].set_visible(False)

f1b.tight_layout()
f1b.savefig('figure_1b.png', dpi=150, bbox_inches='tight')
print("\nFigure 1b saved as 'figure_1b.png'")


# --- Projecting data onto PC1 ---

pc1_scores = centered @ pc1_vec  # dot product of each centered point with PC1
print(f"\nPC1 projection range: [{pc1_scores.min():.4f}, {pc1_scores.max():.4f}]")


# --- Figure 1c: 1D strip plot along PC1 ---

f2, ax_arr = plt.subplots(3, 1, figsize=(8, 4), sharex=True,
                           gridspec_kw={'hspace': 0.3})

# Helper to clean up each strip subplot
def format_strip(ax):
    ax.set_yticks([])
    ax.axhline(y=0, color='gray', linewidth=0.5, linestyle='-')
    for spine in ['top', 'right', 'left']:
        ax.spines[spine].set_visible(False)

# Top row: all patients together
ax_arr[0].scatter(pc1_scores[mask_negative], [0]*mask_negative.sum(),
                  c='black', s=30, marker='s', zorder=2)
ax_arr[0].scatter(pc1_scores[mask_positive], [0]*mask_positive.sum(),
                  c='red', s=30, marker='s', zorder=3)
ax_arr[0].set_ylabel('All', fontsize=12, fontweight='bold')
format_strip(ax_arr[0])

# Middle row: ER- only
ax_arr[1].scatter(pc1_scores[mask_negative], [0]*mask_negative.sum(),
                  c='black', s=30, marker='s', zorder=2)
ax_arr[1].set_ylabel('ER$^-$', fontsize=12, fontweight='bold')
format_strip(ax_arr[1])

# Bottom row: ER+ only
ax_arr[2].scatter(pc1_scores[mask_positive], [0]*mask_positive.sum(),
                  c='red', s=30, marker='s', zorder=3)
ax_arr[2].set_ylabel('ER$^+$', fontsize=12, fontweight='bold')
format_strip(ax_arr[2])
ax_arr[2].set_xlabel('Projection onto PC1', fontsize=13)

f2.suptitle('Figure 1c: Projection onto PC1', fontsize=14, fontweight='bold', y=1.02)
f2.tight_layout()
f2.savefig('figure_1c.png', dpi=150, bbox_inches='tight')
print("Figure 1c saved as 'figure_1c.png'")

print("\n=== Analysis Complete! ===")
