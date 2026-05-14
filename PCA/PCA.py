import gzip
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

def main():
    # 1. Load the classes
    with open('data/class.tsv', 'r') as f:
        classes = [int(line.strip()) for line in f.readlines() if line.strip() != '']

    # 2. Extract XBP1 (4404) and GATA3 (4359) expression levels
    xbp1_idx = -1
    gata3_idx = -1

    xbp1_expr = []
    gata3_expr = []

    with gzip.open('data/filtered.tsv.gz', 'rt') as f:
        # Read header to get column indices
        header_line = f.readline()
        header = [x.strip() for x in header_line.split('\t')]
        
        try:
            xbp1_idx = header.index('4404')
            gata3_idx = header.index('4359')
        except ValueError as e:
            print(f"Error finding gene in header: {e}")
            return
            
        # Read the rest of the samples (105 lines expected)
        for line in f:
            parts = line.split('\t')
            # Extract expression levels and convert to float
            xbp1_expr.append(float(parts[xbp1_idx]))
            gata3_expr.append(float(parts[gata3_idx]))

    xbp1 = np.array(xbp1_expr)
    gata3 = np.array(gata3_expr)
    classes = np.array(classes)
    
    print(f"Loaded {len(xbp1)} samples.")

    # 3. Create Scatter plot (Figure 1a)
    plt.figure(figsize=(5, 5))

    er_pos = classes == 1
    er_neg = classes == 0

    # Nature PCA primer uses red for ER+ and black for ER-
    plt.scatter(gata3[er_pos], xbp1[er_pos], color='red', label='ER+', alpha=0.8, s=15)
    plt.scatter(gata3[er_neg], xbp1[er_neg], color='black', label='ER-', alpha=0.8, s=15)

    plt.xlabel("GATA3 Expression")
    plt.ylabel("XBP1 Expression")
    plt.title("Figure 1a: GATA3 vs XBP1 Expression")
    plt.legend()
    # Remove top/right spines for cleaner look
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.savefig('fig1a.png', dpi=300, bbox_inches='tight')
    plt.close()

    # 4. Run PCA on the 2D matrix
    # Construct a matrix of shape (105, 2)
    X = np.column_stack((gata3, xbp1))

    pca = PCA(n_components=2)
    X_pca = pca.fit_transform(X)
    
    pc1_var = pca.explained_variance_ratio_[0] * 100
    pc2_var = pca.explained_variance_ratio_[1] * 100
    
    print(f"PC1 Variance: {pc1_var:.2f}%")
    print(f"PC2 Variance: {pc2_var:.2f}%")

    # Figure 1b: Scatter plot with Principal Component axes
    plt.figure(figsize=(5, 5))
    
    plt.scatter(gata3[er_pos], xbp1[er_pos], color='red', label='ER+', alpha=0.8, s=15)
    plt.scatter(gata3[er_neg], xbp1[er_neg], color='black', label='ER-', alpha=0.8, s=15)
    
    mean_gata3 = np.mean(gata3)
    mean_xbp1 = np.mean(xbp1)
    
    std_pc1 = np.sqrt(pca.explained_variance_[0])
    std_pc2 = np.sqrt(pca.explained_variance_[1])
    
    pc1_vector = pca.components_[0] * std_pc1 * 2
    pc2_vector = pca.components_[1] * std_pc2 * 2

    # Draw PC1 axis
    plt.annotate('', xy=(mean_gata3 + pc1_vector[0], mean_xbp1 + pc1_vector[1]), 
                 xytext=(mean_gata3 - pc1_vector[0], mean_xbp1 - pc1_vector[1]),
                 arrowprops=dict(arrowstyle="-", color='gray', linewidth=2))
                 
    # Draw PC2 axis
    plt.annotate('', xy=(mean_gata3 + pc2_vector[0], mean_xbp1 + pc2_vector[1]), 
                 xytext=(mean_gata3 - pc2_vector[0], mean_xbp1 - pc2_vector[1]),
                 arrowprops=dict(arrowstyle="-", color='gray', linewidth=2))

    plt.xlabel("GATA3 Expression")
    plt.ylabel("XBP1 Expression")
    plt.title("Figure 1b: GATA3 vs XBP1 with PC Axes")
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.savefig('fig1b.png', dpi=300, bbox_inches='tight')
    plt.close()

    # 5. Create 1D projection plot (Figure 1c)
    plt.figure(figsize=(4, 4))
    
    # We want 3 rows: y=0 for ER+, y=1 for ER-, y=2 for All
    
    # Draw horizontal lines for each row
    plt.axhline(0, color='gray', linewidth=1, zorder=0)
    plt.axhline(1, color='gray', linewidth=1, zorder=0)
    plt.axhline(2, color='gray', linewidth=1, zorder=0)

    # Plot ER+ at y=0
    plt.scatter(X_pca[er_pos, 0], np.zeros(sum(er_pos)), color='red', s=15, zorder=2)
    # Plot ER- at y=1
    plt.scatter(X_pca[er_neg, 0], np.ones(sum(er_neg)), color='black', s=15, zorder=2)
    
    # Plot All at y=2
    plt.scatter(X_pca[er_pos, 0], np.full(sum(er_pos), 2), color='red', s=15, zorder=2)
    plt.scatter(X_pca[er_neg, 0], np.full(sum(er_neg), 2), color='black', s=15, zorder=2)

    plt.xlabel("Projection onto PC1")
    plt.yticks([0, 1, 2], ['ER$^+$', 'ER$^-$', 'All'])
    
    # Remove top and right spines
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    # Also remove left spine to match the look more closely, if desired
    # plt.gca().spines['left'].set_visible(False) 
    
    # Don't need title/legend for this specific replicate if we just want it exactly like the image
    # plt.title("Figure 1c: 1D Projection")
    
    plt.tight_layout()
    plt.savefig('fig1c.png', dpi=300, bbox_inches='tight')
    plt.close()

    print("Plots generated successfully: fig1a.png, fig1b.png, fig1c.png")

if __name__ == "__main__":
    main()
