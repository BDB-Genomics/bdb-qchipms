import matplotlib.pyplot as plt
import numpy as np

fig, ax = plt.subplots(figsize=(8, 5), dpi=300)
fig.patch.set_facecolor('#ffffff')
ax.set_facecolor('#f8f9fa')

np.random.seed(42)
igg = np.random.normal(19.0, 0.4, 500)
trf2 = np.random.normal(28.0, 0.4, 500)

ax.hist(igg, bins=30, alpha=0.6, color='#3182bd', label='IgG Control (Low/Noise)')
ax.hist(trf2, bins=30, alpha=0.6, color='#de2d26', label='TERF2 Target (High/Signal)')

ax.axvline(np.mean(igg), color='#08519c', linestyle='--', linewidth=2, label=f'Mean IgG ({np.mean(igg):.1f})')
ax.axvline(np.mean(trf2), color='#a50f15', linestyle='--', linewidth=2, label=f'Mean TERF2 ({np.mean(trf2):.1f})')

ax.annotate('', xy=(np.mean(trf2), 35), xytext=(np.mean(igg), 35),
            arrowprops=dict(arrowstyle='<->', color='#252525', lw=2))
ax.text((np.mean(igg) + np.mean(trf2))/2, 37, 'SIGNAL (Log2FC Difference = +9.0)',
        ha='center', va='bottom', fontsize=10, fontweight='bold', color='#252525')

ax.set_xlabel('Log2 LFQ Intensity', fontsize=11, fontweight='bold')
ax.set_ylabel('Frequency (Protein Replicates)', fontsize=11, fontweight='bold')
ax.set_title("Welch's t-Test Concept: High Signal vs Low Noise (p < 0.0001)", fontsize=13, fontweight='bold', pad=15)
ax.legend(loc='upper right', frameon=True, facecolor='#ffffff', edgecolor='#cccccc')

plt.tight_layout()
plt.savefig('/home/himanshu/.gemini/antigravity-ide/brain/d6303ffc-29ea-4fc1-9499-09d29b4a7ebc/signal_to_noise_concept.png')
print("Concept diagram saved!")
