"""
Graphs performance of our curve order determination algorithm
vs Sage's point counting algorithms,

use lemma/*-eisenstein-mapping* to generate the raw data file
"""

import matplotlib.pyplot as plt
import pandas as pd

# Load the data
df = pd.read_csv('data/2-vs-sage.csv')

# Calculate the ratio
ratio = df['ours'] / df['sage']

# Create the plot
plt.figure(figsize=(10, 6))
plt.plot(df['bitlen'], ratio)

plt.xlabel('Bit Length')
plt.ylabel('Ratio (Ours/Sage)')
plt.title('Performance Ratio: Ours vs Sage')
plt.axhline(y=1, color='red', linestyle='--', alpha=0.5, label='Equal Performance')
plt.yscale('log')
plt.legend()
plt.grid(True, alpha=0.3)

plt.tight_layout()

#plt.show()
plt.savefig("graphs/deterministic-curve-order-performance.png", dpi=300)
plt.savefig("graphs/deterministic-curve-order-performance.svg", dpi=300)
