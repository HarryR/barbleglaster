import numpy as np
import matplotlib.pyplot as plt
from math import sqrt, log, log2
from scipy.optimize import root_scalar, minimize_scalar

log2_p = 256
p = 2**log2_p
log_p = log(p)
sqrt_p = sqrt(p)

# Define the function whose crossings we want to find
def normalized_value(x):
    d = 2**(log2(p)*x)
    a = log(p) * sqrt(p/d) + sqrt(d)
    return log2(a)/log2(p)

# Function to find a crossing (returns zero when function crosses a threshold)
def find_crossing(x, threshold):
    return normalized_value(x) - threshold

# Find precise crossing points using numerical methods
threshold = 0.475
target_value = 0.5

# First 0.5 crossing
result_05 = root_scalar(lambda x: find_crossing(x, target_value),
                       bracket=[0.05, 0.15], method='brentq')
x_05 = result_05.root

# First threshold crossing
result_first = root_scalar(lambda x: find_crossing(x, threshold),
                          bracket=[0.08, 0.15], method='brentq')
x_first = result_first.root

# Second threshold crossing
result_second = root_scalar(lambda x: find_crossing(x, threshold),
                           bracket=[0.9, 0.98], method='brentq')
x_second = result_second.root

# Find the minimum value
result_min = minimize_scalar(normalized_value, bounds=(0.00001, 0.99999), method='bounded')
x_min = result_min.x
y_min = normalized_value(x_min)

# Round values for display
x_05_rounded = round(x_05, 2)
x_first_rounded = round(x_first, 2)
x_second_rounded = round(x_second, 2)
x_min_rounded = round(x_min, 2)

# Create the plot
# Use exact boundaries from 0 to 1
x_values = np.linspace(0.00001, 0.99999, 1000)
y_values = np.array([normalized_value(x) for x in x_values])

plt.figure(figsize=(5, 3))

# Plot the full line in black
#plt.plot(x_values, y_values, 'black', linewidth=2)

# Fill areas above/below threshold
plt.fill_between(x_values, y_values, threshold,
                where=(y_values >= threshold),
                alpha=0.3, color='green', label=f'Good (≥ {threshold})')

plt.fill_between(x_values, y_values, threshold,
                where=(y_values < threshold),
                alpha=0.3, color='red', label=f'Bad (< {threshold})')

# Add horizontal lines at thresholds
plt.axhline(y=threshold, color='black')

# Add horizontal line at minimum y value
#plt.axhline(y=y_min, color='red', linestyle=':', linewidth=1)

# Mark the minimum point with a star
#plt.plot(x_min, y_min, 'r*', markersize=10, label=f'Minimum ({x_min_rounded}, {round(y_min, 4)})')

# Custom x-ticks to include only 0, 1, the crossing points, and the minimum point
custom_ticks = [0.0, x_first_rounded, x_05_rounded, x_min_rounded, x_second_rounded, 1.0]
plt.xticks(custom_ticks, [f"{x:.2f}" for x in custom_ticks])

# Custom y-ticks to include the minimum value
y_min_rounded = round(y_min, 2)
y_ticks = list(plt.yticks()[0])  # Get current y-ticks
if y_min_rounded not in y_ticks:
    y_ticks.append(y_min_rounded)
    y_ticks.sort()
    plt.yticks(y_ticks)

# Set exact x-axis limits
plt.xlim(0, 1)

# Set labels and title
plt.xlabel('x = log2(d)/log2(p)')
#plt.ylabel('y = log2(O(log(p)*sqrt(p/d)+sqrt(d)))/log2(p)')
#plt.title('Functions vs sqrt(p) where p = 2^256')
plt.grid(True)
#plt.legend()
plt.tight_layout()
plt.savefig(f"graphs/cheon_resistance.png", dpi=300)
plt.savefig(f"graphs/cheon_resistance.svg", dpi=300)
#plt.show()

# Print the exact crossing points and minimum for reference
print(f"Highest y value", y_values[0])
print(f"First 0.5 crossing at x = {x_05:.6f} (rounded to {x_05_rounded})")
print(f"First threshold crossing at x = {x_first:.6f} (rounded to {x_first_rounded})")
print(f"Second threshold crossing at x = {x_second:.6f} (rounded to {x_second_rounded})")
print(f"Minimum value occurs at x = {x_min:.6f} (rounded to {x_min_rounded})")
print(f"The minimum y value is {y_min:.6f}")