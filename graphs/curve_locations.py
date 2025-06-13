import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import sqlite3
import math
import os

def db_open(bitsize):
    db_path = f"data/{bitsize}.sqlite3"
    if not os.path.exists(db_path):
        print(f"Database {db_path} does not exist. Run previous steps first.")
        return None, None
    conn = sqlite3.connect(db_path)
    curves_table = f"curves_2p{bitsize}_m2p32_mx"
    return conn, curves_table

def log2_formatter(x, pos):
    if x <= 0:
        return '0'
    log_val = np.log2(x)
    if log_val == int(log_val):
        return f'2^{int(log_val)}'
    else:
        return f'2^{int(log_val)}'

def plot_cornacchia_data(bitsize):
    # Open database connection
    conn, curves_table = db_open(bitsize)
    if conn is None:
        return

    cornacchia_table = f"cornacchia_2p{bitsize}_m2p32_mx"

    # Query the data
    query = f"""
    SELECT cor.a, cor.b
    FROM {curves_table} ct
    JOIN {cornacchia_table} cor ON cor.mx = ct.mx
    WHERE ct.is_prime = 1
    ORDER BY RANDOM()
    LIMIT 10000
    """

    # Execute query and fetch data
    cursor = conn.execute(query)
    data = cursor.fetchall()
    conn.close()
    if not data:
        print("No data found with the given query.")
        return

    # Convert to numpy arrays
    #from functools import reduce
    #from operator import mul
    #i = reduce(mul, primes_first_n(25))

    a_values = np.array([int(row[0]) for row in data]) # * 3 * 29
    b_values = np.array([int(row[1]) for row in data]) # * 2 * 3 * 5**2

    print(f"Retrieved {len(data)} data points")

    # Create side-by-side plots
    fig = plt.figure(figsize=(5,3))
    ax1 = fig.add_subplot()
    ma,mb = max(a_values), max(b_values)
    ax1.axhline(max(b_values), linewidth=1, alpha=0.2, linestyle='--', color='k')
    ax1.axvline(max(a_values), linewidth=1, alpha=0.2, linestyle='--', color='k')

    print("max a =", ma)
    print("max b =", mb)

    plt.grid(True, alpha=0.3)

    #"""
    ax1.scatter(
        a_values % b_values,
        b_values % a_values,
        marker='.', s=0.5,
        label="x=a%b y=b%a")
    #"""

    """
    hb = ax1.hexbin(
        [math.dist([_[0]], [_[1]]) for _ in zip(a_values, b_values)],
        [math.isqrt(a*b) for a,b in zip(a_values, b_values)],
        bins='log',
    )
    """

    """
    ax1.scatter(
        [abs(a-b) for a,b in zip(a_values, b_values)],
        [math.sqrt(b*a) for a,b in zip(a_values, b_values)],
        marker='.', s=0.5,
        label="x=|b-a| y=sqrt(ab)")
    #"""

    #ax1.hexbin(a_values, b_values, bins='log')
    ax1.scatter(a_values, b_values, s=0.05, label="p = a^2 + 3b^2")

    ax1.set_title(f'Cornacchia p%12=7 {bitsize}-bit GLV curves')
    ax1.set_xlabel('log2(a)')
    ax1.set_ylabel('log2(b)')
    #plt.annotate(f'{round(math.log2(ma),1)}', (ma,mb))

    plt.legend()

    #"""
    ax1.set_yscale('log')
    ax1.set_xscale('log')
    ax1.set_xticklabels([str(int(math.log2(x))) for x in ax1.get_xticks().tolist()])
    ax1.set_yticklabels([str(int(math.log2(y))) for y in ax1.get_yticks().tolist()])
    #"""

    #ax1.xaxis.set_major_formatter(ticker.FuncFormatter(log2_formatter))
    #ax1.yaxis.set_major_formatter(ticker.FuncFormatter(log2_formatter))

    plt.tight_layout()
    plt.savefig("graphs/curve_locations.png", dpi=300)
    plt.savefig("graphs/curve_locations.svg", dpi=300)
    #plt.show()

if __name__ == "__main__":
    plot_cornacchia_data(256)