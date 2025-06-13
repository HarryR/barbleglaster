#!/usr/bin/env python3
import sys
import os
import json
import time
import sqlite3
from math import log2
from sage.all import factor
import numpy as np

ORDER_OFFSETS = set([-1,0])

def analyze_factors(x):
    time_start = time.perf_counter()
    factors = factor(x)
    time_end = time.perf_counter()
    #print("  - ", (time_end-time_start), factors)

    factors_primes = [int(prime) for prime, _ in factors]
    factors_powered = [int(prime**power) for prime, power in factors]
    factors_powers = [int(power) for _, power in factors]
    log2_primes = [log2(prime) for prime in factors_primes]
    log2_powered = [log2(powered) for powered in factors_powered]

    # Fixed entropy calculation
    total_powers = sum(factors_powers)
    entropy = -sum((power/total_powers) * log2(power/total_powers)
                   for power in factors_powers if power > 0)

    factors_json = [[str(prime), int(power)] for prime, power in factors]

    return {
        'entropy': entropy,
        'n_factors': len(factors),
        'factors_json': json.dumps(factors_json),

        # Size metrics
        'largest_prime_powered_log2': max(log2_powered),
        'largest_prime_log2': max(log2_primes),
        'smallest_prime_log2': min(log2_primes),
        'smallest_prime_powered_log2': min(log2_powered),

        # Central tendencies
        'avg_prime_powered_log2': np.mean(log2_powered),
        'avg_prime_log2': np.mean(log2_primes),
        'median_prime_powered_log2': np.median(log2_powered),
        'median_prime_log2': np.median(log2_primes),

        # Variance and standard deviation
        'var_prime_powers_log2': np.var(log2_powered),
        'var_prime_log2': np.var(log2_primes),
        'std_prime_powers_log2': np.std(log2_powered),
        'std_prime_log2': np.std(log2_primes),

        # Power distribution
        'max_prime_power': max(factors_powers),
        'total_prime_powers': total_powers,
        'second_largest_prime_log2': sorted(log2_primes, reverse=True)[1] if len(factors) >= 2 else 0,
    }

def create_curvefactor_table(db_path, bitsize):
    conn = sqlite3.connect(db_path)
    curvefactor_table = f"curvefactor_2p{bitsize}_m2p32_mx"
    conn.execute(f"""
        CREATE TABLE IF NOT EXISTS {curvefactor_table} (
            mx INTEGER,
            generator_power INTEGER NOT NULL,
            order_offset INTEGER NOT NULL,
            factors_json TEXT NOT NULL,
            n_factors INTEGER NOT NULL,

            -- Distribution metrics
            entropy REAL NOT NULL,

            -- Size metrics (all log2 values)
            largest_prime_powered_log2 REAL NOT NULL,
            largest_prime_log2 REAL NOT NULL,
            smallest_prime_log2 REAL NOT NULL,
            smallest_prime_powered_log2 REAL NOT NULL,
            second_largest_prime_log2 REAL NOT NULL,

            -- Central tendencies (all log2 values)
            avg_prime_powered_log2 REAL NOT NULL,
            avg_prime_log2 REAL NOT NULL,
            median_prime_powered_log2 REAL NOT NULL,
            median_prime_log2 REAL NOT NULL,

            -- Variance and standard deviation (log2 space)
            var_prime_powers_log2 REAL NOT NULL,
            var_prime_log2 REAL NOT NULL,
            std_prime_powers_log2 REAL NOT NULL,
            std_prime_log2 REAL NOT NULL,

            -- Power distribution (raw integer values)
            max_prime_power INTEGER NOT NULL,
            total_prime_powers INTEGER NOT NULL,

            PRIMARY KEY (mx, generator_power, order_offset)
        )
    """)
    conn.commit()
    return conn, curvefactor_table

def factor_curve_order(a:int, b:int, offset_eisenstein_c:int, offset_eisenstein_d:int, existing_offsets:set[int]):
    c = a + b
    d = 2 * b
    q_c = c + offset_eisenstein_c
    q_d = d + offset_eisenstein_d
    q = q_c**2 + q_d**2 - (q_c * q_d)
    missing_offsets = ORDER_OFFSETS - existing_offsets
    #print(f"Processing mx={mx}, missing: {missing_offsets}")
    for order_offset in missing_offsets:
        #print("Dorp, curve", g_i, mx, q+order_offset, is_prime(q), existing_offsets, row_count)
        yield order_offset, analyze_factors(q+order_offset)

def find_factoring_work(conn, bitsize, batch_size=6):
    cornacchia_table = f"cornacchia_2p{bitsize}_m2p32_mx"
    curves_table = f"curves_2p{bitsize}_m2p32_mx"
    curvefactor_table = f"curvefactor_2p{bitsize}_m2p32_mx"

    # Get prime families that need more processing
    incomplete_families_query = f"""
        SELECT cu.mx, cu.generator_power, cu.offset_eisenstein_c, cu.offset_eisenstein_d,
                cor.a, cor.b,
               COUNT(cft.order_offset) as existing_count
        FROM {curves_table} cu
        LEFT JOIN {curvefactor_table} cft ON  cft.mx = cu.mx AND cft.generator_power = cu.generator_power
        JOIN {cornacchia_table} cor ON cor.mx = cu.mx
        -- Curve family has at least one prime order
        WHERE cu.mx IN (
            SELECT DISTINCT mx
            FROM {curves_table}
            WHERE is_prime = 1
        )
        GROUP BY cu.mx, cu.generator_power
        HAVING existing_count < {len(ORDER_OFFSETS)}
        ORDER BY RANDOM()
        LIMIT {batch_size}
    """

    # Now for each mx that needs work, get the curve details and missing offsets
    for mx, generator_power, offset_eisenstein_c, offset_eisenstein_d, a, b, existing_count in conn.execute(incomplete_families_query).fetchall():
        a = int(a)
        b = int(b)

        # Get specific missing offsets
        existing_offsets_query = f"""
            SELECT DISTINCT order_offset
            FROM {curvefactor_table}
            WHERE mx = ? AND generator_power = ?
        """
        existing_offsets = {row[0] for row in conn.execute(existing_offsets_query, (mx, generator_power)).fetchall()}
        for order_offset, result in factor_curve_order(a, b, offset_eisenstein_c, offset_eisenstein_d, existing_offsets):
            yield mx, generator_power, order_offset, result

def process_curves(bitsize):
    db_path = f"data/{bitsize}.sqlite3"
    if not os.path.exists(db_path):
        print(f"Database {db_path} does not exist. Run previous steps first.")
        return

    conn, curvefactor_table = create_curvefactor_table(db_path, bitsize)
    sql = f"""INSERT OR IGNORE INTO {curvefactor_table}
        (mx, generator_power, order_offset, factors_json, n_factors,
         entropy, largest_prime_powered_log2, largest_prime_log2,
         smallest_prime_log2, smallest_prime_powered_log2, second_largest_prime_log2,
         avg_prime_powered_log2, avg_prime_log2, median_prime_powered_log2,
         median_prime_log2, var_prime_powers_log2, var_prime_log2,
         std_prime_powers_log2, std_prime_log2, max_prime_power, total_prime_powers)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)"""

    while True:
        batch = []
        batch_size = 1200*3 if bitsize <= 64 else 6
        for mx, generator_power, order_offset, result in find_factoring_work(conn, bitsize, batch_size):
            batch.append([
                mx, generator_power, order_offset,
                result['factors_json'], int(result['n_factors']),
                result['entropy'], result['largest_prime_powered_log2'], result['largest_prime_log2'],
                result['smallest_prime_log2'], result['smallest_prime_powered_log2'], result['second_largest_prime_log2'],
                result['avg_prime_powered_log2'], result['avg_prime_log2'], result['median_prime_powered_log2'],
                result['median_prime_log2'], result['var_prime_powers_log2'], result['var_prime_log2'],
                result['std_prime_powers_log2'], result['std_prime_log2'], int(result['max_prime_power']), int(result['total_prime_powers'])
            ])
        if len(batch) > 0:
            conn.executemany(sql, batch)
            conn.commit()
            print(f"Processed {len(batch)}")
            batch = []
        else:
            break
    conn.close()

def main():
    if len(sys.argv) < 2:
        print("Usage: python 7-curvefactor.py <bitsize>")
        print("Example: python 7-curvefactor.py 256")
        sys.exit(1)

    bitsize = int(sys.argv[1])
    if not (33 <= bitsize <= 512):
        print("Bitsize must be between 33 and 512")
        sys.exit(2)

    process_curves(bitsize)

if __name__ == "__main__":
    main()
