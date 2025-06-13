#!/usr/bin/env python3
import os
import sys
import json
import sqlite3
from math import log2
from sage.all import is_prime
from sage.rings.factorint import factor_trial_division

def create_trial_division_table(conn, bitsize):
    """Create trial division table if it doesn't exist"""
    table_name = f"trial_division_2p{bitsize}_m2p32_mx"
    conn.execute(f"""
        CREATE TABLE IF NOT EXISTS {table_name} (
            mx INTEGER PRIMARY KEY,
            factors_json TEXT NOT NULL,
            remaining_is_prime INTEGER NOT NULL,
            remaining_log2 REAL NOT NULL
        )
    """)
    conn.commit()
    return table_name

def get_pending_trial_division(conn, bitsize):
    """Get mx values that exist in primes table but not in trial division table"""
    primes_table = f"primes_2p{bitsize}_m2p32_mx"
    trial_table = f"trial_division_2p{bitsize}_m2p32_mx"

    query = f"""
        SELECT p.mx FROM {primes_table} p
        LEFT JOIN {trial_table} t ON p.mx = t.mx
        WHERE t.mx IS NULL
    """

    cursor = conn.execute(query)
    return [row[0] for row in cursor.fetchall()]

def analyze_prime_minus_one(p, small_bits=16):
    p_minus_one = p - 1
    trial_factors = factor_trial_division(p_minus_one, 2**small_bits)

    # Convert to list of [prime, power] pairs
    small_factors = []
    remaining = p_minus_one

    for prime_power in trial_factors:
        prime = prime_power[0]
        power = prime_power[1]
        product = prime**power
        product_is_prime = is_prime(product)
        if product_is_prime:
            small_factors.append([str(int(prime)), int(power)])
            remaining //= (prime ** power)
        else:
            break

    # Check if remaining factor is prime (or 1)
    if remaining == 1:
        remaining_is_prime = True  # Vacuously true
    else:
        remaining_is_prime = is_prime(remaining)

    return small_factors, remaining_is_prime, int(remaining)

def process_trial_division(bitsize):
    """Process trial division for pending primes"""

    # Database setup
    db_path = f"data/{bitsize}.sqlite3"
    if not os.path.exists(db_path):
        print(f"Database {db_path} does not exist. Run the prime finder first.")
        return

    conn = sqlite3.connect(db_path)
    trial_table = create_trial_division_table(conn, bitsize)

    # Get pending mx values
    pending_mx = get_pending_trial_division(conn, bitsize)

    if not pending_mx:
        print(f"No pending trial division computations for {bitsize}-bit range")
        conn.close()
        return

    print(f"Processing {len(pending_mx)} pending trial division computations")

    base = 2**bitsize - 2**32
    batch = []
    batch_size = 200000 if bitsize <= 64 else 1000
    trial_div_size = int(bitsize/64)
    processed = 0
    satisfying_condition = 0

    sql = f"INSERT OR IGNORE INTO {trial_table} (mx, factors_json, remaining_is_prime, remaining_log2) VALUES (?, ?, ?, ?)"

    for mx in pending_mx:
        # Calculate prime from mx
        p = base - mx

        # Analyze p-1
        small_factors, remaining_is_prime, remaining = analyze_prime_minus_one(p, trial_div_size)

        # Store results
        factors_json = json.dumps(small_factors)
        remaining_is_prime_int = 1 if remaining_is_prime else 0
        remaining_log2 = log2(remaining)

        record = (mx, factors_json, remaining_is_prime_int, remaining_log2)
        batch.append(record)
        processed += 1

        if remaining_is_prime:
            satisfying_condition += 1

        if len(batch) >= batch_size:
            conn.executemany(sql, batch)
            conn.commit()
            print(f"Processed {processed} / {len(pending_mx)} - Satisfying condition: {satisfying_condition} - Latest: mx={mx}, remaining={remaining}, prime={remaining_is_prime}")
            batch = []

    # Insert remaining batch
    if batch:
        conn.executemany(sql, batch)
        conn.commit()
        print(f"Final batch processed - {len(batch)} items")

    print(f"Trial division processing complete:")
    print(f"  Total processed: {processed}")
    print(f"  Satisfying condition (remaining factor is prime): {satisfying_condition}")
    print(f"  Percentage satisfying: {satisfying_condition/processed*100:.2f}%")

    conn.close()

def query_results(bitsize):
    """Query and display some results"""
    db_path = f"data/{bitsize}.sqlite3"
    conn = sqlite3.connect(db_path)
    trial_table = f"trial_division_2p{bitsize}_m2p32_mx"

    # Get some examples of primes that satisfy the condition
    cursor = conn.execute(f"""
        SELECT mx, factors_json
        FROM {trial_table}
        WHERE remaining_is_prime = 1
        LIMIT 5
    """)

    print("\nSample results (satisfying condition):")
    base = 2**bitsize - 2**32
    for mx, factors_json in cursor:
        p = base - mx
        factors = json.loads(factors_json)
        print(f"  mx={mx}, p={p}")
        print(f"    Small factors: {factors}")
        print()

    conn.close()

def main():
    if len(sys.argv) != 2:
        print("Usage: python trial_division_script.py <bitsize>")
        print("Example: python trial_division_script.py 256")
        sys.exit(1)
    try:
        bitsize = int(sys.argv[1])
        if not (33 <= bitsize <= 512):
            print("Bitsize must be between 33 and 512")
            sys.exit(1)
    except ValueError:
        print("Bitsize must be an integer")
        sys.exit(1)
    process_trial_division(bitsize)
    query_results(bitsize)

if __name__ == "__main__":
    main()
