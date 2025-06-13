#!/usr/bin/env python3
import sys
import sqlite3
import os
from math import isqrt
import gmpy2

gmpy2.get_context().precision = 256

ZETA_gmpy2 = lambda n,x,p: pow(x%p, ((p-1)//n), p)
MODSQRT_gmpy2 = lambda n, p: pow(n%p, (p + 1) // 4, p)

def cornacchia_gmpy2(d, p):
    """Standard Cornacchia algorithm for x^2 + d*y^2 = p"""
    assert p % 12 == 7
    p = gmpy2.mpz(p)
    if d <= 0 or d >= p:
        raise ValueError("invalid input")
    if ZETA_gmpy2(2, -d, p) != 1:
        raise ValueError("no solution")
    x0 = MODSQRT_gmpy2(-d, p)
    # Choose the larger square root
    if x0 < p // 2:
        x0 = p - x0
    # Extended Euclidean algorithm
    a, b = p, x0
    limit = isqrt(p) # Python 3.8+
    while b > limit:
        a, b = b, a % b
        assert a > 0 and b > 0 # guaranteed positive integers
    remainder = p - b * b
    assert remainder % d == 0  # guaranteed congruence
    c = remainder // d
    t = isqrt(c) # Python 3.8+
    assert t * t == c  # guaranteed exact squares
    return b, t

def create_cornacchia_table(db_path, bitsize):
    """Create cornacchia table if it doesn't exist"""
    conn = sqlite3.connect(db_path)
    table_name = f"cornacchia_2p{bitsize}_m2p32_mx"

    conn.execute(f"""
        CREATE TABLE IF NOT EXISTS {table_name} (
            mx INTEGER PRIMARY KEY,
            a TEXT NOT NULL,
            b TEXT NOT NULL
        )
    """)
    conn.commit()
    return conn, table_name

def get_pending_primes(conn, bitsize):
    """Get mx values that exist in primes table but not in cornacchia table"""
    primes_table = f"primes_2p{bitsize}_m2p32_mx"
    cornacchia_table = f"cornacchia_2p{bitsize}_m2p32_mx"

    query = f"""
        SELECT p.mx FROM {primes_table} p
        LEFT JOIN {cornacchia_table} c ON p.mx = c.mx
        WHERE c.mx IS NULL
    """

    cursor = conn.execute(query)
    return [row[0] for row in cursor.fetchall()]

def process_cornacchia(bitsize):
    """Process Cornacchia algorithm for pending primes"""

    # Database setup
    db_path = f"data/{bitsize}.sqlite3"
    if not os.path.exists(db_path):
        print(f"Database {db_path} does not exist. Run the prime finder first.")
        return

    conn, cornacchia_table = create_cornacchia_table(db_path, bitsize)

    # Get pending mx values
    pending_mx = get_pending_primes(conn, bitsize)

    if not pending_mx:
        print(f"No pending Cornacchia computations for {bitsize}-bit range")
        conn.close()
        return

    print(f"Processing {len(pending_mx)} pending Cornacchia computations")

    base = 2**bitsize - 2**32
    batch = []
    batch_size = 1000 if bitsize > 64 else 200000
    processed = 0

    for mx in pending_mx:
        # Calculate prime from mx
        p = base - mx

        # Apply Cornacchia with d=3
        a, b = cornacchia_gmpy2(3, p)

        batch.append((mx, str(int(a)), str(int(b))))
        processed += 1

        if len(batch) >= batch_size:
            # Batch insert
            conn.executemany(f"INSERT OR IGNORE INTO {cornacchia_table} (mx, a, b) VALUES (?, ?, ?)", batch)
            conn.commit()
            print(f"Processed {processed} / {len(pending_mx)} - Latest: mx={mx}, a={a}, b={b}")
            batch = []

    # Insert remaining batch
    if batch:
        conn.executemany(f"INSERT OR IGNORE INTO {cornacchia_table} (mx, a, b) VALUES (?, ?, ?)", batch)
        conn.commit()
        print(f"Final batch processed - {len(batch)} items")

    print(f"Cornacchia processing complete: {processed} items processed")
    conn.close()

def main():
    if len(sys.argv) != 2:
        print("Usage: python 2-cornacchia.py <bitsize>")
        print("Example: python 2-cornacchia.py 256")
        sys.exit(1)

    try:
        bitsize = int(sys.argv[1])
        if not (33 <= bitsize <= 512):
            print("Bitsize must be between 33 and 512")
            sys.exit(1)
    except ValueError:
        print("Bitsize must be an integer")
        sys.exit(1)

    process_cornacchia(bitsize)

if __name__ == "__main__":
    main()
