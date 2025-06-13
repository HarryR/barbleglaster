#!/usr/bin/env python3
import sys
import sqlite3
import time
import os
from sage.all import GF

def create_table(db_path, bitsize):
    conn = sqlite3.connect(db_path)
    table_name = f"generator_2p{bitsize}_m2p32_mx"
    conn.execute(f"""
        CREATE TABLE IF NOT EXISTS {table_name} (
            mx INTEGER PRIMARY KEY,
            g INTEGER NOT NULL
        )
    """)
    conn.commit()
    return conn, table_name

def get_pending_primes(conn, bitsize):
    """Get mx values that exist in primes table but not in cornacchia table"""
    primes_table = f"primes_2p{bitsize}_m2p32_mx"
    generator_table = f"generator_2p{bitsize}_m2p32_mx"
    trial_table = f"trial_division_2p{bitsize}_m2p32_mx"
    query = f"""
        SELECT p.mx FROM {primes_table} p
        LEFT JOIN {generator_table} c ON p.mx = c.mx
        JOIN {trial_table} td ON td.mx = p.mx AND td.remaining_is_prime = 1
        WHERE c.mx IS NULL
    """
    cursor = conn.execute(query)
    return [row[0] for row in cursor.fetchall()]

def process_generator(bitsize):
    db_path = f"data/{bitsize}.sqlite3"
    if not os.path.exists(db_path):
        print(f"Database {db_path} does not exist. Run the prime finder first.")
        return
    conn, generator_table = create_table(db_path, bitsize)
    pending_mx = get_pending_primes(conn, bitsize)
    if not pending_mx:
        print(f"No pending F_p* generator computations for {bitsize}-bit range")
        conn.close()
        return

    print(f"Processing {len(pending_mx)} pending F_p* generator computations")
    base = 2**bitsize - 2**32
    batch = []
    batch_size = 50000 if bitsize <= 64 else 100
    processed = 0
    sql = f"INSERT OR IGNORE INTO {generator_table} (mx, g) VALUES (?, ?)"
    time_start = time.perf_counter()
    for mx in pending_mx:
        p = base - mx
        F = GF(p)
        g = F.multiplicative_generator()
        batch.append((mx, int(g)))
        processed += 1
        if len(batch) >= batch_size:
            conn.executemany(sql, batch)
            conn.commit()
            time_end = time.perf_counter()
            print(f"Processed {processed} / {len(pending_mx)} - Latest: mx={mx}, g={g}, Time: {time_end-time_start} ({round((time_end-time_start)/len(batch),3)} each)")
            time_start = time_end
            batch = []

    if batch:
        conn.executemany(sql, batch)
        conn.commit()
        print(f"Final batch processed - {len(batch)} items")

    print(f"Generator processing complete: {processed} items processed")
    conn.close()

def main():
    if len(sys.argv) != 2:
        print("Usage: python 4-generator.py <bitsize>")
        print("Example: python 4-generator.py 256")
        sys.exit(1)

    try:
        bitsize = int(sys.argv[1])
        if not (33 <= bitsize <= 512):
            print("Bitsize must be between 33 and 512")
            sys.exit(1)
    except ValueError:
        print("Bitsize must be an integer")
        sys.exit(1)

    process_generator(bitsize)

if __name__ == "__main__":
    main()
