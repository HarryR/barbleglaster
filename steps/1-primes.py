#!/usr/bin/env python3
import os
import sys
import time
import sqlite3
from sage.all import next_prime

def create_database_and_table(db_path, bitsize):
    """Create database and table if they don't exist"""
    os.makedirs(os.path.dirname(db_path), exist_ok=True)

    conn = sqlite3.connect(db_path)
    table_name = f"primes_2p{bitsize}_m2p32_mx"

    conn.execute(f"""
        CREATE TABLE IF NOT EXISTS {table_name} (
            mx INTEGER PRIMARY KEY
        )
    """)
    conn.commit()
    return conn, table_name

def get_resume_point(conn, table_name):
    """Get the smallest mx value to resume from, or 2^31 if table is empty"""
    cursor = conn.execute(f"SELECT MIN(mx) FROM {table_name}")
    result = cursor.fetchone()[0]
    return (result - 1) if result is not None else 2**31

def find_primes_mod_7_12(bitsize):
    """Main function to find primes p ≡ 7 (mod 12) in the specified range"""

    # Database setup
    db_path = f"data/{bitsize}.sqlite3"
    conn, table_name = create_database_and_table(db_path, bitsize)

    # Resume point - start from the largest mx (smallest prime)
    current_mx = get_resume_point(conn, table_name)
    min_mx = 1  # We stop when mx reaches 1

    if current_mx < min_mx:
        print(f"Work already complete for {bitsize}-bit range")
        conn.close()
        return

    print(f"Starting from mx = {current_mx}, searching for {bitsize}-bit primes p ≡ 7 (mod 12)")

    # Starting prime: 2^bitsize - 2^32 - current_mx
    base = 2**bitsize - 2**32
    current_prime = base - current_mx

    batch = []
    batch_size = 100 if bitsize > 64 else 100000

    time_start = time.perf_counter()
    while current_mx >= min_mx:
        # Find next prime (searching upwards)
        p = next_prime(current_prime)

        # Calculate mx for this prime
        prime_mx = int(base - p)

        # If mx has gone below our minimum, we're done
        if prime_mx < min_mx:
            break

        # Check if p ≡ 7 (mod 12)
        if p % 12 == 7:
            batch.append((prime_mx,))

            if len(batch) >= batch_size:
                # Batch insert
                time_db = time.perf_counter()
                conn.executemany(f"INSERT OR IGNORE INTO {table_name} (mx) VALUES (?)", batch)
                conn.commit()
                time_end = time.perf_counter()
                print(f"Inserted batch ending at mx = {prime_mx}, prime = {p} took {time_end-time_start}s ({(time_end-time_start)/len(batch)} each), db took {time_end-time_db}")
                time_start = time_end
                print(f"")
                batch = []

        # Skip ahead by 12 to next candidate that could be ≡ 7 (mod 12)
        current_prime = p + 12
        current_mx = base - current_prime

    # Insert remaining batch
    if batch:
        conn.executemany(f"INSERT OR IGNORE INTO {table_name} (mx) VALUES (?)", batch)
        conn.commit()
        print(f"Final batch inserted, {len(batch)} primes")

    print(f"Search complete for {bitsize}-bit range")
    conn.close()

def main():
    if len(sys.argv) != 2:
        print("Usage: python 1-primes.py <bitsize>")
        print("Example: python 1-primes.py 256")
        sys.exit(1)
    try:
        bitsize = int(sys.argv[1])
        if not (33 <= bitsize <= 512):
            print("Bitsize must be between 33 and 512")
            sys.exit(1)
    except ValueError:
        print("Bitsize must be an integer")
        sys.exit(1)
    find_primes_mod_7_12(bitsize)

if __name__ == "__main__":
    main()
