#!/usr/bin/env python3
import sys
import os
import gmpy2
import sqlite3
from sage.all import is_prime

gmpy2.get_context().precision = 256

ZETA_gmpy2 = lambda n,x,p: pow(x%p, ((p-1)//n), p)

def curve_orders_eisenstein_offsets():
    return [(-1,0), (-1,-1), (0,-1), (1,0), (1,1), (0,1)]

def curve_orders_eisenstein_coords(c,d):
    return [(c+i, d+j) for i,j in curve_orders_eisenstein_offsets()]

def make_terms_cd(c,d):
    return [(x_0*c) + (x_1*d) for x_0, x_1 in [
        (-2,1), (-1,-1), (1,-2), (2,-1), (1,1), (-1,2)
    ]]

def make_norms_cd(c,d,p):
    return [p + 1 + _ for _ in make_terms_cd(c,d)]

ALL_SUB_PATTERNS = [
    (0, 1, 2, 3, 4, 5),
    (0, 5, 4, 3, 2, 1),
    (3, 4, 5, 0, 1, 2),
    (3, 2, 1, 0, 5, 4),
]

def calculate_curve_orders(p, g, a, b):
    assert p % 12 == 7
    g = gmpy2.mpz(g)
    p = gmpy2.mpz(p)
    a = gmpy2.mpz(a)
    b = gmpy2.mpz(b)
    assert a%2 == 0 and b%2 == 1
    c, d = a + b, 2 * b
    assert c**2 - c*d + d**2 == p
    assert p + 1 + a - (3*b) == p + 1 + c - (2*d)
    u0 = int(((c+d) % 3) == 2)
    u1 = int(gmpy2.mod(gmpy2.fma(ZETA_gmpy2(3,g,p), c, d), p) == 0)
    idx = (u0 * 2) + u1
    result = [0] * 6
    result_eisenstein_offsets = [0] * 7
    norms_cd = make_norms_cd(c,d,p)
    orders_eisenstein_offsets = curve_orders_eisenstein_offsets()
    for i,j in enumerate(ALL_SUB_PATTERNS[idx]):
        result[i] = norms_cd[j]
        result_eisenstein_offsets[i] = orders_eisenstein_offsets[j]
        assert result_eisenstein_offsets[i] == result[i]
    return result, result_eisenstein_offsets

def generate_curves_from_cornacchia(mx, g, a, b, base):
    # a, b: Cornacchia results where p = a² + 3b²
    # base: 2^bitsize - 2^32 (for calculating p = base - mx)
    p = base - mx
    orders, orders_eisenstein_offsets = calculate_curve_orders(p, g, a, b)
    curves = []
    for i in range(6):
        offset_eisenstein_c, offset_eisenstein_d = orders_eisenstein_offsets[i]
        q = int(orders[i])
        curves.append({
            'mx': mx,
            'generator_power': i,
            'is_prime': is_prime(q),
            'offset_eisenstein_c': offset_eisenstein_c,
            'offset_eisenstein_d': offset_eisenstein_d
        })
    return curves

def create_curves_table(db_path, bitsize):
    """Create curves table if it doesn't exist"""
    conn = sqlite3.connect(db_path)
    table_name = f"curves_2p{bitsize}_m2p32_mx"

    conn.execute(f"""
        CREATE TABLE IF NOT EXISTS {table_name} (
            mx INTEGER,
            generator_power INTEGER NOT NULL,
            is_prime INTEGER NOT NULL,
            offset_eisenstein_c INTEGER NOT NULL,
            offset_eisenstein_d INTEGER NOT NULL,
            PRIMARY KEY (mx, generator_power)
        )
    """)
    conn.commit()
    return conn, table_name

def get_pending_curves(conn, bitsize):
    cornacchia_table = f"cornacchia_2p{bitsize}_m2p32_mx"
    curves_table = f"curves_2p{bitsize}_m2p32_mx"
    trialdiv_table = f"trial_division_2p{bitsize}_m2p32_mx"
    generator_table = f"generator_2p{bitsize}_m2p32_mx"

    query = f"""
        SELECT c.mx, c.a, c.b, gt.g
          FROM {cornacchia_table} c
          JOIN {trialdiv_table} td ON td.mx = c.mx
          JOIN {generator_table} gt ON gt.mx = c.mx
          LEFT JOIN {curves_table} cr ON c.mx = cr.mx
         WHERE cr.mx IS NULL
           AND td.remaining_is_prime = 1
    """

    cursor = conn.execute(query)
    return [(int(row[0]), int(row[1]), int(row[2]), int(row[3])) for row in cursor.fetchall()]

def process_curves(bitsize):
    db_path = f"data/{bitsize}.sqlite3"
    if not os.path.exists(db_path):
        print(f"Database {db_path} does not exist. Run previous steps first.")
        return

    conn, curves_table = create_curves_table(db_path, bitsize)
    pending_curves = get_pending_curves(conn, bitsize)
    if not pending_curves:
        print(f"No pending curve computations for {bitsize}-bit range")
        conn.close()
        return

    print(f"Processing {len(pending_curves)} pending curve computations")

    base = 2**bitsize - 2**32
    batch = []
    batch_size = 200000 if bitsize <= 64 else 1000
    processed = 0
    total_curves = 0
    prime_order_curves = 0

    sql = f"""INSERT OR IGNORE INTO {curves_table}
            (mx, generator_power, is_prime, offset_eisenstein_c, offset_eisenstein_d)
            VALUES
            (?,  ?,               ?,        ?,                   ?)
            """
    for mx, a, b, g in pending_curves:
        for curve in generate_curves_from_cornacchia(mx, g, a, b, base):
            batch.append((
                mx,
                curve['generator_power'],
                curve['is_prime'],
                curve['offset_eisenstein_c'],
                curve['offset_eisenstein_d']
            ))
            total_curves += 1
            if curve['is_prime']:
                prime_order_curves += 1
        processed += 1
        if len(batch) >= batch_size:
            # Batch insert
            conn.executemany(sql, batch)
            conn.commit()
            print(f"Processed {processed} / {len(pending_curves)} - Total curves: {total_curves}, Prime order: {prime_order_curves}")
            batch = []

    # Insert remaining batch
    if batch:
        conn.executemany(sql, batch)
        conn.commit()
        print(f"Final batch processed - {len(batch)} items")

    print(f"Curve processing complete:")
    print(f"  'nice' primes processed: {processed}")
    print(f"  Total curves generated: {total_curves}")
    print(f"  Prime order curves: {prime_order_curves}")

    conn.close()

def main():
    if len(sys.argv) != 2:
        print("Usage: python 5-curves.py <bitsize>")
        print("Example: python 5-curves.py 256")
        sys.exit(1)

    bitsize = int(sys.argv[1])
    if not (33 <= bitsize <= 512):
        print("Bitsize must be between 33 and 512")
        sys.exit(1)

    process_curves(bitsize)

if __name__ == "__main__":
    main()
