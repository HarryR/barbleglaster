#!/usr/bin/env python3
import sys
import os
import gmpy2
import sqlite3
from sage.all import GF, EllipticCurve

gmpy2.get_context().precision = 256

MODSQRT_gmpy2 = lambda n, p: gmpy2.powmod(n, (p + 1)//4, p)
ZETA_gmpy2 = lambda n,x,p: gmpy2.powmod(x, ((p-1)//n), p)

def check_glv_endomorphism(curve, p, n, generator):
    seen_betas = set()
    for beta_i in range(2,1000):
        beta_val = ZETA_gmpy2(3, beta_i, p)
        if beta_val == 1 or beta_val in seen_betas:
            continue
        seen_betas.add(beta_val)
        seen_lambdas = set()
        for lambda_i in range(2,1000):
            lambda_val = ZETA_gmpy2(3, lambda_i, n)
            if lambda_val == 1 or lambda_val in seen_lambdas:
                continue
            seen_lambdas.add(lambda_val)
            endo_point = curve(beta_val * generator[0], generator[1])
            scalar_point = lambda_val * generator
            if endo_point == scalar_point:
                return (beta_i, beta_val, lambda_i, lambda_val)
    assert False

def create_glv_table(db_path, bitsize):
    conn = sqlite3.connect(db_path)
    table_name = f"glv_2p{bitsize}_m2p32_mx"
    conn.execute(f"""
        CREATE TABLE IF NOT EXISTS {table_name} (
            mx INTEGER,
            generator_power INTEGER NOT NULL,
            beta_val TEXT,
            lambda_val TEXT,
            beta_i INTEGER,
            lambda_i INTEGER,
            PRIMARY KEY (mx, generator_power)
        )
    """)
    conn.commit()
    return conn, table_name

def get_pending_curves(conn, bitsize):
    cornacchia_table = f"cornacchia_2p{bitsize}_m2p32_mx"
    curves_table = f"curves_2p{bitsize}_m2p32_mx"
    generator_table = f"generator_2p{bitsize}_m2p32_mx"
    glv_table = f"glv_2p{bitsize}_m2p32_mx"

    query = f"""
        SELECT c.mx, c.a, c.b, gt.g,
               cu.generator_power, cu.offset_eisenstein_c, cu.offset_eisenstein_d
          FROM {curves_table} cu
          JOIN {cornacchia_table} c ON c.mx = cu.mx
          JOIN {generator_table} gt ON gt.mx = c.mx
          LEFT JOIN {glv_table} glv ON c.mx = glv.mx
         WHERE cu.is_prime = 1
           AND glv.mx IS NULL
    """

    cursor = conn.execute(query)
    return [(int(row[0]), int(row[1]), int(row[2]), int(row[3]),
             int(row[4]), int(row[5]), int(row[6]))
             for row in cursor.fetchall()]

def find_generator(g,p,E):
    p = gmpy2.mpz(p)
    g = gmpy2.mpz(int(g))
    x = gmpy2.mpz(1)
    while True:
        yy = (gmpy2.powmod(x,3,p) + g) % p
        y = MODSQRT_gmpy2(yy, p)
        if (y*y) % p == yy:
            if y & 1:
                y = p - y
            if E.point((x,y)).order() == E.order():
                return x,y
        x += 1

def process_curves(bitsize):
    db_path = f"data/{bitsize}.sqlite3"
    if not os.path.exists(db_path):
        print(f"Database {db_path} does not exist. Run previous steps first.")
        return

    conn, glv_table = create_glv_table(db_path, bitsize)
    pending_curves = get_pending_curves(conn, bitsize)
    if not pending_curves:
        print(f"No pending curve computations for {bitsize}-bit range")
        conn.close()
        return

    print(f"Processing {len(pending_curves)} pending curve computations")

    batch = []
    batch_size = 1000 if bitsize <= 64 else 10
    processed = 0
    total_curves = 0
    glv_curves = 0

    sql = f"""INSERT OR IGNORE INTO {glv_table}
            (mx, generator_power, beta_val, lambda_val, beta_i, lambda_i)
            VALUES
            (?,  ?,               ?,        ?,          ?,      ?)
            """
    for mx, a, b, g, g_i, off_c, off_d in pending_curves:
        p = (a**2) + (3*(b**2))
        Fp = GF(p)
        E = EllipticCurve(Fp, [0, g**g_i])
        G = E.point(find_generator(g**g_i, p, E))
        c = a + b
        d = 2 * b
        q_c = c + off_c
        q_d = d + off_d
        q = q_c**2 + q_d**2 - (q_c * q_d)
        beta_i, beta_val, lambda_i, lambda_val = check_glv_endomorphism(E, p, q, G)
        batch.append((
            mx,
            g_i,
            str(int(beta_val)),
            str(int(lambda_val)),
            int(beta_i),
            int(lambda_i)
        ))
        total_curves += 1
        if beta_i is not None:
            glv_curves += 1
        processed += 1
        if len(batch) >= batch_size:
            conn.executemany(sql, batch)
            conn.commit()
            print(f"Processed {processed} / {len(pending_curves)} - Total curves: {total_curves}, GLV?: {glv_curves}")
            batch = []

    if batch:
        conn.executemany(sql, batch)
        conn.commit()
        print(f"Final batch processed - {len(batch)} items")

    print(f"Curve processing complete:")
    print(f"  'nice' primes processed: {processed}")
    print(f"  Total curves generated: {total_curves}")
    print(f"  Curves with GLV: {glv_curves}")

    conn.close()

def main():
    if len(sys.argv) != 2:
        print("Usage: python 6-glv.py <bitsize>")
        print("Example: python 6-glv.py 256")
        sys.exit(1)

    bitsize = int(sys.argv[1])
    if not (33 <= bitsize <= 512):
        print("Bitsize must be between 33 and 512")
        sys.exit(1)

    process_curves(bitsize)

if __name__ == "__main__":
    main()
