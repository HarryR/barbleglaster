import os
import sys
import json
import sqlite3
from sage.all import factor, Zmod, is_prime, prod, gcd, Integer, power_mod
from math import log2

def db_open(bitsize):
    db_path = f"data/{bitsize}.sqlite3"
    if not os.path.exists(db_path):
        print(f"Database {db_path} does not exist. Run previous steps first.")
        return
    conn = sqlite3.connect(db_path)
    curvefactor_table = f"curvefactor_2p{bitsize}_m2p32_mx"
    return conn, curvefactor_table

def show_factored_curves(bitsize):
    conn, curvefactor_table = db_open(bitsize)
    cornacchia_table = f"cornacchia_2p{bitsize}_m2p32_mx"
    curves_table = f"curves_2p{bitsize}_m2p32_mx"

    random_mx_sql = f"""
    SELECT cft.mx, COUNT(cft.order_offset) AS offset_count
      FROM {curvefactor_table} cft
      GROUP BY cft.mx
      HAVING offset_count = 18
      ORDER BY RANDOM()
      LIMIT 1
    """
    random_mx = conn.execute(random_mx_sql).fetchone()[0]
    print(random_mx)

    sql = f"""
    SELECT cft.factors_json, ct.generator_power, ct.offset_eisenstein_c, ct.offset_eisenstein_d,
           cor.a, cor.b, cft.order_offset
    FROM {curvefactor_table} cft
    JOIN {curves_table} ct ON ct.mx = cft.mx
    JOIN {cornacchia_table} cor ON cor.mx = cft.mx
    WHERE cft.mx = {random_mx}
      AND cft.generator_power = ct.generator_power
      AND cft.order_offset = 0
    GROUP BY cft.mx, ct.generator_power, ct.offset_eisenstein_c, ct.offset_eisenstein_d, cft.order_offset
    """
    p = 2**bitsize - 2**32 - random_mx
    print('p = 2^256 - 2^32 -', random_mx, ' =', factor(p-1))
    print()
    for factors_json, g_i, offset_c, offset_d, a, b, order_offset in conn.execute(sql).fetchall():
        a, b = int(a), int(b)
        #p = a**2 + 3 * b**2
        q_c, q_d = a + b + offset_c, (2 * b) + offset_d
        q = q_c**2 + q_d**2 - (q_c * q_d)
        factors_json = [(int(prime),int(power)) for prime,power in json.loads(factors_json)]

        print('q', q, '+', order_offset)
        print(f'y^2 = x^3 + g^{g_i}')
        print('\tis_prime(q) ?', is_prime(q))
        #print("\tp - q =", p-q)
        #print("\tp % q =", p % q)
        #print('\tfactor(p%q) =', factor(p%q))
        #print("\tq % p =", q % p)
        #print('\tfactor(q%p) =', factor(q%p))

        #print("\tp > q ?", p > q)
        print("\tfactor(q) =", ' * '.join([''.join([str(prime), f'^{power}' if power != 1 else '']) for prime, power in factors_json]))
        if order_offset == 0:
            #ed_pq = Zmod(p)(q).multiplicative_order()
            #print("\tq^k mod p = 1", '| log2(k) =', log2(ed_pq), '|', factor(ed_pq))
            ed_qp = Zmod(q)(p).multiplicative_order()
            print("\tp^k mod q = 1", '| log2(k) =', round(log2(ed_qp),2))
            print("\t        k =", factor(ed_qp))
        #print("\tfactor(p-q) =", factor(p - q))
        print()

def main():
    if len(sys.argv) < 2:
        print("Usage: python 8-embedding.py <bitsize>")
        print("Example: python 8-embedding.py 256")
        sys.exit(1)

    bitsize = int(sys.argv[1])
    if not (64 <= bitsize <= 512):
        print("Bitsize must be between 64 and 512")
        sys.exit(2)

    show_factored_curves(bitsize)

if __name__ == "__main__":
    main()
