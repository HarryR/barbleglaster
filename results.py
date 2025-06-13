#!/usr/bin/env python3
import os
import sys
import math
import random
import sqlite3
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
from sage.all import EllipticCurve, GF, is_prime
from lib_glv import Curve256GLV, find_generator, test_curve
from lib_eta import eta, eta_norm, eta_map, factors_load, factors_metrics, factors_metrics_map, minmax, factors_str

def db_open(bitsize) -> sqlite3.Connection:
    # Database setup
    db_path = f"data/{bitsize}.sqlite3"
    if not os.path.exists(db_path):
        print(f"Database {db_path} does not exist. Run the prime finder first.")
        return
    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    return conn

def get_curves_by_mx(conn:sqlite3.Connection, bitsize, mx) -> sqlite3.Row:
    # Retrieve a prime-ordered curve given its
    cornacchia_table = f"cornacchia_2p{bitsize}_m2p32_mx"
    curves_table = f"curves_2p{bitsize}_m2p32_mx"
    generator_table = f"generator_2p{bitsize}_m2p32_mx"
    glv_table = f"glv_2p{bitsize}_m2p32_mx"
    sql = f"""
        SELECT cu.generator_power, cu.offset_eisenstein_c, cu.offset_eisenstein_d,
               cor.a, cor.b,
               gt.g AS prime_gen,
               glv.lambda_val, glv.lambda_i,
               glv.beta_val, glv.beta_i
         FROM {curves_table} cu
         JOIN {cornacchia_table} cor ON cor.mx = cu.mx
         JOIN {generator_table} gt ON gt.mx = cu.mx
         JOIN {glv_table} glv ON glv.mx = cu.mx AND glv.generator_power = cu.generator_power
        WHERE cu.mx = ?
          AND cu.is_prime = 1
         LIMIT 1
    """
    return conn.execute(sql, (mx,)).fetchall()

def twist_factors_by_mx(conn:sqlite3.Connection,bitsize:int,mx):
    curvefactor_table = f"curvefactor_2p{bitsize}_m2p32_mx"
    curves_table = f"curves_2p{bitsize}_m2p32_mx"
    sql = f"""
        SELECT ct.mx, ct.is_prime, ct.generator_power, cft.order_offset, cft.factors_json
        FROM {curvefactor_table} cft
        JOIN {curves_table} ct
             ON ct.mx = cft.mx
        WHERE (cft.order_offset = 0 OR (ct.is_prime AND cft.order_offset == -1))
          AND ct.mx = {mx}
          AND cft.generator_power = ct.generator_power
        GROUP BY ct.mx, ct.generator_power, cft.order_offset, cft.factors_json
        ORDER BY cft.mx, ct.generator_power ASC, order_offset ASC
    """
    return conn.execute(sql).fetchall()

def twist_factors(conn:sqlite3.Connection,bitsize:int):
    curvefactor_table = f"curvefactor_2p{bitsize}_m2p32_mx"
    curves_table = f"curves_2p{bitsize}_m2p32_mx"
    sql = f"""
        SELECT ct.mx, ct.is_prime, cft.generator_power, cft.order_offset, cft.factors_json
        FROM {curvefactor_table} cft
        JOIN {curves_table} ct ON ct.mx = cft.mx
        WHERE cft.generator_power = ct.generator_power
          AND cft.order_offset <= 0
          AND ct.mx IN(SELECT DISTINCT ct2.mx FROM {curves_table} ct2 WHERE ct2.is_prime = 1)
        GROUP BY ct.mx, ct.is_prime, ct.generator_power, cft.order_offset, cft.factors_json
        ORDER BY cft.mx, ct.generator_power ASC, order_offset ASC
    """
    return conn.execute(sql).fetchall()

def curve_metrics(conn:sqlite3.Connection,bitsize:int):
    results = eta_map(eta_norm({
        (mx,is_prime,g_i,order_offset): factors_metrics(factors_load(factors_json), bitsize)
        for mx, is_prime, g_i, order_offset, factors_json in twist_factors(conn,bitsize)
    }))
    # Then split the pime and non-prime orders
    new_results = defaultdict(list)
    for (mx,is_prime,_,_),v in results.items():
        new_results[mx,is_prime] += v
    results = new_results
    new_results = defaultdict(list)
    for (mx,is_prime),v in eta_map(results).items():
        new_results[mx] += v if is_prime else [min(v)]
    return eta_map(new_results)

def rank_primes(conn: sqlite3.Connection, bitsize: int) -> dict[int, dict[str, float]]:
    trial_div_table = f"trial_division_2p{bitsize}_m2p32_mx"
    curves_table = f"curves_2p{bitsize}_m2p32_mx"
    sql = f"""
SELECT td.mx, td.factors_json
 FROM {trial_div_table} td
 JOIN {curves_table} cu ON (cu.mx = td.mx)
 WHERE cu.is_prime = 1
 GROUP BY td.mx
    """
    all_primes = list(conn.execute(sql).fetchall())
    prime_factors = {mx: factors_load(factors_json)
                     for mx, factors_json in all_primes}
    prime_ranks = factors_metrics_map(prime_factors, bitsize)
    return prime_ranks, prime_factors

def show_curve(bitsize, mx, curve:sqlite3.Row, rank, is_interesting):
    p_a, p_b = int(curve['a']), int(curve['b'])
    p = p_a**2 + 3 * p_b**2
    assert is_prime(p)
    F_p = GF(p)
    p_g = F_p(curve['prime_gen'])
    p_c = p_a + p_b
    p_d = 2 * p_b
    q_c = p_c + curve['offset_eisenstein_c']
    q_d = p_d + curve['offset_eisenstein_d']
    q = q_c**2 + q_d**2 - (q_c*q_d)
    g_i = curve['generator_power']
    assert is_prime(q)
    E = EllipticCurve(F_p, [0, p_g**g_i])
    print(f"p = 2^{bitsize} - 2^32 - {mx} = a^2 + 3b^2 = c^2 + d^2 - cd")
    print(f"  =", hex(p))
    #print(f"mx: {mx}")
    print(f"(a,b):", (int(curve['a']), int(curve['b'])))
    #print(f"\tp%i for i in 2..12 = ", [(i, p%i) for i in range(2,13)])
    #print(f"F_p* g:", p_g)
    print(f"curve: y^2 = x^3 + g^{g_i}   note: g={p_g}, g^{g_i} = {p_g**g_i}")
    #print(f"\t CM discriminant:", E.cm_discriminant())
    print(f"|E_p_{g_i}|=q:", hex(q))
    #print(f"\tgcd(p-1,q-1)", math.gcd(p-1,q-1))
    #print(f"\tlog2(lcm(p-1,q-1))", math.log2(math.lcm(p-1,q-1)))
    #print(f"\tq%i for i in 2..12 = ", [(i, q%i) for i in range(2,13)])
    G = find_generator(p_g**g_i, p, E)
    print(f" E_p_{g_i} G: ({hex(G[0])},{hex(G[1])})")
    #print(f"\tGLV endomorphism: lambda * k * G = k * Phi(G) = k * (beta * G.x, G.y)")

    print("embedding degree log2:", round(math.log2(GF(q)(p).multiplicative_order()),2))
    curve:Curve256GLV
    curve, scores = test_curve(p, p_g**g_i)
    #print(f"glv     scores:", scores)
    glv = curve.glv
    print(f"glv     lambda: {glv.lambda_i}^((q-1)/3) =", hex(int(glv.lambda_val)))
    print(f"\t  beta: {glv.beta_i}^((p-1)/3) =", hex(int(glv.beta)))
    print(f"\t   b_1: {glv.b1}")
    print(f"\t   b_2: {glv.b2}")
    print(f"\t   g_1: {glv.g1}")
    print(f"\t   g_2: {glv.g2}")

    #print(f"\t   decomposition score: {scores}")
    k = random.randint(1, curve.n)
    k1,k2 = curve.decompose_scalar(k)
    #print("k1,k2",k1,k2)
    k1_log2 = round(math.log2(abs(k1))) if k1 != 0 else 0
    k2_log2 = round(math.log2(abs(k2))) if k2 != 0 else 0
    print(f"\t   decomposition: log2(k1)={k1_log2} log2(k2)={k2_log2} score={scores[1]}")

    if is_interesting or rank == 0.0 or rank == 1.0:
        # Embedding degrees of each curve order to others
        #for a in []
        pass

def get_scores(conn, bitsize):
    prime_scores, prime_factors = rank_primes(conn, bitsize)
    total_scores = []
    for mx,curve_scores in curve_metrics(conn, bitsize).items():
        normalized_score = eta(curve_scores + prime_scores[mx])
        total_scores.append((mx,normalized_score))
    total_scores = sorted(total_scores, key=lambda _: _[1], reverse=True)
    return total_scores, prime_factors

def process(bitsize, interest:list[int]):
    conn = db_open(bitsize)
    total_scores,prime_factors = get_scores(conn, bitsize)
    tsd = dict(total_scores)
    interest_mx_ids = set(interest)
    interest = [(mx,tsd[mx]) for mx in interest]
    xmin,xmax = minmax([_[1] for _ in total_scores])
    hist = defaultdict(int)
    for mx,score in total_scores:
        hist[int(score*100)] += 1
    for mx,score in total_scores[:3] + total_scores[-3:] + interest:
        rank = (score - xmin) / (xmax - xmin) if xmax != xmin else 0.0
        for curve in get_curves_by_mx(conn, bitsize, mx):
            show_curve(bitsize, mx, curve, rank, mx in interest_mx_ids)
        print("factor(p-1) =", factors_str(prime_factors[mx]))
        for _, is_prime, generator_power, order_offset, factors_json in twist_factors_by_mx(conn,bitsize,mx):
            offset_str = ''
            if order_offset != 0:
                if order_offset < 0:
                    offset_str = f' - {abs(order_offset)}'
                else:
                    offset_str = f' + {order_offset}'
            print(f"twist g^{generator_power}{offset_str} =",
                  factors_str(factors_load(factors_json)),
                  'prime' if (is_prime and order_offset == 0) else '')
        print("score:", score)
        print("rank:", rank)
        print()
    print()
    print(" Score range", xmin, xmax)
    for k in sorted(hist.keys()):
        print(f"\t{k}", hist[k])

    fig, ax = plt.subplots(1, 1, figsize=(5, 3))
    ax.set_xlabel('Score')
    ax.set_ylabel('Count')
    sorted_vals = [hist[k] for k in sorted(hist.keys())]
    ax.bar(sorted(hist.keys()), sorted_vals)
    ax.axhline(np.mean(sorted_vals), label="mean", linewidth=2, linestyle=':')
    ax.axhline(np.median(sorted_vals), label="median", color='k', alpha=0.5, linewidth=1, linestyle='--')
    plt.tight_layout()
    plt.legend()
    filename = "graphs/results" if bitsize == 256 else f"graphs/results-{bitsize}"
    plt.savefig(f"{filename}.png", dpi=300)
    plt.savefig(f"{filename}.svg", dpi=300)

def main():
    if len(sys.argv) < 2:
        print("Usage: python results.py <bitsize>")
        print("Example: python results.py 256")
        sys.exit(1)

    try:
        bitsize = int(sys.argv[1])
        if not (33 <= bitsize <= 512):
            print("Bitsize must be between 33 and 512")
            sys.exit(1)
    except ValueError:
        print("Bitsize must be an integer")
        sys.exit(1)

    interests = []
    if len(sys.argv) > 2:
        interests = [int(_) for _ in sys.argv[2:]]

    process(bitsize, interests)

if __name__ == "__main__":
    main()
