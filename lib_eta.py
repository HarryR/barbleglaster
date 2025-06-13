# General routines for taking averages and normalising values (ranking)

import math
import json

def avg(x):
    return sum(x) / len(x)

def dist(x):
    return math.sqrt((1/len(x)) * sum(_**2 for _ in x))

def eta(x:list[float]):
    return (dist(x) + avg(x)) / 2
    #return dist(x)
    #return avg(x)

def eta_map(x:dict[int,list[float]]):
    return {k: [eta(v)] for k,v in x.items()}

def minmax(a):
    return min(a), max(a)

def normalize_value(value, xmin, xmax):
    return (value - xmin) / (xmax - xmin) if xmax != xmin else 0.0

def eta_norm(x:dict[int,list[float]]):
    keys = list(x.keys())
    n_stats = len(x[keys[0]])
    minmaxes = [minmax([x[k][i] for k in keys]) for i in range(n_stats)]
    fn = lambda facets: [normalize_value(value, xmin, xmax)
                        for (xmin, xmax), value in zip(minmaxes,facets)]
    result = {k: fn(v) for k,v in x.items()}
    for k in result.keys():
        result[k] = [1-f if f < 0 else f for f in result[k]]
    return result

def factors_str(factors:list[tuple[int,int]]):
    return ' * '.join([f'{prime if prime < 1000 else hex(prime)}' + ('' if power == 1 else f'^{power}') for prime,power in factors])

def factors_metrics(factors:list[tuple[int,int]],bitsize:int):
    return [
            (math.log2(factors[-1][0]) * factors[-1][1]) / bitsize,
            #-max([math.log2(prime**power) for prime, power in factors[:-1]]) if len(factors) > 1 else 0,
            #float(len(factors)),
            #bin(q).count('0')
        ]

def factors_metrics_map(factors:dict[int,list[tuple[int,int]]], bitsize:int) -> dict[int,list[float]]:
    return {mx: factors_metrics(factors,bitsize) for mx, factors in factors.items()}

def factors_load(x:str) -> list[tuple[int,int]]:
    return [(int(prime), int(power)) for prime, power in json.loads(x)]

def factors_to_int(factors:list[tuple[int,int]]) -> int:
    product = 1
    for prime,power in factors:
        product *= (prime ** power)
    return product
