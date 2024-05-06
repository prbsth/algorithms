import sys
import copy
import heapq
import math
import random

import numpy as np
# from matplotlib import pyplot as plt

MAX_ITER = 25000


def random_s(a):
    return [random.choice([-1, 1]) for i in range(len(a))]

def random_s_neighbor(s):
    return_s = copy.deepcopy(s)

    i = j = None
    while (i == j):
        i = random.randint(0, len(s)-1)
        j = random.randint(0, len(s)-1)

    return_s[i] = -return_s[i]
    if random.randint(0,1):
        return_s[j] = -return_s[j]
    
    return return_s

def residue(s, a):
    return abs(sum([s_i * a_i for s_i, a_i in zip(s, a)]))

#prepartioning
def random_p(n):
    return [random.randint(1, n) for _ in range(n)]

def random_p_neighbor(p):
    return_p = copy.deepcopy(p)
    i = random.randint(0, len(p)-1)
    j = random.randint(1, len(p))
    while p[i] == j:
        j = random.randint(1, len(p))
    return_p[i] = j
    return return_p

def prepartition_to_a_prime(p, a):
    a_prime = [0] * len(a)
    for j in range(len(p)):
        a_prime[p[j]-1] += a[j]
    return a_prime

def karmarkar_karp(a):
    h = [-i for i in a]  
    heapq.heapify(h)

    while len(h) > 1:
        first = -heapq.heappop(h)
        second = -heapq.heappop(h)
        if first != second:  #for efficiency, pushing 0 makes no difference
            heapq.heappush(h, -(first - second)) 

    return -h[0] if h else 0


def residue_prepartition(p, a):
    a_prime = prepartition_to_a_prime(p, a)
    return karmarkar_karp(a_prime)


def repeated_random(a):
    s = random_s(a)

    for i in range(MAX_ITER):
        s_prime = random_s(a)
        if residue(s_prime, a) < residue(s, a):
            s = s_prime

    return residue(s, a)

def hill_climbing(a):
    s = random_s(a)

    for i in range(MAX_ITER):
        s_prime = random_s_neighbor(s)
        if residue(s_prime, a) < residue(s, a):
            s = s_prime
    
    return residue(s, a)

def cooling_schedule(t):
    return (10.**10.*((0.8)**(t/300)))
    

def simulated_annealing(a):
    s = random_s(a)
    s_double_prime = copy.deepcopy(s)

    for i in range(MAX_ITER):
        s_prime = random_s_neighbor(s)
        if residue(s_prime, a) < residue(s, a):
            s = copy.deepcopy(s_prime)
        else:
            s = [s, s_prime][random.random() < math.exp(-(residue(s_prime, a)-residue(s, a))/cooling_schedule(i))]
        if residue(s, a) < residue(s_double_prime, a):
            s_double_prime = copy.deepcopy(s)
    
    return residue(s_double_prime, a)


def simulated_annealing_prepartition(a):
    p = random_p(len(a))
    p_double_prime = copy.deepcopy(p)
    
    for i in range(MAX_ITER):
        p_prime = random_p_neighbor(p)
        if residue_prepartition(p_prime, a) < residue_prepartition(p, a):
            p = copy.deepcopy(p_prime)
        else:
            p = [p, p_prime][random.random() < math.exp(-(residue_prepartition(p_prime, a)-residue_prepartition(p, a))/cooling_schedule(i))]
        if residue_prepartition(p, a) < residue_prepartition(p_double_prime, a):
            p_double_prime = copy.deepcopy(p)

    return karmarkar_karp(prepartition_to_a_prime(p_double_prime, a))


def repeated_random_prepartition(a):
    p = random_p(len(a))
    for i in range(MAX_ITER):
        p_prime = random_p(len(a))
        if residue_prepartition(p_prime, a) < residue_prepartition(p, a):
            p = p_prime

    return karmarkar_karp(prepartition_to_a_prime(p, a))

def hill_climbing_prepartition(a):
    p = random_p(len(a))
    for i in range(MAX_ITER):
        p_prime = random_p_neighbor(p)
        if residue_prepartition(p_prime, a) < residue_prepartition(p, a):
            p = p_prime

    return karmarkar_karp(prepartition_to_a_prime(p, a))

def generate_random_instance():
    return [random.randint(1, 10**12) for _ in range(100)]

def run_experiment(flag, algorithm, a):
    if algorithm == "0":
        return karmarkar_karp(a)
    elif algorithm == "1":
        return repeated_random(a)
    elif algorithm == "2":
        return hill_climbing(a)
    elif algorithm == "3":
        return simulated_annealing(a)
    elif algorithm == "11":
        return repeated_random_prepartition(a)
    elif algorithm == "12":
        return hill_climbing_prepartition(a)
    elif algorithm == "13":
        return simulated_annealing_prepartition(a)
    else:
        raise ValueError("Invalid algorithm code")

if __name__ == "__main__":
    flag = int(sys.argv[1])
    algorithm = sys.argv[2]
    input_file = sys.argv[3]

    if flag == 0:
        with open(input_file) as fin:
            a = [int(line.strip()) for line in fin.readlines()]

        result = run_experiment(flag, algorithm, a)
        print(result)

    else:
        # print(f"Residue for {algorithm}: {result}")
        # Generate 50 random instances and run experiments
        instances = [generate_random_instance() for _ in range(50)]
        
        kk_results = []
        rr_results = []
        hc_results = []
        sa_results = []
        prr_results = []
        phc_results = []
        psa_results = []
        
        for instance in instances:
            kk_results.append(karmarkar_karp(instance))
            rr_results.append(repeated_random(instance))
            hc_results.append(hill_climbing(instance))
            sa_results.append(simulated_annealing(instance))
            prr_results.append(repeated_random_prepartition(instance))
            phc_results.append(hill_climbing_prepartition(instance))
            psa_results.append(simulated_annealing_prepartition(instance))

        print("Algorithm\t\tMean Residue\t\tStd Dev")
        print(f"Karmarkar-Karp\t\t{np.mean(kk_results):.2f}\t\t{np.std(kk_results):.2f}")
        print(f"Repeated Random\t\t{np.mean(rr_results):.2f}\t\t{np.std(rr_results):.2f}")
        print(f"Hill Climbing\t\t{np.mean(hc_results):.2f}\t\t{np.std(hc_results):.2f}")
        print(f"Simulated Annealing\t{np.mean(sa_results):.2f}\t\t{np.std(sa_results):.2f}")
        print(f"Prep Repeated Random\t{np.mean(prr_results):.2f}\t\t{np.std(prr_results):.2f}")
        print(f"Prep Hill Climbing\t{np.mean(phc_results):.2f}\t\t{np.std(phc_results):.2f}")
        print(f"Prep Simulated Annealing\t{np.mean(psa_results):.2f}\t\t{np.std(psa_results):.2f}")
        
        # Plot graphs
        # plt.figure(figsize=(10, 6))
        # plt.boxplot([kk_results, rr_results, hc_results, sa_results, prr_results, phc_results, psa_results],
        #             labels=["KK", "RR", "HC", "SA", "PRR", "PHC", "PSA"])
        # plt.ylabel("Residue")
        # plt.title("Comparison of Algorithms")
        # plt.show()
