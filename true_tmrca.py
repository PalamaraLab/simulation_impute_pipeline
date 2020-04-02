import math
import os, re
import tskit
import msprime
import pandas as pd


def tmrca(trees, u, v, block_size):
    info = []
    for tree in trees.trees():
        info.append([u, v, tree.interval[1], tree.mrca(u, v), tree.time(tree.mrca(u, v))])
    info = pd.DataFrame(info)
    info.columns = ['id1','id2', 'interval', 'mrca', 'tmrca']
    truth = pd.DataFrame({'pos':pd.cut(info['interval'], [i for i in range(0, 10000001, block_size)], labels=False), 
                          'tmrca':info['tmrca']})
    truth = truth.groupby('pos').mean().reset_index()
    return(truth)


def simulate_test_tree():
    strees = msprime.simulate(sample_size=5, Ne=1000, length=4e4, 
                              mutation_rate=1e-8, random_seed=6, recombination_rate=2e-8)
    strees.dump("test_tmrca.tree")
    for tree in strees.trees():
        print("tree {}: interval = {}".format(tree.index, tree.interval))
        print(tree.draw(format="unicode"))
    return(strees)

def find_all_samples(dir):
    f = os.listdir(dir)
    f = list(filter(lambda x: re.search('map', x), f))
    f.sort()
    f = pd.Series(f)
    f = f.str.split("_")
    id1 = [x[1] for x in f]
    id2 = [x[3] for x in f]
    hapid1 = pd.Series([x[2] for x in f])
    hapid2 = pd.Series([x[4] for x in f])
    hapid1 = hapid1.str.split("-")
    hapid2 = hapid2.str.split(".")
    hapid1 = [x[0] for x in hapid1]
    hapid2 = [x[0] for x in hapid2]

    ff = pd.DataFrame({'id1': id1, 'id2':id2, 'hap1':hapid1, 'hap2':hapid2})
    ff['path'] = dir
    return(ff)


def rmse(data):
    max_est = data['tmrca_est'].max()
    data.loc[data['tmrca_truth'] > max_est, 'tmrca_truth'] = max_est
    e = (data['tmrca_est'] - data['tmrca_truth']) ** 2 / data.shape[0]
    return(math.sqrt(e.sum()))


def tmrca_batch(tree, sample):
    rmse_result = []
    for index, row in sample.iterrows():
        ind1 = int(row['id1']) * 2 + int(row['hap1'])
        ind2 = int(row['id2']) * 2 + int(row['hap2'])
        truth = tmrca(tree, ind1, ind2, 100)
        est = pd.read_csv(row['path'] + "/msp_" + row['id1'] + "_" + 
                          row['hap1'] + "-msp_" + row['id2'] + "_" + row['hap2'] + ".map", 
                          header = None)
        est['pos'] = [i for i in range(0, est.shape[0])]
        est.columns = ['tmrca', 'pos']
        data = est.merge(truth, left_on='pos', right_on='pos', suffixes=['_est', '_truth'])
        rmse_result.append([ind1, ind2, rmse(data)])
        
    rmse_result = pd.DataFrame(rmse_result)
    rmse_result.columns = ['id1','id2', 'rmse']
    rmse_result.to_csv("rmse.csv", index=None, header=True)


if __name__ == "__main__":
##    all_samples = find_all_samples("exp03_chip0.6_sim20k_ref1k_ref1000/imputed_with_ap")
#    tree = tskit.load("exp03_chip0.6_sim20k_ref1k_ref1000/sim.tree")
    all_samples = find_all_samples(".")
    tree = tskit.load("../sim.tree")
    tmrca_batch(tree, all_samples)
    

#msp_2925
#msp_19083
#    strees = simulate_test_tree()
#    result = tmrca(strees, 5850, 5851)
#    result = result.append(tmrca(strees, 38166, 38167))
#    result = result.append(tmrca(strees, 5850, 38166))
#    result = result.append(tmrca(strees, 5850, 38167))
#    result = result.append(tmrca(strees, 5851, 38166))
#    result = result.append(tmrca(strees, 5851, 38167))
#    result.to_csv("truth_tmrca.csv", index= None, header=True)


