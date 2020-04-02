import msprime
import numpy as np
import pandas
import math
import gzip

class Simulator:
    def __init__(self, mapfile, demofile, sample_size, mu, rho = 1.2e-8, check = False):
        self.mu = mu
        self.sample_size = sample_size
        self.mapfile = mapfile
        self.demography = self.read_demo(demofile)
        self.rho = rho
        self.configure_demography()
        print("demographics:", demofile)
        print("sample size:", self.sample_size)
        if check: check_demographics()


    def read_demo(self, demofile):
        df = pandas.read_csv(demofile, sep=",", header=None)
        df.columns = ['generation', 'size']
        return df

    def configure_demography(self):
        self.demographic_events = []
        self.pc = [msprime.PopulationConfiguration(self.sample_size)]
        for index in self.demography.index:
            if index == self.demography.shape[0] - 1: break
            forward_time = self.demography['generation'][index + 1]
            forward_size = self.demography['size'][index + 1]
            now_time = self.demography['generation'][index]
            now_size = self.demography['size'][index]
            g = (math.log(now_size) - math.log(forward_size)) / (forward_time - 
                                                                 now_time)
            self.demographic_events.append(
            msprime.PopulationParametersChange(now_time, now_size, growth_rate=0))

    def check_demographics(self):
        dp = msprime.DemographyDebugger(
            population_configurations=self.pc, 
            demographic_events=self.demographic_events)
        dp.print_history()

    def simulation(self, len, seed, output=None):
        print("Simulating the trees ...")
        if (mapfile==None):
            print("constant recombination rate: ", self.rho)
            print("mutation rate: ", self.mu)
            tree_seq = msprime.simulate(
                population_configurations = self.pc,
                demographic_events = self.demographic_events,
                mutation_rate = self.mu,
                length=len,
                recombination_rate=self.rho,
                random_seed=seed,
                Ne=20000)
        else:
            print("Reading " + self.mapfile, "...")
            self.recomb_map = msprime.RecombinationMap.read_hapmap(mapfile)
            print("recombination map:", self.mapfile)
            tree_seq = msprime.simulate(
                population_configurations = self.pc,
                demographic_events = self.demographic_events,
                mutation_rate = self.mu,
#                length=len,
                recombination_map=self.recomb_map,
                random_seed=seed,
                Ne=20000)


            # tree_seq = msprime.simulate(
        #     sample_size=10, Ne=2e4, length=5e3, recombination_rate=1.2e-8,
        #     mutation_rate=1.65e-8, random_seed=10)
        if output != None:
            with gzip.open(output + ".vcf.gz", "wt") as vcf_file:
                tree_seq.write_vcf(vcf_file, ploidy=2)
            tree_seq.dump(output + ".tree")
        return tree_seq

    # def simulation(self, output):
    #     print("Simulating the trees ...")
    #     tree_seq = msprime.simulate(
    #         population_configurations = self.pc, 
    #         demographic_events = self.demographic_events,
    #         mutation_rate = self.mu, 
    #         recombination_map = self.recomb_map
    #     )
    #     # with open(output + ".vcf", "w") as vcf_file:
    #         # tree_seq.write_vcf(vcf_file, 2)
    #     # tree_seq.dump(output + ".tree")
    #     return tree_seq

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Simulate trees')
    parser.add_argument('-demo', required=True)
    parser.add_argument('-sample_size', default=2)
    parser.add_argument('-out', default="output")
    parser.add_argument('-n', default=64e6)
    parser.add_argument('-seed', default=10)
    parser.add_argument('-map', default=None)
    parser.add_argument('-rho', default=1.2e-8)
    args = vars(parser.parse_args())
    mapfile = args['map']
    demofile = args['demo']
    sim = Simulator(mapfile, demofile, sample_size=int(args['sample_size']), mu=1.65e-8, rho=1.2e-8)
    sim.simulation(len=int(args['n']), output=args['out'], seed=args['seed'])
