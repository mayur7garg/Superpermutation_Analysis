from typing import List, Tuple
from random import choices, randint, random, sample
from functools import reduce

def join_perms(left: str, right: str):
    for i in range(min(len(left), len(right)), 0, -1):
        if left[-i:] == right[:i]:
            return left + right[i:]
    return left + right

class StringSuperPerm_GA():
    def __init__(self, permutations : List[str], population_size: int, mutation_rate: float, elitism: int, 
    max_mutations: int = 1, weighted_selection: bool = False, first_n_unchanged: int = 0):

        self.permutations = permutations
        self.population_size = max(1, population_size)
        self.mutation_rate = max(0, mutation_rate)
        self.elitism = max(0, min(elitism, population_size))
        self.max_mutations = max(1, max_mutations)
        self.weighted_selection = weighted_selection
        self.first_n_unchanged = max(0, first_n_unchanged)

        self.current_generation: int = 0
        self.total_perms: int = len(self.permutations)

        base_DNA = [i for i in range(self.total_perms)]
        self.current_population = [base_DNA[: self.first_n_unchanged] + sample(base_DNA[self.first_n_unchanged :], self.total_perms - self.first_n_unchanged) for _ in range(self.population_size)]

        self.current_population_scores: List[Tuple[int, int]] = self.get_current_population_scores()
        self.best_dna_per_generation: List[List[int]] = [self.current_population[self.current_population_scores[0][1]]]
        self.best_scores_per_generation: List[int] = [self.current_population_scores[0][0]]

    def get_superpermutation(self, dna: List[int]):
        return reduce(join_perms, [self.permutations[gene] for gene in dna])
    
    def get_dna_score(self, id: int):
        return (len(self.get_superpermutation(self.current_population[id])), id)

    def get_current_population_scores(self):
        return sorted(list(map(self.get_dna_score, range(self.population_size))))

    def get_dna_weights(self):
        total_score = sum([score[0] for score in self.current_population_scores])
        return [total_score / score[0] for score in self.current_population_scores]

    def crossover(self, left_dna: List[int], right_dna: List[int]):
        child_dna = left_dna[: randint(0, self.total_perms)]

        for gene in right_dna:
            if gene not in child_dna:
                child_dna.append(gene)

        return child_dna

    def mutate(self, dna: List[int]):
        for i in range(randint(1, self.max_mutations)):

            if random() < self.mutation_rate:
                i = randint(self.first_n_unchanged, self.total_perms - 1)
                j = randint(self.first_n_unchanged, self.total_perms - 1)

                dna[i], dna[j] = dna[j], dna[i]

            else:
                break

        return dna

    def get_current_best_superpermutation(self):
        best_dna_score = min(self.best_scores_per_generation)
        best_dna_gen = self.best_scores_per_generation.index(best_dna_score)

        return (best_dna_gen, self.get_superpermutation(self.best_dna_per_generation[best_dna_gen]), best_dna_score)

    def get_new_population(self):
        new_population: List[List[int]] = []

        for i in range(self.elitism):
            new_population.append(self.current_population[self.current_population_scores[i][1]])

        if self.weighted_selection:
            for i in range(self.population_size - self.elitism):
                new_population.append(self.mutate(self.crossover(*choices(population = self.current_population, k = 2, weights = self.get_dna_weights()))))
        else:
            for i in range(self.population_size - self.elitism):
                new_population.append(self.mutate(self.crossover(*sample(self.current_population, 2))))

        return new_population

    def iterate_generation(self):
        self.current_population = self.get_new_population()
        self.current_population_scores = self.get_current_population_scores()
        self.best_dna_per_generation.append(self.current_population[self.current_population_scores[0][1]])
        self.best_scores_per_generation.append(self.current_population_scores[0][0])
        self.current_generation += 1

    def iterate_generations(self, generations):
        for _ in range(generations):
            self.iterate_generation()
