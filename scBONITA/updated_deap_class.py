from deap import tools, base, creator

from statistics import stdev, mean
import math
import copy
import matplotlib.pyplot as plt
import numpy as np
from itertools import chain, repeat, product
import operator
from random import sample, choice, randint, random
from collections import defaultdict
from scipy.sparse import csr_matrix
import logging
from alive_progress import alive_bar
import os

# -------- Genetic Algorithm Code --------
class CustomDeap:
    def __init__(self,network,network_name,dataset_name,binMat,nodeList,nodes,individual_length,nodeDict,successorNums):
        # Pass in the node objects
        self.nodes = nodes

        # Genetic algorithm parameters
        self.mutate_percent_pop = 0.25
        self.generations = 15
        self.starting_population_size = 50
        self.parent_population_size = 25
        self.child_population_size = 25
        self.crossover_probability = 0.1
        self.mutation_probability = 0.9
        self.bitFlipProb = 0.5
        self.nodeDict = nodeDict
        self.successorNums = successorNums

        # General parameters
        self.network = network
        self.network_name = network_name
        self.dataset_name = dataset_name
        self.binMat = binMat
        self.nodeList = nodeList
        self.stats = tools.Statistics(key=self.get_fitness_values)
        self.individual_length = individual_length
        self.individualParse = [0] + [node.rule_end_index if node.rule_end_index else nodes[i-1].rule_end_index for i, node in enumerate(nodes)]
        self.size = len(self.individualParse)
        self.node_combination_length = [len(node.possibilities) for node in nodes]
        self.total_combinations = [node.possibilities for node in nodes]
        self.total_inversions = [node.inversions for node in nodes]
        self.make_toolbox()

        logging.basicConfig(format='%(message)s', level=logging.INFO)

        _, num_columns = np.shape(self.binMat)

        # Chunk the data matrix to reduce noise and put into numpy array to speed up processing
        self.num_chunks = 100

        # Chunk if there are more cells than columns, otherwise just use the columns
        if num_columns > self.num_chunks:
            self.chunked_data_numpy = np.array(self.chunk_data(num_chunks=self.num_chunks))
            self.coarse_chunked_dataset = np.array(self.chunk_data(num_chunks=round(self.num_chunks / 2, 1)))
        else:
            self.num_chunks = num_columns
            self.chunked_data_numpy = np.array(self.chunk_data(num_chunks=num_columns))
            self.coarse_chunked_dataset = np.array(self.chunk_data(num_chunks=num_columns))

        logging.debug(f'self.chunked_data_numpy: {self.chunked_data_numpy}')

    # 1. Main Genetic Algorithm Function
    def genetic_algorithm(self):
        logbook = tools.Logbook()

        logging.info(f'\n-----CHUNKING DATASET-----')
        num_rows, num_columns = np.shape(self.binMat)
        logging.info(f"\tOriginal Data Shape: {num_rows} rows, {num_columns} columns")

        chunked_rows, chunked_columns = np.shape(self.chunked_data_numpy)
        logging.info(f"\tChunked Data Shape: {chunked_rows} rows, {chunked_columns} columns")

        coarse_chunked_rows, coarse_chunked_columns = np.shape(self.coarse_chunked_dataset)
        logging.info(f"\tCoarse Chunked Data Shape: {coarse_chunked_rows} rows, {coarse_chunked_columns} columns")

        logging.info(f'\n-----GENETIC ALGORITHM-----')
        population = self.toolbox.population(n=self.starting_population_size)

        total_fitnesses = []

        logbook.header = ["gen", "nevals"] + (self.stats.fields if self.stats else [])
        lastcheck = []
        modellist = []
        fitnesslist = []
        popList = []

        # Evaluate the individuals with an invalid fitness
        invalid_ind = [ind for ind in population if not ind.fitness.valid]

        raw_fitnesses, fitnesses = self.fitness_calculation(invalid_ind, self.chunked_data_numpy, current_gen=0, max_gen=self.generations)

        for ind, fit in zip(invalid_ind, fitnesses):
            ind.fitness.values = fit

        # Append the fitnesses for graphing
        total_fitnesses.append(raw_fitnesses)

        fitnesslist.append([list(ind.fitness.values) for ind in population])
        popList.append([list(inder[1]) for inder in population])
        modellist.append(
            [
                [
                    (modeler[0].size),
                    list(modeler[0].nodeList),
                    list(modeler[0].individualParse),
                    list(modeler[0].total_combinations),
                    list(modeler[0].total_inversions),
                    list(modeler[0].node_combination_length),
                    list(modeler[0].nodeList),
                    dict(modeler[0].nodeDict),
                ]
                for modeler in population
            ]
        )

        # Calculate values to display to the terminal for the initial generation
        average_fitness = round(sum(raw_fitnesses) / len(raw_fitnesses), 3)
        min_fitness = round(min(raw_fitnesses), 3)
        max_fitness = round(max(raw_fitnesses), 3)
        stdev_fitness = round(stdev(raw_fitnesses), 3)

        logging.info(f'ngen\tnevals\tavg\tstd\tmin\tmax')
        logging.info(f'{0}\t{len(raw_fitnesses)}\t{average_fitness}\t{stdev_fitness}\t{min_fitness}\t{max_fitness}')

        # Begin the generational process
        for gen in range(1, self.generations + 1):

            # Perform mutation and crossover
            offspring = self.__varOrAdaptive(
                population,
                self.toolbox,
                (1.0 * gen / self.generations),
                self.mutate_percent_pop
            )

            invalid_offspring = [ind for ind in offspring if not ind.fitness.valid]

            # Calculate the fitness for each individual
            raw_fitnesses, fitnesses = self.fitness_calculation(invalid_offspring, self.chunked_data_numpy, gen, self.generations+1)

            for ind, fit in zip(invalid_offspring, fitnesses):
                ind.fitness.values = fit

            # Append the fitnesses for graphing
            total_fitnesses.append(raw_fitnesses)

            if gen == self.generations:
                # self.graph_results(total_fitnesses)
                with open("deap_results.txt", "w") as temp_file:

                    joined_generations = ",".join([str(i) for i in range(self.generations+1)])
                    joined_fitness = ",".join([str(i) for i in total_fitnesses])
                    temp_file.write(str(joined_generations))
                    temp_file.write(str(joined_fitness))
                return raw_fitnesses, invalid_offspring, logbook
            else:

                # Select the next generation population
                population[:] = self.toolbox.select(offspring, self.parent_population_size)

                # Get the fitness
                fitnesslist.append([list(ind.fitness.values) for ind in population])
                popList.append([list(inder[1]) for inder in population])
                modellist.append(
                    [
                        [
                            (modeler[0].size),
                            list(modeler[0].nodeList),
                            list(modeler[0].individualParse),
                            list(modeler[0].total_combinations),
                            list(modeler[0].total_inversions),
                            list(modeler[0].node_combination_length),
                            list(modeler[0].nodeList),
                            dict(modeler[0].nodeDict),
                        ]
                        for modeler in population
                    ]
                )

                # Calculate values to display to the terminal for the initial generation
                average_fitness = round(sum(raw_fitnesses) / len(raw_fitnesses), 3)
                min_fitness = round(min(raw_fitnesses), 3)
                max_fitness = round(max(raw_fitnesses), 3)
                stdev_fitness = round(stdev(raw_fitnesses), 3)

                logging.info(f'{gen}\t{len(population)+len(raw_fitnesses)}\t{average_fitness}\t{stdev_fitness}\t{min_fitness}\t{max_fitness}')


    def graph_genetic_algorithm_results(self, total_fitnesses):
        x = [i for i in range(self.generations + 1)]
        y = total_fitnesses
        plt.xlabel("Generations")
        plt.ylabel("Error")
        plt.title("Error over time")
        for i in range(len(y[0])):
            plt.plot(x, [pt[i] for pt in y])
        plt.ylim([0, 0.5])
        plt.legend()
        plt.show()


    def find_best_individual(self, population, raw_fitnesses):
        # Set the fitness values to it's unweighted fitness
        for i, ind in enumerate(population):
            ind.fitness.values = (raw_fitnesses[i],)

        # Refine the rules
        best_rulesets, ruleset_errors = self.refine_rules(population, raw_fitnesses)

        logging.debug(f'ruleset errors: {ruleset_errors}')

        # If there is more than one ruleset, sort the list by the minimum error
        if len(best_rulesets) > 1:
            sorted_rulesets, sorted_ruleset_errors = zip(*sorted(zip(best_rulesets, ruleset_errors)))     

            logging.debug(f'sorted ruleset errors: {sorted_ruleset_errors}')
        else:
            sorted_rulesets, sorted_ruleset_errors = best_rulesets, ruleset_errors
        
        # Write out the rulesets to the output file and terminal
        individual_num = 1
        for ruleset in sorted_rulesets:
            logging.info(f'\nEquivalent Ruleset {individual_num}')
            self.write_ruleset(ruleset, sorted_ruleset_errors[individual_num-1], individual_num)
            individual_num += 1
        
        # Set the best rules for the nodes based on which ruleset has the lowest error
        best_ruleset = sorted_rulesets[0]
        for node_index, node in enumerate(self.nodes):
            node_best_rule = best_ruleset[node_index][0]
            node_best_rule_index = best_ruleset[node_index][1]

            node.best_rule = node_best_rule
            node.best_rule_index = node_best_rule_index
            logging.debug(f'Node {node.name}')
            logging.debug(f'\tnode best rule {node.best_rule}')
            logging.debug(f'\tnode best rule index {node.best_rule_index}')
        
        return best_ruleset

    def write_ruleset(self, ruleset, error, individual_num):
        # Write out the rules to an output file
        rule_path = f'rules_output/{self.dataset_name}_rules/{self.network_name}_{self.dataset_name}_ind_{individual_num}_rules.txt'
        os.makedirs(f'rules_output/{self.dataset_name}_rules/', exist_ok=True)
        with open(rule_path, 'w') as rule_file:
            for rule, _, _ in ruleset:
                rule_name = rule[0]
                incoming_nodes = rule[1]
                logic = rule[2]
                logic = logic.split("_")
                if len(logic) == 5:
                    NOT_1 = self.invert_rules(rule, 0)
                    NOT_2 = self.invert_rules(rule, 1)
                    NOT_3 = self.invert_rules(rule, 2)
                    line = f'{rule_name} = {NOT_1}{incoming_nodes[0]} {logic[1]} {NOT_2}{incoming_nodes[1]} {logic[3]} {NOT_3}{incoming_nodes[2]}'
                    logging.info(f'\t{line}')
                    rule_file.write(line)
                    rule_file.write('\n')
                elif len(logic) == 3:
                    NOT_1 = self.invert_rules(rule, 0)
                    NOT_2 = self.invert_rules(rule, 1)
                    line = f'{rule_name} = {NOT_1}{incoming_nodes[0]} {logic[1]} {NOT_2}{incoming_nodes[1]}'
                    logging.info(f'\t{line}')
                    rule_file.write(line)
                    rule_file.write('\n')
                else:
                    NOT_1 = self.invert_rules(rule, 0)
                    line = f'{rule_name} = {NOT_1}{incoming_nodes[0]}'
                    logging.info(f'\t{line}')
                    rule_file.write(line)
                    rule_file.write('\n')
            logging.info(f'Refined Error: {error}\n')
            rule_file.write(f'Refined_error:\t{error}')

    def invert_rules(self, rule, index):
        if rule[3][index] == True:
            return 'NOT '
        else:
            return ''

    def get_fitness_values(self, ind):
        """
        Gets around needing to use a lambda function in tools.Statistics so that the ruleset object can be pickled
        """
        return ind.fitness.values

    def refine_rules(self, population, raw_fitnesses):

        logging.info(f'\n-----REFINING RULESETS-----')
        min_error_individuals = [ind for ind in population if ind.fitness.values == min(raw_fitnesses)]

        # Set to keep track of unique value combinations
        seen_values = set()

        # List comprehension to filter out duplicates
        unique_min_error_individuals = []
        for ind in min_error_individuals:
            # Convert the individual's values to a tuple (or another immutable type if necessary)
            values_tuple = tuple(ind[1])

            # Check if we've already seen these values
            if values_tuple not in seen_values:
                unique_min_error_individuals.append(ind)
                seen_values.add(values_tuple)

        # Now, unique_min_error_individuals contains only individuals with unique values
        all_best_rules = []
        ruleset_error = []


        # Iterate through each ruleset with the minimum error
        for i, individual_ruleset in enumerate(unique_min_error_individuals):
            logging.info(f'Equivalent Ruleset {i+1} / {len(unique_min_error_individuals)}')
            total_rules = []
            rules_per_node = []
            best_rules = []
            max_iterations = 3

            # Loop through each node
            with alive_bar(len(self.nodes), title='Refining Gene Rules') as bar:
                for index, node in enumerate(self.nodes):
                    continue_to_next_node = False
                    iteration = 1
                    passed_error_cutoff = True
                    individual = individual_ruleset[1]

                    # Check to see if we should move to the next node
                    while continue_to_next_node == False:
                        node.reset_state() # Reset the state of the node predictions for this ruleset
                        if iteration == 1:
                            logging.debug(f'\nRefining rules for node: {node.name} ({index+1} / {len(self.nodes)})')

                        # If this is the first iteration through the refinement process for this node, use the ruleset as is
                        if iteration != 1:
                            logging.debug(f'\tIteration {iteration}')
                            individual = [1 for _ in individual_ruleset[1]]

                        # Get the rules for each of the lowest error individuals
                        if node.rule_start_index is not None: # Make sure this is not an empty node
                            rule_length = node.rule_end_index - node.rule_start_index

                            # Make the predictions based on which rules are predicted by the ruleset
                            if rule_length >= 1:
                                node.rule_predictions = individual[node.rule_start_index:node.rule_end_index]
                            elif rule_length == 0:
                                node.rule_predictions = individual[node.rule_start_index]

                            # Find which rules for the node are predicted by the individual
                            node.find_rule_predictions()
                            
                            # Make and format the rules for the node
                            node.make_AND_OR_NOT_rules()

                            # Append the rules for the node to a list of the total rules and a list of the number of rules predicted for this node
                            total_rules.append(node.node_rules)
                            rules_per_node.append(len(node.node_rules))

                            # rules = (None, None)
                            if passed_error_cutoff == True:
                                logging.debug(f'\tPassed Error cutoff')
                                rules = self.calculate_refined_errors(node, self.chunked_data_numpy)
                                logging.debug(f'\t\tRules: {rules}')
                                
                            # If the node did not pass the error cutoff, look at every rule combination coarsely,
                            # increasing for each iteration
                            else:
                                rules = self.calculate_refined_errors(node, self.coarse_chunked_dataset)
                                for rule, index, error in rules:
                                    logging.debug(f'\tCOARSE SEARCH: rule {rule}, error = {error}')
                                # If every node is wrong, try with different inversion rules
                                for rule, index, error in rules:
                                    if error == 1:
                                        logging.debug(f'\tError is 1, flipping inversion rule for incoming nodes')
                                        for key, value in node.inversions.items():
                                            logging.debug(f'\tInversion rule before: {node.inversions[key]}')
                                            if value == True:
                                                node.inversions[key] = False
                                            elif value == False:
                                                node.inversions[key] = True
                                            logging.debug(f'\tInversion rule after: {node.inversions[key]}')

                                        # Re-make the rules for the node
                                        node.make_AND_OR_NOT_rules()

                                        rules = self.calculate_refined_errors(node, self.coarse_chunked_dataset)
                                        for rule, index, error in rules:
                                            logging.debug(f'\tAFTER INVERSION: rule {rule}, error = {error}')
                                
                                # random_equivalent_rule = rules[random.choice(len(rules))]
                            # Check if the error is below the cutoff
                            for rule, index, error in rules:
                                # print(f'Checking rule {rule}, error = {error}')
                                if error < 0.10:
                                    logging.debug(f'\tError is under 0.10: {error}')
                                    logging.debug(f'\trule: {rule}, index: {index}, error: {error}')
                                    best_rules.append((rule, index, error))
                                    passed_error_cutoff = True
                                    continue_to_next_node = True
                                    break

                                else:
                                    passed_error_cutoff = False
                                    continue_to_next_node = False

                            # If the error is too high, iterate again
                            if not continue_to_next_node:
                                logging.debug(f'\tDid not pass error cutoff, increasing iteration')
                                iteration += 1

                            # Break if no best rule is found (changed to leave in high error node rules)
                            if iteration >= max_iterations:
                                logging.debug(f'\tMax iterations reached, finding the minimum error rule')

                                # Gets around using lambda with sorted so that the ruleset object can be pickled
                                def get_third_element(item):
                                    return item[2]

                                min_error_rule = sorted(rules, key=get_third_element)
                                logging.debug(f'\t\tMin error rule: {min_error_rule}')
                                best_rules.append(min_error_rule[0])
                                continue_to_next_node = True
                        
                        # If the node has no predictions, set the nodes expression to itself
                        else:
                            logging.debug(f'\tException: {node.name} has no incoming nodes, signaling rule set to itself')
                            if len(node.node_rules) == 0:
                                # Add in self signaling for the node
                                logging.debug(f'\t[{node.name}, [{node.name}], [False]]')
                                node.node_rules = [node.name, [node.name], 'A', [False]]
                                node.predecessors[node.index] = node.name
                                node.inversions[node.index] = False

                                # Set the node's incoming node rules to itself
                                rule = [node.name, [node.name], 'A', [False]]
                                index = 0
                                error = 0
                                
                                best_rules.append((rule, index, error))
                                continue_to_next_node = True
                                # print(f'\tNo predictions\n')

                    bar()

            all_best_rules.append(best_rules)

            # print('FINAL BEST RULES:')
            final_error = 0
            for rule, index, error in best_rules:
                logging.debug(f'{rule}: Error = {error}')
                final_error += error
            ruleset_error.append(final_error / len(best_rules))


        return (all_best_rules, ruleset_error)

    def calculate_error(self, node, predicted_rule, dataset):
        # Get the row in the dataset for the node being evaluated
        node_evaluated = dataset[node.index]

        # Get the dataset row indices for the incoming nodes included in this rule
        incoming_node_indices = [index for index, name in node.predecessors.items() if name in predicted_rule[1]]

        # Check if an index exists
        def index_exists(lst, index):
            try:
                _ = lst[index]
                return True
            except IndexError:
                return False
        
        # Set A, B, and C to 0
        A = 0
        B = 0
        C = 0
        
        # Update the values of A, B, and C based on the gene index in the dataset
        if index_exists(incoming_node_indices, 0):
            A = dataset[incoming_node_indices[0]]
        if index_exists(incoming_node_indices, 1):
            B = dataset[incoming_node_indices[1]]
        if index_exists(incoming_node_indices, 2):
            C = dataset[incoming_node_indices[2]]

        # Allows for arrays to be passed into the function rather than one value at a time
        predicted = np.vectorize(eval(predicted_rule))

        # Sum the number of differences between the predicted output and the actual data
        difference = np.sum(predicted != node_evaluated)

        # Store the number of total predictions for error calculation (different / total)
        count = len(predicted)

        return difference, count

    def calculate_refined_errors(self, node, chunked_dataset):
        # Calculate the error for each prediction based on the entire dataset (not chunked)
        prediction_errors = []
        most_incoming_node_rules = []
        best_rule_errors = []
        best_rule_indices = []
        best_rules = []

        logging.debug(f'\tPredicted rule errors:')
        for predicted_rule in node.node_rules:
            logging.debug(f'\t\t{predicted_rule}')
            logging.debug(f'\t\tNode rules: {node.node_rules}')

            difference, count = self.calculate_error(node, predicted_rule, chunked_dataset)
            prediction_error = difference / count
            prediction_errors.append(prediction_error)
            logging.debug(f'\t\t\tError: {prediction_error}')

        # Finding indices of all occurrences of the minimum value
        logging.debug(f'\t\tPrediction_errors: {prediction_errors}')
        if prediction_errors: # Excludes nodes with no predictions
            logging.debug(f'\t\tPrediction errors: {prediction_errors}')
            minimum_error = min(prediction_errors)
            minimum_error_indices = [index for index, value in enumerate(prediction_errors) if
                                     value == minimum_error]

            # Calculate which rules have the lowest index, prioritize nodes with more incoming nodes
            # (If three incoming nodes have the same error as two incoming nodes, the two incoming node rule is
            # probably a subset of the three incoming node rule)

            logging.debug(f'\n\tMinimum error indices: {minimum_error_indices}, minimum error: {minimum_error}')
            # Append the length of the rules and the index to a list of node lengths
            num_incoming_nodes = [len(node.node_rules[value][1]) for value in minimum_error_indices]

            # Find the maximum number of incoming nodes
            max_incoming_nodes = max(num_incoming_nodes)

            for index in minimum_error_indices:  # iterate through each node
                num_incoming_nodes = len(node.node_rules[index][1])
                logging.debug(f'\t\tminimum error index: {index}')
                logging.debug(f'\t\t\tmax_incoming_nodes = {max_incoming_nodes}')
                logging.debug(f'\t\t\tnum_incoming_nodes = {num_incoming_nodes}')
                # One incoming node rules in the best ruleset
                if max_incoming_nodes == 1:  # If there is one incoming node possibilities in the lowest error rules
                    if num_incoming_nodes == 1:  # Only print the rules with one incoming nodes
                        logging.debug(f'\t\tbest_rule: {node.node_rules[index]}, Index: {index}, error: {prediction_errors[index]}')
                        most_incoming_node_rules.append(node.node_rules[index])
                        best_rule_indices.append(index)
                        best_rule_errors.append(prediction_errors[index])

                # If there are no minimum error rules with one incoming node, look for rules with two incoming nodes
                elif max_incoming_nodes == 2:
                    if num_incoming_nodes == 2:  # Only print the rules with two incoming nodes
                        logging.debug(f'\t\tbest_rule: {node.node_rules[index]}, Index: {index}, error: {prediction_errors[index]}')
                        most_incoming_node_rules.append(node.node_rules[index])
                        best_rule_indices.append(index)
                        best_rule_errors.append(prediction_errors[index])

                # If no three or two incoming node rules, print the rule and the number of incoming nodes
                elif max_incoming_nodes == 3:
                    if num_incoming_nodes == 3:
                        logging.debug(f'\t\tbest_rule: {node.node_rules[index]}, Index: {index}, error: {prediction_errors[index]}')
                        most_incoming_node_rules.append(node.node_rules[index])
                        best_rule_indices.append(index)
                        best_rule_errors.append(prediction_errors[index])
            
            logging.debug(f'\t\tMost incoming node rules: {most_incoming_node_rules}')
            # Find the best rules
            min_rules = min([rule[2].count("AND") for rule in most_incoming_node_rules])

            logging.debug(f'\tmost_incoming_node_rules: {most_incoming_node_rules}')
            for i, rule in enumerate(most_incoming_node_rules):
                if rule[2].count("AND") == min_rules:
                    best_rules.append((most_incoming_node_rules[i], best_rule_indices[i], best_rule_errors[i]))
            
            logging.debug('\tbest rules:')
            for i in best_rules:
                logging.debug(f'\t\t{i[0]}')
                logging.debug(f'\t\tError = {i[1]}')


        return best_rules

    def chunk_data(self, num_chunks):
        """
        Chunks the data by breaking the data into num_chunks number of chunks and taking the average value of all
        columns (cells) within the chunk for each row (genes). For each row, if the average value of all cells in
        the chunk is > 0.5, the chunk is set to 1 for that row. If the average is < 0.5, the chunk is set to 0.
        :param binMat:
        :param num_chunks:
        :return chunked_data:
        """

        num_chunks = int(num_chunks)

        # Shuffle the columns to randomize the order of cells in the chunks
        np.random.seed(42)  # Optional: for reproducibility. Will shuffle in the same way every time
        # chunked_data = self.binMat.todense()

        column_permutation = np.random.permutation(self.binMat.shape[1])
        shuffled_binMat = self.binMat[:, column_permutation]

        # Get the shape of the data
        num_rows, num_columns = np.shape(shuffled_binMat)

        # Define the chunk size by the number of cells and number of chunks
        # Ensure num_chunks is not greater than num_columns
        if num_chunks > num_columns:
            num_chunks = num_columns

        # Calculate chunk size
        chunk_size = max(1, num_columns // num_chunks)  # Ensure chunk_size is at least 1

        # Initiate a blank matrix filled with 0's of size: number of genes (rows) x number of chunks (columns)
        chunked_data = np.zeros((int(num_rows), int(num_chunks)))

        # Iterate through each of the chunks and take the average value of all columns within the chunks
        for chunk_index in range(int(num_chunks)):
            # Define the start and end column for this chunk
            start_column = chunk_index * chunk_size  # Current chunk * chunk size = start
            end_column = start_column + chunk_size  # End column defined by the start + chunk size

            # Find the columns between the start and end columns
            subset = shuffled_binMat[:, start_column:end_column]

            # Get the average row value for all columns in the chunk
            row_avg = np.mean(subset, axis=1)

            # Binarize the data
            row_avg[row_avg >= 0.5] = 1
            row_avg[row_avg < 0.5] = 0

            # Flattens the chunked array so that it is one column, add it to the chunked_data at the right index
            chunked_data[:, chunk_index] = row_avg.flatten()

        return chunked_data

    # 1.1 Sets up the toolbox from Deap to monitor the progress of the genetic algorithm
    def make_toolbox(self):
        """sets up GA toolbox from deap"""

        toolbox = base.Toolbox()
        weightTup = (-1.0,)  # specify weights of the errors
        for i in range(len(self.nodes) - 1):
            weightTup += (-1.0,)

        # MAKE TYPES
        creator.create(
            "FitnessMin", base.Fitness, weights=[-1.0]  # MODIFIED
        )  # make a fitness minimization function #the objective function has to be MINIMIZED

        creator.create(
            "individual", list, fitness=creator.FitnessMin
        )  # create a class of individuals that are lists of floats

        # INITIALIZATION
        # register our bitstring generator and how to create an individual, population
        toolbox.register("genRandomBitString", self.__genBits)  # , model=self)
        toolbox.register(
            "individual",
            tools.initIterate,
            creator.individual,
            toolbox.genRandomBitString,
        )
        toolbox.register("population", tools.initRepeat, list, toolbox.individual)

        # REGISTER STATISTICS
        # create statistics toolbox and give it functions
        stats = tools.Statistics(key=self.get_fitness_values)
        stats.register("avg", np.mean)
        stats.register("std", np.std)
        stats.register("min", np.min)
        stats.register("max", np.max)

        # REGISTER CROSSOVER, MUTATION, AND SELECTION FUNCTIONS
        # finish registering the toolbox functions
        toolbox.register("mate", tools.cxTwoPoint)
        toolbox.register("mutate", tools.mutFlipBit, indpb=self.bitFlipProb)
        toolbox.register("select", self.__selNSGA2)
        toolbox.register("similar", np.array_equal)

        self.toolbox = toolbox
        self.stats = stats

    # 1.1.1 Generate a random bitlist for the individual
    def __genBits(self):
        """Generates a random bitlist for the individual, makes sure the bitlist is not 0 and there are
        between 1 and 5 ones in the random list
        :return: [copy.deepcopy(self), startInd]
        """
        startInd = list(self.__genRandBits())
        counter = 0

        # make sure bitlist isn't zero
        while np.sum(startInd) == 0 and counter < float("Inf"):
            startInd = list(self.__genRandBits())
            counter += 1

        # go through nodes and make sure that there are 1-5 ones in the random list
        for index, node in enumerate(self.nodes):
            if node.rule_start_index:
                end = node.rule_end_index
                start = node.rule_start_index
                if (end - start) > 1:
                    counter = 0

                    while np.sum(startInd[start:end]) > 5 and counter < float("Inf"):
                        chosen = math.floor(random() * (end - start))
                        startInd[start + int(chosen)] = 0
                        counter += 1
                    if np.sum(startInd[start:end]) == 0:
                        chosen = math.floor(random() * (end - start))
                        startInd[start + int(chosen)] = 1
                elif (end - start) == 1:
                    startInd[start] = 1

        return [copy.deepcopy(self), startInd]

    # 1.1.1.1 Generates a random string of 1's and 0's for the individuals with length = number of possible combinations
    def __genRandBits(self):
        """Generates a random bitstring for the individual with a length equal to totalLenList (all possible
        incoming node combinations for each node in the network)"""
        bitstring = np.random.randint(2, size=(int(self.individual_length),))
        return list(bitstring)

    # 1.1.2 Calculates the fitness for an individual
    def __selNSGA2(self, individuals, k):
        """
        Calculate  fitness for an individual.

        NSGA2 selection taken from deap Apply NSGA-II selection operator on
        the *individuals*. Usually, the size of *individuals* will be larger
        than *k* because any individual present in *individuals* will appear
        in the returned list at most once.

        Having the size of *individuals* equals to *k* will have no effect other
        than sorting the population according to their front rank. The
        list returned contains references to the input *individuals*.

        :param individuals: A list of individuals to select from.
        :param k: The number of individuals to select.
        :returns: A list of selected individuals.
        """
        pareto_fronts = self.__sortNondominatedAdapt(individuals, k)
        for front in pareto_fronts:
            self.__assignCrowdingDist(front)

        chosen = list(chain(*pareto_fronts[:-1]))
        k = k - len(chosen)
        if k > 0:
            sorted_front = sorted(
                pareto_fronts[-1], key=operator.attrgetter("fitness.crowding_dist"), reverse=True
            )
            chosen.extend(sorted_front[:k])

        return chosen

    # 1.1.2.1 Sort individuals into Pareto fronts based on their fitness values
    def __sortNondominatedAdapt(self, individuals, k, first_front_only=False):
        """
        Taken from deap and modified slightly to make pareto sorting less strict.

         Sort the first *k* *individuals* into different nondomination levels
        using the "Fast Nondominated Sorting Approach"

        Time complexity of :math:`O(MN^2)`
            - :math:`M` is the number of objectives
            - :math:`N` the number of individuals.

        A **Pareto front** is a set of solutions that are not dominated
        by any other solutions in the population. The first Pareto front is the
        set of best solutions, and subsequent fronts are progressively less
        optional.

        :param individuals: A list of individuals to select from.
        :param k: The number of individuals to select.
        :param first_front_only: If :obj:`True` sort only the first front and
                                exit.
        :returns: A list of Pareto fronts (lists), the first list includes
                nondominated individuals.
        """
        if k == 0:
            return []

        # 1. Prepares mapping from fitness values to individuals sharing that fitness
        map_fit_ind = defaultdict(list)
        for ind in individuals:
            map_fit_ind[ind.fitness].append(ind)

        # 2. Initialize data structures to keep track of dominaitng and dominated fitness values
        fits = list(map_fit_ind)
        current_front = []
        next_front = []
        dominating_fits = defaultdict(int)
        dominated_fits = defaultdict(list)

        # 3. Rank first Pareto front
        # First Pareto front identified by finding individuals not dominated by any others. Added to current_front
        for i, fit_i in enumerate(fits):
            for fit_j in fits[i + 1:]:
                if self.__dominated(fit_i, fit_j):
                    dominating_fits[fit_j] += 1
                    dominated_fits[fit_i].append(fit_j)
                elif self.__dominated(fit_j, fit_i):
                    dominating_fits[fit_i] += 1
                    dominated_fits[fit_j].append(fit_i)
            if dominating_fits[fit_i] == 0:
                current_front.append(fit_i)

        fronts = [[]]
        for fit in current_front:
            fronts[-1].extend(map_fit_ind[fit])
        pareto_sorted = len(fronts[-1])

        # 4. Rank the next front until all individuals are sorted or
        # the given number of individual are sorted.
        if not first_front_only:
            N = min(len(individuals), k)
            while pareto_sorted < N:
                # Find the next Pareto front
                fronts.append([])
                for fit_p in current_front:
                    for fit_d in dominated_fits[fit_p]:
                        dominating_fits[fit_d] -= 1
                        if dominating_fits[fit_d] == 0:
                            next_front.append(fit_d)
                            pareto_sorted += len(map_fit_ind[fit_d])
                            fronts[-1].extend(map_fit_ind[fit_d])
                current_front = next_front
                next_front = []
        return fronts

    # 1.1.2.1.1 Looks to see if the weight values for ind1 are greater than the weight values of ind2
    def __dominated(self, ind1, ind2):
        """
        Taken from deap and modified slightly to make pareto sorting less strict.
        Return true if each objective of *self* is not strictly worse than the
        corresponding objective of *other* and at least one objective is strictly better.

        Slice indicating on which objectives the domination is tested. The default value
        is `slice(None)`, representing every objectives.

         :return: not_equal (True if weight values for ind1 > ind2 else False)
        """

        not_equal = False
        mean1 = np.mean(ind1.wvalues)
        mean2 = np.mean(ind2.wvalues)

        if mean1 > mean2:
            not_equal = True
        elif mean1 < mean2:
            return False
        return not_equal

    # 1.1.2.2 Assigns a crowding distance to each individual's fitness
    def __assignCrowdingDist(self, individuals):
        """taken from deap. Assign a crowding distance to each individual's fitness. The
        crowding distance can be retrieve via the :attr:`crowding_dist`
        attribute of each individual's fitness.
        """
        if len(individuals) == 0:
            return

        distances = [0.0] * len(individuals)
        crowd = [(ind.fitness.values, i) for i, ind in enumerate(individuals)]

        nobj = len(individuals[0].fitness.values)

        for i in range(nobj):

            # Gets around having to use a lambda funciton in crowd.sort so that the ruleset object can be pickled
            def get_key_function(i):
                def key_function(element):
                    return element[0][i]
                return key_function


            crowd.sort(key=get_key_function(i))
            distances[crowd[0][1]] = float("inf")
            distances[crowd[-1][1]] = float("inf")
            if crowd[-1][0][i] == crowd[0][0][i]:
                continue
            norm = nobj * float(crowd[-1][0][i] - crowd[0][0][i])
            for prev, cur, next in zip(crowd[:-2], crowd[1:-1], crowd[2:]):
                distances[cur[1]] += 1.0 * (next[0][i] - prev[0][i]) / norm

        for i, dist in enumerate(distances):
            individuals[i].fitness.crowding_dist = dist

    # 1.2 Calculates the fitness of the individual based on how closely the predicted rules match the data
    def fitness_calculation(self, invalid_ind, dataset, current_gen, max_gen):
        # ----- Variables -----

        fitnesses = []
        raw_fitnesses = []

        for individualRuleset in invalid_ind:
            # --------------------- Create the rulesets from the predictions ------------------------

            total_rules = []
            rules_per_node = []

            total_difference = 0
            total_count = 0

            # Calculate the error for each node in the individual
            for node in self.nodes:

                # Reset the state of the node
                node.reset_state()

                individual_bitstring = individualRuleset[1]

                node.find_multiple_rule_predictions(individual_bitstring)
                
                # If the node is in the individual
                if node.rule_start_index is not None:

                    rule_predictions = node.find_multiple_rule_predictions(individual_bitstring)

                    # For each predicted rule
                    for predicted_rule in rule_predictions:
                        # Calculate the error of that rule based on the dataset
                        difference, count = self.calculate_error(node, predicted_rule, dataset)
                    
                        total_difference += difference
                        total_count += count
                    
                    total_rules.append(node.node_rules)
                    rules_per_node.append(len(node.node_rules))

            # Calculate the fitness score based on the total error of all predicted rules
            fitness_score = total_difference / total_count

            raw_fitnesses.append(fitness_score)

            # Add a complexity penalty that increases over the generations
            complexity_penalty = sum(rules_per_node) * (current_gen / max_gen) * 0.002

            # Adjust the fitness by the complexity penalty
            adjusted_fitness = fitness_score + complexity_penalty

            fitnesses.append((adjusted_fitness,))

        return raw_fitnesses, fitnesses

    # 1.3 Generates a list of offspring, decides to do crossover or mutation
    def __varOrAdaptive(self, population, toolbox, genfrac, mutModel):
        """generates list of offspring to be compared... decides to do crossover or mutation"""
        # algorithm for generating a list of offspring... copied and pasted from DEAP with modification for adaptive mutation
        assert (self.crossover_probability + self.crossover_probability) <= 1.0, (
            "The sum of the crossover and mutation "
            "probabilities must be smaller or equal to 1.0."
        )
        offspring = []
        for _ in range(self.child_population_size):
            op_choice = random()
            if op_choice < self.crossover_probability:  # Apply crossover
                inds = []
                for samp in sample(population, 2):
                    ind = toolbox.clone(samp)
                    inds.append(ind)
                ind1, ind2 = inds
                ind1, ind2 = self.__cxTwoPointNode(ind1, ind2)
                del ind1.fitness.values
                offspring.append(ind1)
            elif op_choice < self.crossover_probability + self.mutation_probability:  # Apply mutation
                ind = toolbox.clone(choice(population))
                (ind,) = self.__mutFlipBitAdapt(ind, genfrac, mutModel)
                del ind.fitness.values
                offspring.append(ind)
            else:  # shouldn't happen... clone existing individual
                offspring.append(choice(population))
        return offspring

    # 1.3.1 Executes a two-point crossover between two individuals
    def __cxTwoPointNode(self, ind1, ind2):
        """Executes a two-point crossover on the input :term:`sequence`
        individuals. The two individuals are modified in place and both keep
        their original length.
        :returns: A tuple of two individuals.
        This function uses the :func:`~random.randint` function from the Python
        base :mod:`random` module.

        Modified from deap to cross over between rules = needed to account for bistring only being one of two components of individual
        """
        size = len(ind1[0].nodeList)
        cxpointer1 = randint(1, size)
        cxpointer2 = randint(1, size - 1)
        # make sure pointers are in right order
        if cxpointer2 >= cxpointer1:
            cxpointer2 += 1
        else:  # Swap the two cx points
            cxpointer1, cxpointer2 = cxpointer2, cxpointer1
        cxpoint1 = ind1[0].individualParse[cxpointer1]
        cxpoint2 = ind1[0].individualParse[cxpointer2]
        # cross over both bitlists and the total_combinationss (as well as total_inversionss)
        ind1[1][cxpoint1:cxpoint2], ind2[1][cxpoint1:cxpoint2] = (
            ind2[1][cxpoint1:cxpoint2],
            ind1[1][cxpoint1:cxpoint2],
        )
        (
            ind1[0].total_combinations[cxpointer1:cxpointer2],
            ind2[0].total_combinations[cxpointer1:cxpointer2],
        ) = (
            ind2[0].total_combinations[cxpointer1:cxpointer2],
            ind1[0].total_combinations[cxpointer1:cxpointer2],
        )
        (
            ind1[0].total_inversions[cxpointer1:cxpointer2],
            ind2[0].total_inversions[cxpointer1:cxpointer2],
        ) = (
            ind2[0].total_inversions[cxpointer1:cxpointer2],
            ind1[0].total_inversions[cxpointer1:cxpointer2],
        )
        # update the arrays seen by C code updateBool
        # ind1[0]._ruleMaker__updateCpointers()
        # ind2[0]._ruleMaker__updateCpointers()
        return ind1, ind2

    # 1.3.2 Mutation algorithm
    def __mutFlipBitAdapt(self, indyIn, genfrac, mutModel):
        """mutation algorithm
        :return: (indyIn,)
        """
        errors = list(indyIn.fitness.values)  # get errors of the individual
        individual = indyIn[1]  # Individual's data (bitstring)
        model = indyIn[0]  # Information about the individual
        # print(errors)

        # --- Get rid of errors in nodes that can't be changed ---
        errorNodes = 0
        for j in range(len(errors)):
            if model.node_combination_length[j] < 2:
                errors[j] = 0
            else:
                errorNodes = errorNodes + 1

        if np.sum(errors) < 0.05 * errorNodes or errorNodes == 0:
            # condition selection on number of incoming edges + downstream edges
            for i in range(len(model.nodeList)):
                # If there are no successors, pseudoerrors = length of possibilities for the node
                if model.successorNums[i] == 0:
                    pseudoerrors = [len(model.total_combinations[i])]
                # If there are successors, pseudoerrors = length of possibilities for the node * number of successors
                else:
                    pseudoerrors = [len(model.total_combinations[i]) * model.successorNums[i]]
            # zero out nodes that can't be changed
            for j in range(len(pseudoerrors)):
                if model.node_combination_length[j] < 2:
                    pseudoerrors[j] = 0
            focusNode = self.__selectMutNode(pseudoerrors)
        else:
            # if errors are relatively high, focus on nodes that fit the worst and have highest in-degree
            # calculate probabilities for mutating each node
            for i in range(len(errors)):
                temper = model.successorNums[i]
                if temper == 0:
                    errors[i] = errors[i] * len(model.total_combinations[i])
                else:
                    errors[i] = errors[i] * len(model.total_combinations[i]) * temper
            focusNode = self.__selectMutNode(errors)
        node = self.nodes[focusNode]

        # Perform mutation for the current node
        # Number of possible combinations for that node
        if model.node_combination_length[focusNode] >= 1:

            # find ends of the node of interest in the individual
            start = node.rule_start_index
            end = node.rule_end_index

            # mutate the inputs some of the time
            if len(model.total_combinations[focusNode]) > 3 and random() < mutModel:
                temppermup = []  # temporary upstream nodes
                upstreamAdders = list(model.total_combinations[focusNode])
                rvals = list(node.rvalues)
                while len(temppermup) < 3:
                    randy = random()  # randomly select a node to mutate
                    tempsum = sum(rvals)
                    if tempsum == 0:
                        addNoder = randint(
                            0, len(rvals) - 1
                        )  # int(math.floor(random()*len(upstreamAdders)))
                        # print(addNoder)
                    else:
                        recalc = np.cumsum([1.0 * rval / tempsum for rval in rvals])
                        # print(recalc)
                        addNoder = next(
                            i for i in range(len(recalc)) if recalc[i] > randy
                        )
                        # print(addNoder)
                    temppermup.append(upstreamAdders.pop(addNoder))
                    # print(rvals)
                    rvals.pop(addNoder)
                # model._ruleMaker__update_upstream(focusNode, temppermup)
            for i in range(start, end):
                # print("i: " + str(i))
                if random() < 2 / (end - start + 1):
                    individual[i] = 1
                else:
                    individual[i] = 0
            # ensure that there is at least one shadow and node turned on
            if np.sum(individual[start:end]) == 0:
                individual[start] = 1
            indyIn[0] = model
            indyIn[1] = individual

        else:
            next
            # msg = f'Did not perform mutation, node combination length < 1 ({model.node_combination_length[focusNode]})'
            # print(f'Node {[node.name for node in self.nodes if node.index == focusNode]} model.total_combinations:  {[model.total_combinations[focusNode]]}')
            # raise Exception(msg)
        return (indyIn,)

    # 1.3.2.1 Selects a node to mutate randomly, but nodes with higher errors have a greater chance of being selected
    def __selectMutNode(self, errors):
        """Select a node to mutate, the higher the error the greater the chance of being selected."""
        total_error = np.sum(errors)

        # Check if total_error is zero (or very close to zero)
        if np.isclose(total_error, 0):
            # Handle the case where total_error is zero
            # Example: select a node randomly
            selected_node = np.random.randint(0, len(errors))
        else:
            # Normalize errors to get a probability that the node is modified
            normerrors = [1.0 * error / total_error for error in errors]

            # Calculate cumulative probabilities
            probabilities = np.cumsum(normerrors)

            # Select the node
            selected_node = next(i for i in range(len(probabilities)) if probabilities[i] > random())

        return selected_node