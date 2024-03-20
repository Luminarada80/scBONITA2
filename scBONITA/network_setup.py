from deap import base, creator, gp, tools
from deap import algorithms as algo
import numpy as np
import networkx as nx
from sklearn import preprocessing
from scipy.stats.stats import spearmanr
import ctypes as ctypes
import itertools as itertool
import copy
import pickle
from random import random, randint, sample, choice
import math
from collections import defaultdict
from itertools import chain
from operator import attrgetter
import gc
import pandas as pd
import error_calculation
import time
import copy
from setup.parameters import Params
from parse_nodes import Node
import logging
from alive_progress import alive_bar

class NetworkSetup:
    def __init__(self, graph, gene_list, removeSelfEdges=False, restrictIncomingEdges=True,
                 maxIncomingEdges=3, groundTruth=False, graphName=""):
        """Initialize a NetworkSetup object for rule inference with scBONITA - RD"""

        # Initialize lists to store information about nodes and connections
        self.totalNodeList = []
        self.andNodeList = []
        self.andNodeInvertList = []
        self.andLenList = []
        self.totalLenList = []
        self.permList = []
        self.rvalues = []
        self.predecessors = []
        self.successorNums = []
        self.graph = graph
        self.params = Params()

        # Initialize node attributes
        self.nodeList = list(graph.nodes)  # List of nodes in the graph
        self.nodeDict = {self.nodeList[i]: i for i in range(len(self.nodeList))}  # Dictionary for node lookup

        self.combinations = {"three_incoming_nodes": {}, "two_incoming_nodes": {}, "one_incoming_nodes": {}}

        # Remove self-loops from the graph if required
        # self.remove_self_edges(removeSelfEdges, self.nodeList, graph)
        
        # Initialize an empty directed graph
        self.ruleGraph = nx.empty_graph(0, create_using=nx.DiGraph)

        # Calculate the node information
        self.calculate_node_information()

        # Create Node objects containing the calculated information for each node in the network
        self.nodes, self.deap_individual_length = self.create_nodes()

        # Print information about the network (optional)
        self.print_network_information()

        self.print_node_information()


    def print_node_information(self):
        for node in self.nodes:
            logging.debug(f'\nNode {node.name}:')
            logging.debug(f'\tpredecessors {node.predecessors}')
            logging.debug(f'\tthree_node_combos {node.three_node_combos}')
            logging.debug(f'\ttwo_node_AND_combos {node.two_node_AND_combos}')
            logging.debug(f'\ttwo_node_OR_combos {node.two_node_OR_combos}')
            logging.debug(f'\tone_node_combos {node.one_node_combos}')
            logging.debug(f'\tinversions {node.inversions}')
            logging.debug(f'\trule_length {node.rule_length}')
            logging.debug(f'\trule_start_index {node.rule_start_index}')
            logging.debug(f'\trule_end_index {node.rule_end_index}')

    def print_network_information(self):
        """Prints information about the network."""
        logging.debug(f'\nThree_incoming_nodes: {self.combinations["three_incoming_nodes"]}')
        logging.debug(f'\nTwo_incoming_nodes: {self.combinations["two_incoming_nodes"]}')
        logging.debug(f'\nOne_incoming_nodes: {self.combinations["one_incoming_nodes"]}')
        logging.debug("\nNodelist: " + str(self.nodeList))
        # logging.debug(f"\nsuccessorNums: {self.successorNums}")
        logging.debug("\nNode positions: " + str(self.node_positions))
        # logging.debug("\nPossibilityList: " + str(self.possibilityList))
        logging.debug(f"\ntotalNodeList: {self.totalNodeList}")

    # -------- Parse Node Information for Rule Inference --------
    # 1. Calculates the node information and creates Node objects storing the information for each node
    def calculate_node_information(self):
        """
        Calculates the information for each node in the network and stores the information as object of class Node
        from parse_nodes.py
        """
        # Iterate over all nodes to find predecessors and calculate possible connections
        for i, node in enumerate(self.nodeList):
            predecessors_final = self.find_predecessors(self.ruleGraph, self.nodeList, self.graph, self.nodeDict, [], i)
            possibilities = self.calculate_possible_connections(predecessors_final, self.nodeList)
            node_predecessors = [self.nodeList.index(corr_tuple[0]) for corr_tuple in predecessors_final]
            self.calculate_combinations(self.permList, i, possibilities)
            self.predecessors.append(node_predecessors)

    # 1.1 Removes self-loops from the nodes in the graph
    def remove_self_edges(self, removeSelfEdges, nodeList, graph):
        """
        Remove self loops from the graph
        """
        if removeSelfEdges:
            for node in nodeList:
                repeat = True
                while repeat:
                    repeat = False
                    if node in list(graph.successors(node)):
                        graph.remove_edge(node, node)
                        repeat = True

    # 1.2 Finds the predecessors of each node and stores the top 3
    def find_predecessors(self, ruleGraph, nodeList, graph, nodeDict, possibilityLister, node):
        """
        Find the incoming nodes for each ndoe in the graph, store the top 3 connections as calculated by a spearman
        correlation
        Parameters
        ----------
        ruleGraph
        nodeList
        graph
        nodeDict
        possibilityLister
        node

        Returns
        -------

        """
        # --- Find the predecessors of each node ---
        predecessors_final, predCorr_temp = self.parse_connections(node, nodeList, graph, nodeDict, possibilityLister)

        # Store the correlations between incoming nodes in "rvalues"
        top_three_incoming_node_correlations = sorted(predCorr_temp, reverse=True)[:3]
        self.rvalues.append(top_three_incoming_node_correlations)

        # Append the permanent list with the top 3 predecessors for this node
        self.permList.append([pred[0] for pred in predecessors_final])

        # Add the incoming nodes and their properties to the newly created ruleGraph
        self.store_predecessors(graph, nodeList, predecessors_final, ruleGraph, node)
        return predecessors_final

    # 1.2.1 Finds the top 3 incoming nodes using Spearman correlation for the current node and stores that info
    def parse_connections(self, node, nodeList, graph, nodeDict, possibilityLister):
        """
        Find the top 3 incoming nodes for the node of interest and store the information (part of find_predecessors)
        """
        # Get NAMES of predecessors and successors of the node from original graph
        
        predecessors_temp = list(graph.predecessors(nodeList[node]))
        
        successors_temp = list(graph.successors(nodeList[node]))
        self.successorNums.append(len(successors_temp))

        possibilitytemp = [nodeDict[predder] for predder in predecessors_temp]
        possibilityLister.append(list(possibilitytemp))

        # Calculate the Spearman correlation for the incoming nodes
        predCorr_temp = self.calculate_spearman_correlation(node, predecessors_temp)

        # Select the top 3 predecessors of the node
        predecessors_final = sorted(
            zip(predecessors_temp, predCorr_temp),
            reverse=True,
            key=lambda corrs: corrs[1], )[:3]
        return predecessors_final, predCorr_temp
    
    def calculate_spearman_correlation(self, node, predecessors_temp):
        """
        Calculate the Spearman correlation between incoming nodes to find the top three with the
        highest correlation, used to reduce the dimensionality of the calculations.

        1. calculate_node_information
            1.2 find_predecessors
                1.2.1 parse_connections
                    1.2.1.1 calculate_spearman_correlation
        """
        # Find correlation between the predecessors and the node
        nodeData = (
            self.binarized_matrix[self.node_positions[node], :].todense().tolist()[0]
        )  # find binarized expression data for node "i"
        predCorr_temp = (
            []
        )  # temporarily store correlations between node "i" and all its predecessors

        for predecessor_gene in predecessors_temp:
            # find index of predecessor in the gene_list from the data
            predIndex = self.nodeList.index(predecessor_gene)

            # find binarized expression data for predecessor
            predData = (self.binarized_matrix[predIndex, :].todense().tolist()[0])
            mi, pvalue = spearmanr(nodeData, predData)
            if np.isnan(mi):
                predCorr_temp.append(0)
            else:
                predCorr_temp.append(mi)  # store the calculated correlation
        return predCorr_temp


    # 1.2.2 Stores the top 3 interactions in a new graph
    def store_predecessors(self, graph, nodeList, predecessors_final, ruleGraph, node):
        """
        Stores the information about each incoming node for each node in the graph to a new
        graph.
        """
        for parent in predecessors_final:
            if "interaction" in list(graph[parent[0]][nodeList[node]].keys()):
                ruleGraph.add_edge(
                    parent[0],
                    nodeList[node],
                    weight=parent[1],
                    activity=graph[parent[0]][nodeList[node]]["interaction"],
                )
            if "signal" in list(graph[parent[0]][nodeList[node]].keys()):
                ruleGraph.add_edge(
                    parent[0],
                    nodeList[node],
                    weight=parent[1],
                    activity=graph[parent[0]][nodeList[node]]["signal"],
                )

    # 1.3 Calculates the possible connections for the incoming nodes
    def calculate_possible_connections(self, predecessors_final, nodeList):
        """
        Calculate the possible combinations of incoming nodes
        """
        # Write the node names from predecessors_final and repeat('empty')
        withNones = zip(
            [nodeList.index(corr_tuple[0]) for corr_tuple in predecessors_final],
            itertool.repeat("empty"),
        )

        # Iterate through every possibility for the node, with 'empty' as the filler for 2 or 1 incoming AND node possibilites
        possibilities = list(itertool.product(*withNones))  # Make a list of 0's

        # Trim out any of the "empty" values
        for j in range(0, len(possibilities)):
            possibilities[j] = list(possibilities[j])
            while "empty" in possibilities[j]:
                possibilities[j].remove("empty")
            while [] in possibilities[j]:
                possibilities[j].remove([])
        while [] in possibilities:
            possibilities.remove([])

        return possibilities

    # 1.4 Calculates the possible combinations of three, two, and one incoming nodes for determining possible rules
    def calculate_combinations(self, permList, node, possibilities):
        """
        Takes in the possible combinations for each node and stores them in a combinations dictionary.

        For nodes that have three possible incoming nodes, this stores the possibilities for having signaling from all
        three incoming nodes, two of the three nodes, or just one of the three nodes.

        For nodes that only have two possible incoming nodes, this stores the possibilities for having signaling from
        either both nodes or just one of the two nodes.

        For nodes that only have one possible incoming node, this stores that possibility.
        """

        """
        combination = {
            "three_incoming_nodes": 
                {1: [], 2: [], 3: []}
            "two_incoming_nodes":
                {1: [], 2: [], 3: []}
            "one_incoming_nodes":
                {1: [], 2: [], 3: []}
        }
        """

        three_node_possibilities = []
        two_node_possibilities = []
        one_node_possibilities = []

        # Three possible incoming nodes
        if len(permList[node]) == 3:
            for possibility in possibilities:
                if len(possibility) == 3:
                    three_node_possibilities.append(possibility)

                # Add the two incoming node possibilities to a list
                elif len(possibility) == 2:
                    two_node_possibilities.append(possibility)

                # Add the one incoming node possibilities to a list
                elif len(possibility) == 1:
                    one_node_possibilities.append(possibility)

            self.combinations["three_incoming_nodes"][
                node] = three_node_possibilities + two_node_possibilities + one_node_possibilities
            self.combinations["two_incoming_nodes"][node] = two_node_possibilities
            self.combinations["one_incoming_nodes"][node] = one_node_possibilities

        # Two possible incoming nodes
        elif len(permList[node]) == 2:
            for possibility in possibilities:
                # Add the two incoming node possibilities to a list
                if len(possibility) == 2:
                    two_node_possibilities.append(possibility)

                # Add the one incoming node possibilities to a list
                elif len(possibility) == 1:
                    one_node_possibilities.append(possibility)
            self.combinations["two_incoming_nodes"][node] = two_node_possibilities
            self.combinations["one_incoming_nodes"][node] = one_node_possibilities

        # One possible incoming node
        elif len(permList[node]) == 1:
            for possibility in possibilities:
                if len(possibility) == 1:
                    one_node_possibilities.append(possibility)
            self.combinations["one_incoming_nodes"][node] = one_node_possibilities

    # 1.5.1 Calculates the inversion rules for each combination
    def calculate_inversion_rules(self, node):
        """
        Calculates the inversion rules for a node based on the graph interactions or signal for each incoming node
        Parameters
        ----------
        node

        Returns
        -------
        inversion_rules
        """

        inversion_rules = {}
        for incoming_node in list(node.predecessors.keys()):
            edge_attribute = list(self.graph[self.nodeList[incoming_node]][self.nodeList[node.index]].keys())
            # check the 'interaction' edge attribute
            if "interaction" in edge_attribute:
                if self.graph[self.nodeList[incoming_node]][self.nodeList[node.index]]["interaction"] == "i":
                    inversion_rules[incoming_node] = True
                else:
                    inversion_rules[incoming_node] = False

            # check the 'signal' edge attribute
            elif "signal" in edge_attribute:
                if self.graph[self.nodeList[incoming_node]][self.nodeList[node.index]]["signal"] == "i":
                    inversion_rules[incoming_node] = True
                else:
                    inversion_rules[incoming_node] = False
        return inversion_rules

    # 1.6 Creates nodes containing the information calculated from the graph
    def create_nodes(self):
        """
        Creates Node class objects using the information calculated in the rest of the calculate_node_information
        function
        """
        inverted_nodeDict = {v: k for k, v in self.nodeDict.items()}
        gene_name_to_index = {gene_name: gene_index for gene_index, gene_name in enumerate(self.gene_names)}


        nodes = []
        rule_index = 0
        with alive_bar(len(self.nodeList)) as bar:
            for node_index, node_name in enumerate(self.nodeList):
                name = node_name

                # Use .get() with default value [] for missing combinations or inversions
                three_node_combos = self.combinations.get("three_incoming_nodes", {}).get(node_index, [])
                two_node_combos = self.combinations.get("two_incoming_nodes", {}).get(node_index, [])
                one_node_combos = self.combinations.get("one_incoming_nodes", {}).get(node_index, [])

                # Safely retrieve predecessors and put them into a dictionary where key = node index, value = node name
                predecessor_indices = self.predecessors[node_index] if node_index < len(self.predecessors) else []
                predecessors = {}
                for index in predecessor_indices:
                    inverted_nodeDict = {v: k for k, v in self.nodeDict.items()}
                    predecessors[index] = inverted_nodeDict[index]
                
                # Create a new Node object
                node = Node(name, node_index, three_node_combos, two_node_combos, one_node_combos, predecessors)

                # Find the dataset row index of the gene
                node.dataset_index = gene_name_to_index.get(node_name)

                node.rvalues = self.rvalues[node_index]
                node.inversions = self.calculate_inversion_rules(node)

                # Find the start and end indices for where the rule combinations start and stop for this node in the
                # bitstring of the individual (1110010100010010101001110101 etc.)
                rule_length = node.find_rule_length()
                if rule_length > 0:
                    node.rule_start_index = rule_index
                    rule_index += rule_length
                    node.rule_end_index = rule_index
                else:
                    node.rule_start_index = None
                    node.rule_end_index = None
                nodes.append(node)
                bar()
        return nodes, rule_index

    # -------- Find the Best Ruleset --------
    # 1. Find the lowest error individual in a population
    def __findPopBest(self, population):
        """finds the lowest error individual in a population"""
        saveVal = -1
        minny = float("Inf")
        for i in range(len(population)):
            if np.sum(population[i].fitness.values) < minny:
                minny = np.sum(population[i].fitness.values)
                saveVal = i
        ultimate = population[saveVal]
        minvals = population[saveVal].fitness.values
        return minvals, ultimate[1], ultimate[0]

    # -------- Checking Node Possibilities --------
    # 1. Evaluates rule possibilities for a node in a model and identifies rule with lowest error
    def __checkNodePossibilities(self, node, individual, knockout_list, knockin_list, sc_sync_bool_count):
        """
        Evaluates different rule possibilities for a given node in a model and identifies the rule that yields
        the minimum error, along with other rules that have errors within a certain tolerance of the minimum error.

        The function explores all possible combinations of rules for the specified node, calculates the errors
        associated with each rule, and then returns the best rule and a set of equivalent rules.

        **Parameters**
         - node: The node for which the rule possibilities are to be evaluated.
         - individual: A list representing the current set of rules for all nodes in the model.
         - knockout_list: A list representing knockout conditions for the model.
         - knockin_list: A list representing knockin conditions for the model.
         - sc_sync_bool_count: A parameter used in the evaluation of the model.

        **Returns**
         - best_rule: A list representing the best rule found for the node.
         - equivalent_rules: A list of rules that are equivalent to the best rule within a certain tolerance.
         - min_error: The minimum error associated with the best rule.
         - rule_errors: A list of errors associated with each rule possibility evaluated.

        The function is part of a genetic algorithm optimization process and is used to refine the rules
        associated with each node in the model.
        """
        # Set the model to the current instance
        model = self

        # Define a tolerance level for considering rules as equivalent
        tolerance = 0.0

        # Find the start and end indices of the rules corresponding to the given node
        end_index = model._NetworkSetup__findEnd(node)
        start_index = model.individualParse[node]

        # Extract the current rule for the node from the individual
        current_rule = list(individual[start_index:end_index])

        # Initialize a list to store equivalent rules, starting with the current rule
        equivalent_rules = [current_rule]

        # If there are no rules for the node, return early
        if (end_index - start_index) == 0:
            return current_rule, equivalent_rules, equivalent_rules, 0.0

        # Initialize lists to store rule options and their corresponding errors
        rule_options = []
        rule_errors = []

        # Iterate over all possible rule combinations for this node
        for i in range(1, 2 ** (end_index - start_index)):
            # Create a copy of the current individual
            temp_individual = list(individual)

            # Generate a binary representation of the rule being checked
            temp_rule = model._NetworkSetup__bitList(i, len(current_rule))

            # Replace the current rule with the new rule in the temporary individual
            temp_individual[start_index:end_index] = temp_rule

            # Evaluate the error for the temporary individual
            node_errors_temp = self._NetworkSetup__evaluateByNode(
                temp_individual, knockout_list, knockin_list, sc_sync_bool_count, localSearch=True
            )

            # Extract the error corresponding to the node being checked
            current_error = node_errors_temp[node]

            # Store the rule and its error
            rule_options.append(temp_rule)
            rule_errors.append(current_error)

            # Clear the memory
            gc.collect()

        # Find the minimum error among all rules
        min_error = min(rule_errors)
        equivalent_rules = []

        # Identify rules that are equivalent to the best rule within the tolerance
        for i in range(len(rule_options)):
            if rule_errors[i] <= min_error + tolerance:
                equivalent_rules.append(rule_options[i])

        # The best rule is the first one in the list of equivalent rules
        best_rule = equivalent_rules[0]

        # Return the best rule, equivalent rules, minimum error, and all rule errors
        return (best_rule, equivalent_rules, min_error, rule_errors)

    # 1.1 Finds where the current node ends in invalid_ind
    def __findEnd(self, node):
        """Finds where the current node ends in the bitstring of an individual"""
        if node == len(self.nodeList) - 1:
            end = self.size
        else:
            end = self.individualParse[node + 1]
        return end

    # 1.2 Convert integer to a binary string, then to a reversed list of 0s and 1s
    def __bitList(self, number, desired_length):
        # Convert the integer to a binary string and then to a reversed list of 0s and 1s
        binary_digits = [1 if digit == "1" else 0 for digit in bin(number)[::-1]]

        # Ensure the list of binary digits has the desired length
        while len(binary_digits) < desired_length:
            binary_digits.append(0)
        while len(binary_digits) > desired_length:
            binary_digits.pop()

        return binary_digits

    # 1.3 __evaluateByNode


    # -------- Calculating Importance Scores --------

    # 1. Calculates Importance Scores
    def __calcImportance(self, equivs, model, importanceScore, graphName):
        # Create holder for importance scores
        importanceScoresDict = {}
        importanceScoreStdev = {}


        temp_list = []
        # List of node indices
        for node in self.nodes:
            temp_list.append(node.index)

        # Set up the importance scores for each node index
        for node in range(0, len(self.nodeList)):
            importanceScoresDict[self.nodeList[node]] = []
            importanceScoreStdev[self.nodeList[node]] = 0.0

        

        # Try 3 randomly sampled rule sets
        i = 0
        while i < 3:

            # Find an ERS
            individual = self._NetworkSetup__processERS(graphName + "_equivs1.pickle")

            # For each node
            for node in temp_list:

                # Print the node and its position
                print(
                    "Node: "
                    + str(self.nodeList[node])
                    + ", Node Position: "
                    + str(node)
                )

                # Evaluate the importance scores
                temp = self._NetworkSetup__evaluateByNode(
                    individual,
                    [node],
                    [node],
                    importanceScore,
                    localSearch=False,
                    importanceScores=True,
                )

                # Append the node's importance score to the dictionary
                print("Trial: " + str(i) + " Unprocessed IS: " + str(temp))
                importanceScoresDict[self.nodeList[node]].append(temp[0])
            i = i + 1

        print(importanceScoresDict)

        # Find maximum node importance score
        maxScore = max(importanceScoresDict.values())
        print("Max IS: " + str(maxScore))
        minScore = min(importanceScoresDict.values())
        print("Min IS: " + str(maxScore))

        # Rescaling to [0,1] using featureReScale
        for node in range(0, len(self.nodeList)):
            importanceScoresDict[self.nodeList[node]] = (
                importanceScoresDict[self.nodeList[node]][0] - minScore[0]
            ) / (maxScore[0] - minScore[0])

        print(importanceScoreStdev)

        ersFile = open(str(graphName + "_equivs1.pickle"), "rb")
        ers = pickle.load(ersFile)
        obsERS = {}
        maxERS = {}
        inDegreeNet = nx.read_graphml(graphName)
        # Normalize by number of rule sets that were tried
        for node in range(0, len(self.nodeList)):
            obsERS[self.nodeList[node]] = len(ers[node])
            inDegree = inDegreeNet.in_degree(self.nodeList[node])
            if inDegree == 0:
                maxERS[self.nodeList[node]] = 1
            else:
                inDegree = min(inDegree, 3)
                maxERS[self.nodeList[node]] = (
                    2 ** (len(ers[node][0])) - 1
                )  # 2**(inDegree+1) - 1 #
            importanceScoresDict[self.nodeList[node]] = np.mean(
                importanceScoresDict[self.nodeList[node]]
            )
            importanceScoresDict[self.nodeList[node]] = importanceScoresDict[
                self.nodeList[node]
            ] * (
                (maxERS[self.nodeList[node]] - obsERS[self.nodeList[node]] + 1)
                / maxERS[self.nodeList[node]]
            )

        # Print out file of importance scores
        IS_df = pd.DataFrame(
            importanceScoresDict.items(), columns=["Node", "Importance Score"]
        )
        IS_df["ObsERS"] = IS_df["Node"].map(obsERS)
        IS_df["MaxERS"] = IS_df["Node"].map(maxERS)

        IS_df.to_csv(
            str(graphName + "_importanceScores.csv"),
            sep=",",
            encoding="utf-8",
            index=False,
        )

        # Make graphml with importance scores as attributes
        net = self.ruleGraph
        nx.set_node_attributes(net, values=importanceScoresDict, name="importanceScore")
        nx.set_node_attributes(net, values=maxERS, name="maxERS")
        nx.set_node_attributes(net, values=obsERS, name="obsERS")

        # add abundance as attribute to graph
        binarized_matrix2 = self.binarized_matrix.A
        abundance = {}
        abundance_sd = {}
        numZeros = {}
        numOnes = {}

        for node in list(importanceScoresDict.keys()):
            node_index = self.gene_list.index(node)
            expression = binarized_matrix2[node_index, :].tolist()
            abundance[node] = np.mean(expression)
            abundance_sd[node] = np.std(expression)
            expression = np.array(expression)
            numZeros[node] = (expression == 0).sum()
            numOnes[node] = (expression == 1).sum()

        nx.set_node_attributes(net, values=abundance, name="abundanceMean")
        nx.set_node_attributes(net, values=abundance_sd, name="abundanceStdev")
        nx.set_node_attributes(net, values=numZeros, name="abundanceZeros")
        nx.set_node_attributes(net, values=numOnes, name="abundanceOnes")

        nx.write_graphml_lxml(net, str(graphName) + "_IS.graphml")

        return importanceScoresDict


    # 1.1 Calculates equavalent rulesets
    def __processERS(self, equivsName):
        """Create an individual from the ERS generated by the local search, for importance score calculation"""
        ersFile = open(str(equivsName), "rb")
        ers = pickle.load(ersFile)
        ersFile.close()

        # randomly sample the ers to make an individual
        individual = []
        for i in range(len(ers)):
            individual.extend(ers[i][randint(0, len(ers[i]) - 1)])
        return individual

    # 1.2 Simulates KO and KI through the network to calculate node importance
    def __evaluateByNode(
        self,
        individual,
        KOlist,
        KIlist,
        cFunction,
        localSearch=False,
        importanceScores=False,
    ):
        """Includes Network Propagation"""
        model = self
        cellArray = []
        knockins = np.zeros(len(model.nodeList), dtype=np.intc, order="C")
        knockouts = np.zeros(len(model.nodeList), dtype=np.intc, order="C")

        for knocker in KOlist:
            knockouts[knocker] = 1
        for knocker in KOlist:
            knockins[knocker] = 1

        # put objects in correct format for passing to C
        nodeIndividual = np.array(individual, dtype=np.intc, order="C")
        indLen = len(nodeIndividual)
        andNodes = np.array(model.andNodes, dtype=np.intc, order="C")
        nodeNum = len(model.nodeList)
        andNodeInvert = np.array(model.andNodeInvert, dtype=object, order="C")
        individualParse = np.array(model.individualParse, dtype=np.intc, order="C")
        andLenList = np.array(model.andLenList, dtype=np.intc, order="C")
        node_positions1 = model.node_positions
        node_positionsC = np.array(node_positions1, dtype=np.intc, order="C")
        simSteps = self.params.simSteps
        lenSamples1 = len(model.sampleList)
        binarized_matrixC3 = np.array(
            copy.deepcopy(self.binarized_matrix.toarray(order="C")), order="C", dtype=np.intc
        )
        binarized_matrixCPointer = ctypes.c_void_p(
            binarized_matrixC3.ctypes.data
        )  # put input array as C pointer

        # convert objects into C pointers
        nodeIndividual1 = ctypes.c_void_p(nodeIndividual.ctypes.data)
        indLen1 = ctypes.c_void_p(indLen)
        andNodes1 = ctypes.c_void_p(andNodes.ctypes.data)
        individualParse1 = ctypes.c_void_p(individualParse.ctypes.data)
        andLenList1 = ctypes.c_void_p(andLenList.ctypes.data)
        andNodeInvertList1 = ctypes.c_void_p(andNodeInvert.ctypes.data)
        nodeNum1 = ctypes.c_void_p(nodeNum)

        simSteps1 = ctypes.c_void_p(simSteps)
        knockouts1 = ctypes.c_void_p(knockouts.ctypes.data)
        knockins1 = ctypes.c_void_p(knockins.ctypes.data)

        node_positionsCPointer = ctypes.c_void_p(node_positionsC.ctypes.data)

        vals = np.full(
            shape=(100, self.maxNodes), fill_value=0, dtype=np.intc, order="C"
        )  # simData[STEP][NODE]
        valsubmit = ctypes.c_void_p(vals.ctypes.data)
        lenSamples = ctypes.c_void_p(lenSamples1)

        localSearchC = ctypes.c_void_p(int(localSearch))
        importanceScoresC = ctypes.c_void_p(int(importanceScores))

        # errors = np.array(np.full(10000, fill_value=0, dtype=np.intc, order='C'))
        # errorsSubmit=ctypes.c_void_p(errors.ctypes.data)

        if localSearch:
            # look at errors node wise
            errors = np.array(
                np.full(self.maxNodes, fill_value=0, dtype=np.intc, order="C")
            )
            errorsSubmit = ctypes.c_void_p(errors.ctypes.data)
            cFunction(
                valsubmit,
                nodeIndividual1,
                indLen1,
                nodeNum1,
                andLenList1,
                individualParse1,
                andNodes1,
                andNodeInvertList1,
                simSteps1,
                knockouts1,
                knockins1,
                lenSamples,
                binarized_matrixCPointer,
                node_positionsCPointer,
                errorsSubmit,
                localSearchC,
                importanceScoresC,
            )  # in this case scSyncBoolC
            errors = errors.tolist()
            errors = errors[:nodeNum]
            return errors
        else:
            if importanceScores:
                importanceScores = np.array(
                    np.full(1, fill_value=0.0, dtype=np.float64, order="C")
                )
                importanceScoresC = ctypes.c_void_p(importanceScores.ctypes.data)
                cFunction(
                    valsubmit,
                    nodeIndividual1,
                    indLen1,
                    nodeNum1,
                    andLenList1,
                    individualParse1,
                    andNodes1,
                    andNodeInvertList1,
                    simSteps1,
                    knockouts1,
                    knockins1,
                    lenSamples,
                    binarized_matrixCPointer,
                    node_positionsCPointer,
                    importanceScoresC,
                )  # in this case importanceScore
                return importanceScores.tolist()
            else:
                # look at errors by sample
                errors = np.array(
                    np.full(self.maxSamples, fill_value=0, dtype=np.intc, order="C")
                )
                errorsSubmit = ctypes.c_void_p(errors.ctypes.data)
                cFunction(
                    valsubmit,
                    nodeIndividual1,
                    indLen1,
                    nodeNum1,
                    andLenList1,
                    individualParse1,
                    andNodes1,
                    andNodeInvertList1,
                    simSteps1,
                    knockouts1,
                    knockins1,
                    lenSamples,
                    binarized_matrixCPointer,
                    node_positionsCPointer,
                    errorsSubmit,
                    localSearchC,
                    importanceScoresC,
                )  # in this case scSyncBoolC
                errors = errors.tolist()
                return [sum(errors)]

