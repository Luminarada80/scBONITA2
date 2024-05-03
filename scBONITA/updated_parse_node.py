import logging
from boolean_rule_functions import *
from itertools import combinations, product
import numpy as np
import random

class Node:
    def __init__(self,
                 name,
                 index,
                 predecessors,
                 inversions):

        self.name = name # Gene name
        self.index = index # Index in the network
        self.dataset_index = None

        # Combining all of the combinations
        self.rvalues = None

        # Inversion rules for each of the incoming node combinations
        self.inversions = inversions

        # Upstream nodes signaling to this node (dictionary where keys = node index, value = node name)
        self.predecessors = predecessors

        # Finds the possible rules based on the number of predecessors and chooses one using one-hot encoding
        self.possibilities = self.enumerate_possibilities()
        self.bitstring, self.selected_rule = self.choose_rule()
        self.node_rules = [self.name, [i for i in predecessors], self.selected_rule, self.inversions]

        # Indices in the individuals where the rules for this node start and end
        self.rule_start_index = None
        self.rule_end_index = None

        # Finds the best rule
        self.best_rule_index = None
        self.best_rule = None
        self.best_rule_error = None

        self.logic_function = None

        # Compressed best rule information
        self.calculation_function = None

        # Importance scores
        self.importance_score = 0
        self.importance_score_stdev = 0.0

        # Relative abundance
        self.relative_abundance = None
    
    def enumerate_possibilities(self):
        cache = dict()
        or0,or1,and0,and1 = "  or  "," or ","  and  "," and "  
        
        def cacheResult(keys,result=None):
            if not result:
                return [ r.format(*keys) for r in cache.get(len(keys),[]) ]   
            cache[len(keys)] = resFormat = []
            result = sorted(result,key=lambda x:x.replace("  "," "))
            for r in result:
                r = r.replace("and","&").replace("or","|")
                for i,k in enumerate(keys):
                    r = r.replace(k,"{"+str(i)+"}")
                r = r.replace("&","and").replace("|","or")
                resFormat.append(r)
            return result

        def boolCombo(keys):
            if len(keys)==1: 
                return list(keys)
            
            result = cacheResult(keys) or set()
            if result: 
                return result
            
            def addResult(left,right):
                OR = or0.join(sorted(left.split(or0)+right.split(or0)))
                result.add(OR.replace(and0,and1))
                if or0 in left:  
                    left  = f"({left})"
                if or0 in right: 
                    right = f"({right})"
                AND = and0.join(sorted(left.split(and0)+right.split(and0)))
                result.add(AND.replace(or0,or1))
                    
            seenLeft  = set()
            for leftSize in range(1,len(keys)//2+1):
                for leftKeys in combinations(keys,leftSize):
                    rightKeys = [k for k in keys if k not in leftKeys]
                    if len(leftKeys)==len(rightKeys):
                        if tuple(rightKeys) in seenLeft: continue
                        seenLeft.add(tuple(leftKeys))
                    for left,right in product(*map(boolCombo,(leftKeys,rightKeys))):
                        addResult(left,right)
            return cacheResult(keys,result)

        possibilities = None
        num_incoming_nodes = len(self.predecessors)
        if num_incoming_nodes == 3:
            possibilities = np.array(boolCombo("ABC") + boolCombo("AB") + boolCombo("AC") + boolCombo("BC") + ["A"] + ["B"] + ["C"])
        elif num_incoming_nodes == 2:
            possibilities = np.array(boolCombo("AB") + ["A"] + ["B"])
        elif num_incoming_nodes == 1:
            possibilities = np.array(["A"])
        else:
            assert IndexError('Num incoming nodes out of range')

        return possibilities

    def choose_rule(self):
        random_rule_index = random.choice(range(len(self.possibilities)))
        # print(f'Random rule index: {random_rule_index}\n')

        # Creates a bitstring for the rule
        bitstring = np.array([1 if i == random_rule_index else 0 for i, _ in enumerate(self.possibilities)])

        # Finds the predicted rule from the bitstring
        selected_rule = self.possibilities[bitstring == 1][0].replace("  "," ")

        # Formats the rule with not values
        for i, invert_node in enumerate(self.inversions):
            if invert_node == 1:
                if i == 0:
                    selected_rule = selected_rule.replace('A', 'not A')
                if i == 1:
                    selected_rule = selected_rule.replace('B', 'not B')
                elif i == 2:
                    selected_rule = selected_rule.replace('C', 'not C')

        return bitstring, selected_rule
    
    def find_multiple_rule_predictions(self, individual_bitstring):
        
        rule_predictions = []

        # Extract the bitstring for this node from the individual
        bitstring_length = self.rule_end_index - self.rule_start_index

        if bitstring_length >= 1:
            bitstring = np.array(individual_bitstring[self.rule_start_index:self.rule_end_index])

        elif bitstring_length == 0:
            bitstring = np.array(individual_bitstring[self.rule_start_index])
        
        # Replace all 1 values in the bitstring with the rule name
        selected_rules = self.possibilities[bitstring == 1][0].replace("  "," ")

        # Add in the inversion rules
        for rule in selected_rules:
            for i, invert_node in enumerate(self.inversions):
                if invert_node == 1:
                    if i == 0:
                        rule = rule.replace('A', 'not A')
                    if i == 1:
                        rule = rule.replace('B', 'not B')
                    elif i == 2:
                        rule = rule.replace('C', 'not C')
            rule_predictions.append(rule)

        # Update the predicted rules for this node
        for rule in rule_predictions:
            self.node_rules = [self.name, [i for i in self.predecessors], rule, self.inversions]
        
        return rule_predictions


    # Print information about the node
    def print_info(self):
        logging.info(f'\nNode {self.name} Length: {self.rule_length}')
        logging.info(f'\tindex {self.index}')
        logging.info(f'\tpredecessors: {self.predecessors}')
        logging.info(f'\tinversions: {self.inversions}')
        logging.info(f'\tNode rules: {self.node_rules}')

    # Reset the state of the node
    def reset_state(self):
        """
        Reset the state of the Node instance to clear any data that is specific to an individual's rule set.
        This method should be called before processing a new individual.
        """
        self.rule_predictions = None
        self.node_rules = []  # Clearing the list of rules

        # Any other attributes that need to be reset can be added here



