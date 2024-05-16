from itertools import combinations, product
import random
import numpy as np
from parse_node import Node

A = 1
B = 0
C = 1

incoming_nodes = ["A", "B", "C"]

def enumerate_possibilities(incoming_nodes):
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

    num_incoming_nodes = len(incoming_nodes)
    if num_incoming_nodes == 3:
        possibilities = np.array(boolCombo("ABC") + boolCombo("AB") + boolCombo("AC") + boolCombo("BC") + ["A"] + ["B"] + ["C"])
    elif num_incoming_nodes == 2:
        possibilities = np.array(boolCombo("AB") + ["A"] + ["B"])
    elif num_incoming_nodes == 1:
        possibilities = np.array(["A"])
    else:
        assert IndexError('Num incoming nodes out of range')
    return possibilities

def choose_rule(possibilities):
    random_rule_index = random.choice(range(len(possibilities)))
    # print(f'Random rule index: {random_rule_index}\n')

    bitstring = np.array([1 if i == random_rule_index else 0 for i, _ in enumerate(possibilities)])

    selected_rule = possibilities[bitstring == 1][0].replace("  "," ")

    return bitstring, selected_rule

possibilities = enumerate_possibilities(incoming_nodes)
bitstring, selected_rule = choose_rule(possibilities)

# print(f'Bitstring: {bitstring}\n')
# print(f'Selected rule: {selected_rule}\n')

# for i, combo in enumerate(possibilities):
#     clean_combo = combo.replace("  "," ")
#     print(i,clean_combo)

# print(f'Evaluation:')
# print(f'\tA = {A}')
# print(f'\tB = {B}')
# print(f'\tC = {C}')
# print(f'Eval {selected_rule} = {eval(selected_rule)}')

predecessors = {1:'pred1',
                2:'pred2',
                3:'pred3'}

inversions = [1, 0, 1]

test_node = Node('test_node', 0, predecessors, inversions)

print(f'Possible rules: {possibilities}')
print(f'Bitstring: {test_node.bitstring}\n')
print(f'Selected rule: {test_node.selected_rule}\n')

# for i, combo in enumerate(test_node.possibilities):
#     clean_combo = combo.replace("  "," ")
#     print(i,clean_combo)

print(f'Evaluation:')
print(f'\tA = {A}')
print(f'\tB = {B}')
print(f'\tC = {C}')
print(f'Eval {test_node.selected_rule} = {eval(test_node.selected_rule)}')

print(f'Node Rule: {test_node.node_rules}')