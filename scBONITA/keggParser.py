from bioservices.kegg import KEGG
import networkx as nx
import requests
import re
import string
from bs4 import BeautifulSoup
import itertools
from rule_inference import *


def parse_kegg_dict():
    """makes a dictionary to convert ko numbers from KEGG into real gene names
    #this is all file formatting. it reads a line, parses the string into the gene name and ko # then adds to a dict that identifies the two."""
    converterDict = {}
    pathwayFile = requests.get("http://rest.kegg.jp/get/br:ko00001", stream=True)
    for line in pathwayFile.iter_lines():
        line = line.decode("utf-8")
        if (
            len(line) > 1 and line[0] == "D"
        ):  # lines which begin with D translate kegg codes to gene names
            # to split into kegg code, gene names
            converter = re.split(r"\s+", re.split(r";", line)[0])
            converterDict[converter[1].upper()] = converter[2].upper()
    return converterDict


def deconvolute_groups(node_id, groups):
    """
    Creates a list with each group containing the node

    node_id: a node ID that may be a group
    groups: store group IDs and list of sub-ids
    return value: a list that contains all group IDs deconvoluted
    """
    node_list = []
    if node_id in groups.keys():
        # Deconvolute each of the groups that the node is in
        for component_id in groups[node_id]:
            node_list.extend(deconvolute_groups(component_id, groups))
    else:
        node_list.extend([node_id])
    return node_list

def read_kegg(lines, graph, KEGGdict, hsaDict):
    # read all lines into a bs4 object using libXML parser
    soup = BeautifulSoup("".join(lines), "xml")
    groups = {}  # store group IDs and list of sub-ids
    id_to_name = {}  # map id numbers to names

    # Look at each entry in the kgml file. Info: (https://www.kegg.jp/kegg/xml/)
    for entry in soup.find_all("entry"):
        print(f'entry: {entry}\n')

        gene_names_in_entry = entry["name"]
        print(f'\tgene_names_in_entry: {gene_names_in_entry}')

        # Name of each gene in the entry
        entry_split = gene_names_in_entry.split(":")
        print(f'\tentry_split: {entry_split}')
        print(f'\tlen(entry_split) : {len(entry_split)}')
        # If the entry is part of a group (in the network coded by a group containing lots of related genes)
        if len(entry_split) > 2:

            # Choose which dictionary to use based on whether the entries are hsa or kegg elements
            if entry_split[0] == "hsa" or entry_split[0] == "ko":
                if entry_split[0] == "hsa":
                    database = hsaDict
                else:
                    database = KEGGdict

                # Remove hsa or ko from the name of the genes in the group to isolate the gene ID
                gene_number = entry_split.pop(0)
                gene_number = gene_number.split()[0]
                print(f'gene_number: {gene_number}')

                # Format the entry name based on if the gene ID is in the database dictionary
                entry_name = ""
                if gene_number in database.keys():
                    entry_name = entry_name + database[gene_number]
                else:
                    entry_name += gene_number

                gene_number_list = []
                for i, gene_number in enumerate(entry_split):
                    gene_number_list.append(gene_number.split()[0])

                for gene_number in gene_number_list:
                    entry_name = (
                        entry_name + "-" + database[gene_number]
                        if gene_number in database.keys()
                        else entry_name + "-" + gene_number
                    )
                entry_type = entry["type"]
            else:
                entry_name = entry["name"]
                entry_type = entry["type"]
        else:
            # hsa:100132463
            if entry_split[0] == "hsa":
                entry_name = entry_split[1]
                entry_type = entry["type"]
                entry_name = (
                    hsaDict[entry_name] if entry_name in hsaDict.keys() else entry_name
                )
            elif entry_split[0] == "ko":
                entry_name = entry_split[1]
                entry_type = entry["type"]
                entry_name = (
                    KEGGdict[entry_name]
                    if entry_name in KEGGdict.keys()
                    else entry_name
                )
            elif entry_split[0] == "path":
                entry_name = entry["name"]
                entry_type = "path"
            else:
                entry_name = entry["name"]
                entry_type = entry["type"]

        entry_id = entry["id"]
        entry_name = re.sub(",", "", entry_name)
        id_to_name[entry_id] = entry_name

        if entry_type == "group":
            group_ids = []
            for component in entry.find_all("component"):
                group_ids.append(component["id"])
            groups[entry_id] = group_ids
        else:
            graph.add_node(entry_name, name=entry_name, type=entry_type)

    for relation in soup.find_all("relation"):
        (color, signal) = ("black", "a")

        relation_entry1 = relation["entry1"]
        relation_entry2 = relation["entry2"]
        relation_type = relation["type"]

        subtypes = []

        for subtype in relation.find_all("subtype"):
            subtypes.append(subtype["name"])

        if (
                ("activation" in subtypes)
                or ("expression" in subtypes)
                or ("glycosylation" in subtypes)
        ):
            color = "green"
            signal = "a"
        elif ("inhibition" in subtypes) or ("repression" in subtypes):
            color = "red"
            signal = "i"
        elif ("binding/association" in subtypes) or ("compound" in subtypes):
            color = "purple"
            signal = "a"
        elif "phosphorylation" in subtypes:
            color = "orange"
            signal = "a"
        elif "dephosphorylation" in subtypes:
            color = "pink"
            signal = "i"
        elif "indirect effect" in subtypes:
            color = "cyan"
            signal = "a"
        elif "dissociation" in subtypes:
            color = "yellow"
            signal = "i"
        elif "ubiquitination" in subtypes:
            color = "cyan"
            signal = "i"
        else:
            print("color not detected. Signal assigned to activation arbitrarily")
            print(subtypes)
            signal = "a"

        entry1_list = deconvolute_groups(relation_entry1, groups)
        entry2_list = deconvolute_groups(relation_entry2, groups)

        for (entry1, entry2) in itertools.product(entry1_list, entry2_list):
            node1 = id_to_name[entry1]
            node2 = id_to_name[entry2]
            graph.add_edge(
                node1,
                node2,
                color=color,
                subtype="/".join(subtypes),
                type=relation_type,
                signal=signal,
            )

    return graph


def find_kegg_pathways(
    gene_list, preDefList=[], write_graphml=True, organism="hsa", minimumOverlap=1
):

    """
    geneList = the list of genes included in dataset
    write_graphml = whether or not to write out a graphml (usually true)
    organism = organism code from kegg. Eg human = 'hsa', mouse = 'mus'
    """
    koDict = parse_kegg_dict()  # parse the dictionary of ko codes
    try:  # try to retrieve and parse the dictionary containing organism gene names to codes conversion
        url = requests.get("http://rest.kegg.jp/list/" + organism, stream=True)
        # reads KEGG dictionary of identifiers between numbers and actual protein names and saves it to a python dictionary
        aliasDict = {}
        orgDict = {}
        for line in url.iter_lines():
            line = line.decode("utf-8")
            line_split = line.split("\t")
            k = line_split[0].split(":")[1]
            nameline = line_split[1].split(";")
            name = nameline[0]
            if "," in name:
                nameline = name.split(",")
                name = nameline[0]
                for entry in range(1, len(nameline)):
                    aliasDict[nameline[entry].strip()] = name.upper()
            orgDict[k] = name
    except:
        print("Could not get library: " + organism)
    k = KEGG()  # read KEGG from bioservices
    k.organism = organism
    if len(preDefList) == 0:
        pathwayList = list(k.pathwayIds)
    else:
        pathwayList = list(preDefList)
    print(pathwayList)
    pathwayDict = {}
    for x in pathwayList:
        x = x.replace("path:", "")
        code = str(x)
        code = re.sub(
            "[a-zA-Z]+", "", code
        )  # eliminate org letters - retain only numbers from KEGG pathway codes
        origCode = code
        coder = str("ko" + code)  # add ko
        graph = nx.DiGraph()  # open a graph object
        # get ko pathway
        for code in [coder]:
            print(code)
            try:
                url = requests.get(
                    "http://rest.kegg.jp/get/" + code + "/kgml", stream=True
                )
                text = [line.decode("utf-8") for line in url.iter_lines()]
            except:
                print("could not read code: " + code)
                continue
            # read kegg
            graph = readKEGG(text, graph, koDict, orgDict)
        coder = str(organism + origCode)  # set up with org letters
        # get org pathway
        text = []
        for code in [coder]:
            try:
                url = requests.get(
                    "http://rest.kegg.jp/get/" + code + "/kgml", stream=True
                )
                text = []
                for line in url.iter_lines():
                    line = line.decode("utf-8")
                    text.append(line)
            except:
                print("could not read code: " + code)
                continue
            # read kegg
            graph = readKEGG(text, graph, koDict, orgDict)

        # remove complexes and rewire components
        removeNodeList = [x for x in graph.nodes() if "-" in x]
        for rm in removeNodeList:
            for start in graph.predecessors(rm):
                edge1 = graph.get_edge_data(start, rm)["signal"]
                if edge1 == "i":
                    for element in rm.split("-"):
                        graph.add_edge(start, element, signal="i")
                else:
                    for element in rm.split("-"):
                        graph.add_edge(start, element, signal="a")
            for finish in graph.successors(rm):
                edge2 = graph.get_edge_data(rm, finish)["signal"]
                if edge2 == "i":
                    for element in rm.split("-"):
                        graph.add_edge(element, finish, signal="i")
                else:
                    for element in rm.split("-"):
                        graph.add_edge(element, finish, signal="a")
            graph.remove_node(rm)
        # remove dependence of nodes on complexes that include that node
        for node in list(graph.nodes()):
            predlist = graph.predecessors(node)
            for pred in predlist:
                if "-" in pred:
                    genes = pred.split("-")
                    flag = True
                    for gene in genes:
                        if not gene in predlist:
                            flag = False
                    if flag:
                        graph.remove_edge(pred, node)
        # remove self edges
        for edge in list(graph.edges()):
            if edge[0] == edge[1]:
                graph.remove_edge(edge[0], edge[1])
        # check to see if there is a connected component, simplify graph and print if so
        allNodes = set(graph.nodes())
        test = len(allNodes.intersection(gene_list))
        print("Pathway: ", x, " Overlap: ", test, " Edges: ", len(graph.edges()))
        # nx.write_graphml(graph,coder+'.graphml')
        if (
            test > minimumOverlap and len(graph.edges()) > 0
        ):  # if there are at least minimumOverlap genes shared between the network and the genes in the dataset
            # nx.write_graphml(graph,coder+'_processed.graphml') # write graph out
            nx.write_graphml(graph, coder + ".graphml")
            print(
                "nodes: ",
                str(len(graph.nodes())),
                ",   edges:",
                str(len(graph.edges())),
            )
            # print(graph.nodes())
            pathwayDict[code] = graph
    return pathwayDict
