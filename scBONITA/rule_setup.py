import numpy as np
from scipy.stats.stats import spearmanr

def parse_node_activity(graph, nodeList, possibilities, i):
    # create a list of the activities of each node and store alongside the contributors to each and node for easy reference later
    activities = []  # list to store activities of nodes (a vs i)
    activity = []
    for sequence in possibilities:
        activity = []
        for node in sequence:
            print(node)
            print(graph[nodeList[node]][nodeList[i]].keys())
            # check the 'interaction' edge attribute
            if "interaction" in list(graph[nodeList[node]][nodeList[i]].keys()):
                if graph[nodeList[node]][nodeList[i]]["interaction"] == "a":
                    activity.append(False)
                else:
                    if graph[nodeList[node]][nodeList[i]]["interaction"] == "i":
                        activity.append(True)
                    else:
                        if (
                                graph[nodeList[node]][nodeList[i]]["interaction"]
                                == "u"
                        ):
                            print(
                                "Unknown interaction type " + graph[nodeList[node]][nodeList[i]][
                                    "interaction"] + " assigning activation..."
                            )
                            activity.append(False)
                        else:
                            if (
                                    graph[nodeList[node]][nodeList[i]][
                                        "interaction"
                                    ]
                                    == "g"
                            ):
                                print(
                                    "Group edge/interaction type, assigning activation..."
                                )
                                activity.append(False)
                            else:
                                print(
                                    "Unknown interaction type " + graph[nodeList[node]][nodeList[i]][
                                        "interaction"] + " assigning activation..."
                                )
                                activity.append(False)
            # check the 'signal' edge attribute
            if "signal" in list(graph[nodeList[node]][nodeList[i]].keys()):
                if graph[nodeList[node]][nodeList[i]]["signal"] == "a":
                    activity.append(False)
                else:
                    if graph[nodeList[node]][nodeList[i]]["signal"] == "i":
                        activity.append(True)
                    else:
                        if graph[nodeList[node]][nodeList[i]]["signal"] == "u":
                            print(
                                "Unknown interaction type " + graph[nodeList[node]][nodeList[i]][
                                    "signal"] + " assigning activation..."
                            )
                            activity.append(False)
                        else:
                            if (
                                    graph[nodeList[node]][nodeList[i]]["signal"]
                                    == "g"
                            ):
                                print(
                                    "Group edge/interaction type, assigning activation..."
                                )
                                activity.append(False)
                            else:
                                print(
                                    "Unknown interaction type " + graph[nodeList[node]][nodeList[i]][
                                        "signal"] + " assigning activation..."
                                )
                                activity.append(False)
            # If neither edge attribute is present, assign activation
            if not "interaction" in list(
                    graph[nodeList[node]][nodeList[i]].keys()
            ) and not "signal" in list(
                graph[nodeList[node]][nodeList[i]].keys()
            ):
                print("Group edge/interaction type, assigning activation...")
                activity.append(False)

        activities.append(activity)
    return activities, activity

def parse_node_activity2(graph, nodeList, possibilities, i):
    # create a list of the activities of each node and store alongside the contributors to each and node for easy reference later
    activities = []  # list to store activities of nodes (a vs i)
    activity = []
    for sequence in possibilities:
        activity = []
        for node in sequence:
            # check the 'interaction' edge attribute
            if "interaction" in list(graph[nodeList[node]][nodeList[i]].keys()):
                if graph[nodeList[node]][nodeList[i]]["interaction"] == "a":
                    activity.append(False)
                else:
                    if graph[nodeList[node]][nodeList[i]]["interaction"] == "i":
                        activity.append(True)
                    else:
                        if (
                                graph[nodeList[node]][nodeList[i]]["interaction"]
                                == "u"
                        ):
                            print(
                                "Unknown interaction type " + graph[nodeList[node]][nodeList[i]][
                                    "interaction"] + " assigning activation..."
                            )
                            activity.append(False)
                        else:
                            if (
                                    graph[nodeList[node]][nodeList[i]][
                                        "interaction"
                                    ]
                                    == "g"
                            ):
                                print(
                                    "Group edge/interaction type, assigning activation..."
                                )
                                activity.append(False)
                            else:
                                print(
                                    "Unknown interaction type " + graph[nodeList[node]][nodeList[i]][
                                        "interaction"] + " assigning activation..."
                                )
                                activity.append(False)
            # check the 'signal' edge attribute
            if "signal" in list(graph[nodeList[node]][nodeList[i]].keys()):
                if graph[nodeList[node]][nodeList[i]]["signal"] == "a":
                    activity.append(False)
                else:
                    if graph[nodeList[node]][nodeList[i]]["signal"] == "i":
                        activity.append(True)
                    else:
                        if graph[nodeList[node]][nodeList[i]]["signal"] == "u":
                            print(
                                "Unknown interaction type " + graph[nodeList[node]][nodeList[i]][
                                    "signal"] + " assigning activation..."
                            )
                            activity.append(False)
                        else:
                            if (
                                    graph[nodeList[node]][nodeList[i]]["signal"]
                                    == "g"
                            ):
                                print(
                                    "Group edge/interaction type, assigning activation..."
                                )
                                activity.append(False)
                            else:
                                print(
                                    "Unknown interaction type " + graph[nodeList[node]][nodeList[i]][
                                        "signal"] + " assigning activation..."
                                )
                                activity.append(False)
            # If neither edge attribute is present, assign activation
            if not "interaction" in list(
                    graph[nodeList[node]][nodeList[i]].keys()
            ) and not "signal" in list(
                graph[nodeList[node]][nodeList[i]].keys()
            ):
                print("Group edge/interaction type, assigning activation...")
                activity.append(False)

        activities.append(activity)
    return activities, activity



