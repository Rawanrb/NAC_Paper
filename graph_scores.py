from networkx.algorithms import tree
import networkx as nx
import numpy as np
from graph_analysis import graph_creating as gc
from graph_analysis import graph_processing as gp
import glob
import os
import csv
from graph_analysis import graph_communities as gcom


def minimum_spanning_edges_function(G, algorithm="kruskal"):
    # "prim"
    mst = tree.minimum_spanning_edges(G, algorithm=algorithm, data=False)
    edgelist = list(mst)
    return sorted(sorted(e) for e in edgelist)


def calculate_props_whole_graph(G):
    """calculate the props for the graph of specific class"""
    graph_data = []
    edges_mst = minimum_spanning_edges_function(G)
    print('# of edges of MST: {}'.format(len(edges_mst)))

    graph_data.append(len(edges_mst))
    # print('# of edges: {}'.format(G.number_of_edges()))
    # print('# of nodes: {}'.format(G.number_of_nodes()))
    # print(nx.algorithms.components.number_connected_components(G))

    # print('Number of nodes with connectivity:{}'.format(nx.node_connectivity(G)))
    # graph_data.append(nx.node_connectivity(G))

    # dic = nx.algorithms.centrality.degree_centrality(G)
    #
    # c = [dic[v] for v in dic if v > 0]
    #
    # if len(c) > 0:
    #     print('The average of degree_centrality: {}'.format(sum(c) / len(c)))
    #     graph_data.append(sum(c) / len(c))
    # else:
    #     graph_data.append(0)

    # dic = nx.algorithms.centrality.betweenness_centrality(G)
    #
    # c = [dic[v] for v in dic if v > 0]
    #
    # if len(c) > 0:
    #     print('The average of betweenness_centrality: {}'.format(sum(c) / len(c)))
    #     graph_data.append(sum(c) / len(c))
    # else:
    #     graph_data.append(0)

    # edge_betweenness_centrality
    # dic = nx.algorithms.centrality.edge_betweenness_centrality(G)
    #
    # c = [v for v in dic.values() if v > 0]
    #
    # if len(c) > 0:
    #     print('The average of edge_betweenness_centrality: {}'.format(sum(c) / len(c)))
    #     graph_data.append(sum(c) / len(c))
    # else:
    #     graph_data.append(0)

    # dic = nx.algorithms.centrality.closeness_centrality(G)
    #
    # c = [dic[v] for v in dic if v > 0]
    #
    # if len(c) > 0:
    #     print('The average of closeness_centrality: {}'.format(sum(c) / len(c)))
    #     graph_data.append(sum(c) / len(c))
    # else:
    #     graph_data.append(0)

    # Select a list of influential nodes in a graph using VoteRank algorithm
    # print('voterank: {}'.format(
    #     len(nx.algorithms.centrality.voterank(G))))
    # if len(G.nodes()) > 0:
    #
    #     graph_data.append(len(nx.algorithms.centrality.voterank(G)))
    #
    #
    #
    # # print('Average clustering coefficient for nodes: {}'.format(
    # #     nx.algorithms.cluster.average_clustering(G)))
    #     graph_data.append(nx.algorithms.cluster.average_clustering(G))
    # else:
    #     graph_data.append(0)
    #     graph_data.append(0)

    return graph_data


def calculate_props_whole_graph_connected_components(G, classes=[], advance=False):
    props_data = []
    props_data_advance = []
    tissue_frag = 0
    tissue_sheet = 0
    tissue_TB = 0
    subgraphs_with_other_classes = 0
    subgraphs_with_other_classes_perc_per_subgraph = 0
    subgraphs_with_other_classes_perc_per_graph = 0
    other_nodes_total = 0
    destory_power = 0
    MST_average = 0
    MST_average_Whole = 0
    MST_weights_sum_avrg = 0
    MST_weights_sum_perc = 0
    edges_mst_num = 0
    edges_mst_num_local = 0
    edges_mst_num_whole = 0
    node_connectivity_sum = 0
    average_clustering_coefficient = 0
    average_cut_value = 0
    average_shortest_path_length_sum = 0
    shared_edges_perc = 0
    other_nodes_perc = 0
    wiener_index_aver_sum = 0
    voterank_list = []

    print('# of connected_components: {}'.format(nx.algorithms.components.number_connected_components(G)))

    props_data.append(nx.algorithms.components.number_connected_components(G))
    props_data.append(G.number_of_nodes())
    props_data.append(G.number_of_edges())

    eccentricity_averg = 0
    diameter_averg = 0
    radius_averg = 0
    degree_centrality_averg = 0
    betweenness_centrality_averg = 0
    edge_betweenness_centrality = 0
    closeness_centrality_averg = 0
    eccentricity = 0
    diameter = 0
    radius = 0
    wiener_index_sum = 0

    if advance:
        eccentricity = nx.algorithms.distance_measures.eccentricity(G)
        props_data_advance.append(eccentricity)
        print('Eccentricity: {}'.format(
            eccentricity))
        diameter = nx.algorithms.distance_measures.diameter(G)
        props_data_advance.append(diameter)
        print('Diameter: {}'.format(
            diameter))
        radius = nx.algorithms.distance_measures.radius(G)
        props_data_advance.append(radius)
        print('Radius: {}'.format(
            radius))

    for component in nx.connected_components(G):
        comp = G.subgraph(list(component))
        if len(component) < 3:
            tissue_TB += 1
        elif len(component) < 300:
            tissue_frag += 1
        else:
            tissue_sheet += 1

        edges_mst = minimum_spanning_edges_function(comp)

        if len(edges_mst) > 0:
            print('# of edges of MST: {}'.format(len(edges_mst)))
            edges_mst_num += len(edges_mst)
            edges_mst_num_local += len(edges_mst) / comp.number_of_edges()
            edges_mst_num_whole += len(edges_mst) / G.number_of_edges()
            print('Percetage of edges of MST of comps edges: {}'.format((len(edges_mst) / comp.number_of_edges())))
            MST_average += (len(edges_mst) / comp.number_of_edges())
            MST_average_Whole += (
                        (len(edges_mst) / comp.number_of_edges()) * (comp.number_of_edges() / G.number_of_edges()))
            weights_sum = 0
            for edge in edges_mst:
                weights_sum += G.edges[edge]['weight']

            weights_sum_comp = 0
            for edge in comp.edges():
                weights_sum_comp += G.edges[edge]['weight']

            print('Average weights sum of MST edges: {}'.format(weights_sum / len(edges_mst)))

            MST_weights_sum_avrg += weights_sum / len(edges_mst)
            MST_weights_sum_perc += weights_sum / weights_sum_comp

            W = nx.wiener_index(comp, weight='weight')
            print('Wiener index: {}'.format(W))
            wiener_index_sum += W
            wiener_index_aver_sum += (W / weights_sum_comp)

        print('# of edges: {}'.format(comp.number_of_edges()))
        print('# of nodes: {}'.format(comp.number_of_nodes()))

        # print('Number of nodes with connectivity:{}'.format(nx.node_connectivity(comp)))
        # node_connectivity_sum += nx.node_connectivity(comp)

        if advance:
            dic = nx.algorithms.centrality.degree_centrality(comp)
            c = [dic[v] for v in dic if v > 0]
            if len(c) > 0:
                print('The average of degree_centrality: {}'.format(sum(c) / len(c)))
                degree_centrality_averg += (sum(c) / len(c)) * (comp.number_of_nodes() / G.number_of_nodes())

            dic = nx.algorithms.centrality.betweenness_centrality(comp)

            c = [dic[v] for v in dic if v > 0]
            if len(c) > 0:
                print('The average of betweenness_centrality: {}'.format(sum(c) / len(c)))
                betweenness_centrality_averg += (sum(c) / len(c)) * (comp.number_of_nodes() / G.number_of_nodes())
            # edge_betweenness_centrality
            dic = nx.algorithms.centrality.edge_betweenness_centrality(comp)
            c = [v for v in dic.values() if v > 0]
            if len(c) > 0:
                print('The average of edge_betweenness_centrality: {}'.format(sum(c) / len(c)))
                edge_betweenness_centrality += (sum(c) / len(c)) * (comp.number_of_edges() / G.number_of_edges())

            dic = nx.algorithms.centrality.closeness_centrality(comp)

            c = [dic[v] for v in dic if v > 0]
            if len(c) > 0:
                print('The average of closeness_centrality: {}'.format(sum(c) / len(c)))
                closeness_centrality_averg += (sum(c) / len(c)) * (comp.number_of_nodes() / G.number_of_nodes())

            # print('all shortest path length= {}'.format(Pairs[1]))
            # print('SD of all shortest path length= {}'.format(np.std(Pairs[1])))
            # print('Mean of all shortest path length= {}'.format(np.mean(Pairs[1])))

            e = nx.algorithms.distance_measures.eccentricity(comp)
            print('Eccentricity: {}'.format(
                e))
            eccentricity_averg += (e / eccentricity)

            d = nx.algorithms.distance_measures.diameter(comp)
            print('Diameter: {}'.format(
                d))
            diameter_averg += (d / diameter)
            print('Periphery: {}'.format(
                nx.algorithms.distance_measures.periphery(comp)))
            r = nx.algorithms.distance_measures.radius(comp)
            print('Radius: {}'.format(
                r))
            radius_averg += (r / radius)

        # Select a list of influential nodes in a graph using VoteRank algorithm
        # voterank_list.append(nx.algorithms.centrality.voterank(comp))
        # print('voterank list: {}'.format(
        #     nx.algorithms.centrality.voterank(comp)))

        #
        # print('Average clustering coefficient for nodes: {}'.format(
        #     nx.algorithms.cluster.average_clustering(comp)))
        #
        # average_clustering_coefficient += nx.algorithms.cluster.average_clustering(comp)
        # # plot  clustering coefficient and average path length  using watts strogatz or Barabási–Albert
        # the prop of lym attacking tumour

        # The sum of weights of edges in a minimum cut
        # if len(component) > 3:
        #     cut_value, partition = nx.stoer_wagner(comp)
        #     print('cut_value: {}'.format(cut_value))
        #     average_cut_value += (cut_value/weights_sum_comp)

        # number of edges whose deletion increases the graph's number
        # for com in nx.algorithms.connectivity.edge_kcomponents.bridge_components(comp):
        #   print(com)
        # print(nx.algorithms.connectivity.edge_kcomponents.bridge_components(comp))

        # print('average_shortest_path_length: {}'.format(nx.average_shortest_path_length(comp)))
        # average_shortest_path_length_sum += nx.average_shortest_path_length(comp)

        # if I have more than one class in the sent G
        # more props are analysed

        shared_edges = 0
        other_nodes = 0
        if len(classes) > 1:

            for u, v, a in comp.edges(data=True):
                if comp.edges[u, v]['label'] != classes[0] + '_' + classes[0]:
                    shared_edges += 1

                if comp.nodes[u]['class_name'] != classes[0]:
                    other_nodes += 1
                    other_nodes_total += 1

                if comp.nodes[v]['class_name'] != classes[0]:
                    other_nodes += 1
                    other_nodes_total += 1

            if other_nodes > 0:
                if comp.number_of_edges() > 0:
                    shared_edges_perc += shared_edges / comp.number_of_edges()
                other_nodes_perc += other_nodes / comp.number_of_nodes()
                subgraphs_with_other_classes += 1
                subgraphs_with_other_classes_perc_per_subgraph += (other_nodes / comp.number_of_nodes())

                subgraphs_with_other_classes_perc = (other_nodes / comp.number_of_nodes())

                subgraphs_with_other_classes_perc_per_graph += subgraphs_with_other_classes_perc * (
                            comp.number_of_nodes() / G.number_of_nodes())
                if comp.number_of_edges() > 0:
                    destory_power += (shared_edges / (comp.number_of_edges() * (comp.number_of_edges() + 1)))
                    print('Percentage shared edges: {}'.format((shared_edges / comp.number_of_edges())))

            print('Percentage of nodes from other classes: {}'.format((other_nodes / comp.number_of_nodes())))

    props_data.append(tissue_frag)
    props_data.append(tissue_sheet)
    props_data.append(tissue_TB)
    props_data.append(edges_mst_num_local)
    props_data.append(edges_mst_num_whole)
    props_data.append(MST_average)
    props_data.append(MST_average_Whole)
    props_data.append(MST_weights_sum_avrg)
    props_data.append(MST_weights_sum_perc)
    props_data.append(node_connectivity_sum)
    props_data.append(average_clustering_coefficient)
    props_data.append(average_cut_value)
    props_data.append(average_shortest_path_length_sum)
    props_data.append(shared_edges_perc)
    props_data.append(other_nodes_perc)
    props_data.append(wiener_index_sum)
    props_data.append(wiener_index_aver_sum)

    voterank_nodes_number = 0
    voterank_nodes_number_other = 0

    # for item in voterank_list:
    #     if len(item)>0:
    #         for i in item:
    #             voterank_nodes_number += 1
    #             if G.nodes[i]['class_name'] != classes[0]:
    #                 voterank_nodes_number_other +=1

    voterank_nodes_number_other = 1

    props_data.append(voterank_nodes_number)
    props_data.append(voterank_nodes_number_other)
    props_data.append(voterank_nodes_number / voterank_nodes_number_other)

    if len(classes) > 1:
        print('Percentage subgraphs with other classes: {}'.format(
            (subgraphs_with_other_classes / nx.algorithms.components.number_connected_components(G))))
        props_data.append((subgraphs_with_other_classes / nx.algorithms.components.number_connected_components(G)))
        print('Percentage subgraphs with other classes: {}'.format(subgraphs_with_other_classes_perc_per_graph))
        props_data.append(subgraphs_with_other_classes_perc_per_graph)
        print('Percentage other classes nodes: {}'.format(other_nodes_total / G.number_of_nodes()))
        props_data.append(other_nodes_total / G.number_of_nodes())
        print('Percentage main class nodes: {}'.format((G.number_of_nodes() - other_nodes_total) / G.number_of_nodes()))
        props_data.append((G.number_of_nodes() - other_nodes_total) / G.number_of_nodes())
        # how much time does tissue needs to destory the current one by exploring MST
        # it is like tree with injections
        print('Destory power: {}'.format(destory_power))
        props_data.append(destory_power)
        props_data.append(subgraphs_with_other_classes_perc_per_subgraph)

    if advance:
        props_data_advance.append(eccentricity_averg)
        props_data_advance.append(diameter_averg)
        props_data_advance.append(radius_averg)
        props_data_advance.append(degree_centrality_averg)
        props_data_advance.append(betweenness_centrality_averg)
        props_data_advance.append(edge_betweenness_centrality)
        props_data_advance.append(closeness_centrality_averg)
    # number of cluster per unit
    # til cluster size
    # svm of the scores
    # correlation

    # plot  clustering coefficient and average path length  using watts strogatz or Barabási–Albert
    # they can be applied to each cluster to find out how can we connect them (small world)
    # chain #-1 ,1 depends on the label # the strenght of the attack

    if advance:
        return props_data, props_data_advance

    return props_data


def calculate_props_of_communities(coms, G, classes=[], advance=False):
    coms_dic = {}
    for key, value in coms.items():
        if value[0] in coms_dic.keys():
            coms_dic[value[0]].append(key)
        else:
            coms_dic[value[0]] = [key]

    props_data = []
    props_data_advance = []
    tissue_frag = 0
    tissue_sheet = 0
    tissue_TB = 0
    subgraphs_with_other_classes = 0
    subgraphs_with_other_classes_perc_per_subgraph = 0
    subgraphs_with_other_classes_perc_per_graph = 0
    other_nodes_total = 0
    destory_power = 0
    MST_average = 0
    MST_average_Whole = 0
    MST_weights_sum_avrg = 0
    MST_weights_sum_perc = 0
    edges_mst_num = 0
    edges_mst_num_local = 0
    edges_mst_num_whole = 0
    node_connectivity_sum = 0
    average_clustering_coefficient = 0
    average_cut_value = 0
    average_shortest_path_length_sum = 0
    shared_edges_perc = 0
    other_nodes_perc = 0
    wiener_index_aver_sum = 0
    voterank_list = []

    print('# of connected_components: {}'.format(nx.algorithms.components.number_connected_components(G)))

    props_data.append(nx.algorithms.components.number_connected_components(G))
    props_data.append(G.number_of_nodes())
    props_data.append(G.number_of_edges())

    eccentricity_averg = 0
    diameter_averg = 0
    radius_averg = 0
    degree_centrality_averg = 0
    betweenness_centrality_averg = 0
    edge_betweenness_centrality = 0
    closeness_centrality_averg = 0
    eccentricity = 0
    diameter = 0
    radius = 0

    if advance:
        eccentricity = nx.algorithms.distance_measures.eccentricity(G)
        props_data_advance.append(eccentricity)
        print('Eccentricity: {}'.format(
            eccentricity))
        diameter = nx.algorithms.distance_measures.diameter(G)
        props_data_advance.append(diameter)
        print('Diameter: {}'.format(
            diameter))
        radius = nx.algorithms.distance_measures.radius(G)
        props_data_advance.append(radius)
        print('Radius: {}'.format(
            radius))

    for key, value in coms_dic.items():
        comp = G.subgraph(value)
        if len(value) < 3:
            tissue_TB += 1
        elif len(value) < 300:
            tissue_frag += 1
        else:
            tissue_sheet += 1

        edges_mst = minimum_spanning_edges_function(comp)

        if len(edges_mst) > 0:
            print('# of edges of MST: {}'.format(len(edges_mst)))
            edges_mst_num += len(edges_mst)
            edges_mst_num_local += len(edges_mst) / comp.number_of_edges()
            edges_mst_num_whole += len(edges_mst) / G.number_of_edges()
            print('Percetage of edges of MST of comps edges: {}'.format((len(edges_mst) / comp.number_of_edges())))
            MST_average += (len(edges_mst) / comp.number_of_edges())
            MST_average_Whole += (
                    (len(edges_mst) / comp.number_of_edges()) * (comp.number_of_edges() / G.number_of_edges()))
            weights_sum = 0
            for edge in edges_mst:
                weights_sum += G.edges[edge]['weight']

            weights_sum_comp = 0
            for edge in comp.edges():
                weights_sum_comp += G.edges[edge]['weight']

            print('Average weights sum of MST edges: {}'.format(weights_sum / len(edges_mst)))

            MST_weights_sum_avrg += weights_sum / len(edges_mst)
            MST_weights_sum_perc += weights_sum / weights_sum_comp

            W = nx.wiener_index(comp, weight='weight')
            print('Wiener index: {}'.format(W))
            wiener_index_aver_sum += (W / weights_sum_comp)

        print('# of edges: {}'.format(comp.number_of_edges()))
        print('# of nodes: {}'.format(comp.number_of_nodes()))

        print('Number of nodes with connectivity:{}'.format(nx.node_connectivity(comp)))
        node_connectivity_sum += nx.node_connectivity(comp)

        if advance:
            dic = nx.algorithms.centrality.degree_centrality(comp)
            c = [dic[v] for v in dic if v > 0]
            if len(c) > 0:
                print('The average of degree_centrality: {}'.format(sum(c) / len(c)))
                degree_centrality_averg += (sum(c) / len(c)) * (comp.number_of_nodes() / G.number_of_nodes())

            dic = nx.algorithms.centrality.betweenness_centrality(comp)

            c = [dic[v] for v in dic if v > 0]
            if len(c) > 0:
                print('The average of betweenness_centrality: {}'.format(sum(c) / len(c)))
                betweenness_centrality_averg += (sum(c) / len(c)) * (comp.number_of_nodes() / G.number_of_nodes())
            # edge_betweenness_centrality
            dic = nx.algorithms.centrality.edge_betweenness_centrality(comp)
            c = [v for v in dic.values() if v > 0]
            if len(c) > 0:
                print('The average of edge_betweenness_centrality: {}'.format(sum(c) / len(c)))
                edge_betweenness_centrality += (sum(c) / len(c)) * (comp.number_of_edges() / G.number_of_edges())

            dic = nx.algorithms.centrality.closeness_centrality(comp)

            c = [dic[v] for v in dic if v > 0]
            if len(c) > 0:
                print('The average of closeness_centrality: {}'.format(sum(c) / len(c)))
                closeness_centrality_averg += (sum(c) / len(c)) * (comp.number_of_nodes() / G.number_of_nodes())

            # print('all shortest path length= {}'.format(Pairs[1]))
            # print('SD of all shortest path length= {}'.format(np.std(Pairs[1])))
            # print('Mean of all shortest path length= {}'.format(np.mean(Pairs[1])))

            e = nx.algorithms.distance_measures.eccentricity(comp)
            print('Eccentricity: {}'.format(
                e))
            eccentricity_averg += (e / eccentricity)

            d = nx.algorithms.distance_measures.diameter(comp)
            print('Diameter: {}'.format(
                d))
            diameter_averg += (d / diameter)
            print('Periphery: {}'.format(
                nx.algorithms.distance_measures.periphery(comp)))
            r = nx.algorithms.distance_measures.radius(comp)
            print('Radius: {}'.format(
                r))
            radius_averg += (r / radius)

        # Select a list of influential nodes in a graph using VoteRank algorithm
        voterank_list.append(nx.algorithms.centrality.voterank(comp))
        # print('voterank list: {}'.format(
        #     nx.algorithms.centrality.voterank(comp)))

        print('Average clustering coefficient for nodes: {}'.format(
            nx.algorithms.cluster.average_clustering(comp)))

        average_clustering_coefficient += nx.algorithms.cluster.average_clustering(comp)
        # plot  clustering coefficient and average path length  using watts strogatz or Barabási–Albert
        # the prop of lym attacking tumour

        # The sum of weights of edges in a minimum cut
        if len(value) > 3:
            cut_value, partition = nx.stoer_wagner(comp)
            print('cut_value: {}'.format(cut_value))
            average_cut_value += (cut_value / weights_sum_comp)

        # number of edges whose deletion increases the graph's number
        # for com in nx.algorithms.connectivity.edge_kcomponents.bridge_components(comp):
        #   print(com)
        # print(nx.algorithms.connectivity.edge_kcomponents.bridge_components(comp))

        print('average_shortest_path_length: {}'.format(nx.average_shortest_path_length(comp)))
        average_shortest_path_length_sum += nx.average_shortest_path_length(comp)

        # if I have more than one class in the sent G
        # more props are analysed

        shared_edges = 0
        other_nodes = 0
        if len(classes) > 1:

            for u, v, a in comp.edges(data=True):
                if comp.edges[u, v]['label'] != classes[0] + '_' + classes[0]:
                    shared_edges += 1

                if comp.nodes[u]['class_name'] != classes[0]:
                    other_nodes += 1
                    other_nodes_total += 1

                if comp.nodes[v]['class_name'] != classes[0]:
                    other_nodes += 1
                    other_nodes_total += 1

            if other_nodes > 0:
                if comp.number_of_edges() > 0:
                    shared_edges_perc += shared_edges / comp.number_of_edges()
                other_nodes_perc += other_nodes / comp.number_of_nodes()
                subgraphs_with_other_classes += 1
                subgraphs_with_other_classes_perc_per_subgraph += (other_nodes / comp.number_of_nodes())

                subgraphs_with_other_classes_perc = (other_nodes / comp.number_of_nodes())

                subgraphs_with_other_classes_perc_per_graph += subgraphs_with_other_classes_perc * (
                        comp.number_of_nodes() / G.number_of_nodes())
                if comp.number_of_edges() > 0:
                    destory_power += (shared_edges / (comp.number_of_edges() * (comp.number_of_edges() + 1)))
                    print('Percentage shared edges: {}'.format((shared_edges / comp.number_of_edges())))

            print('Percentage of nodes from other classes: {}'.format((other_nodes / comp.number_of_nodes())))

    props_data.append(tissue_frag)
    props_data.append(tissue_sheet)
    props_data.append(tissue_TB)
    props_data.append(edges_mst_num_local)
    props_data.append(edges_mst_num_whole)
    props_data.append(MST_average)
    props_data.append(MST_average_Whole)
    props_data.append(MST_weights_sum_avrg)
    props_data.append(MST_weights_sum_perc)
    props_data.append(node_connectivity_sum)
    props_data.append(average_clustering_coefficient)
    props_data.append(average_cut_value)
    props_data.append(average_shortest_path_length_sum)
    props_data.append(shared_edges_perc)
    props_data.append(other_nodes_perc)
    props_data.append(wiener_index_aver_sum)

    voterank_nodes_number = 0
    voterank_nodes_number_other = 0

    for item in voterank_list:
        if len(item) > 0:
            for i in item:
                voterank_nodes_number += 1
                if G.nodes[i]['class_name'] != classes[0]:
                    voterank_nodes_number_other += 1

    props_data.append(voterank_nodes_number)
    props_data.append(voterank_nodes_number_other)
    props_data.append(voterank_nodes_number / voterank_nodes_number_other)

    if len(classes) > 1:
        print('Percentage subgraphs with other classes: {}'.format(
            (subgraphs_with_other_classes / nx.algorithms.components.number_connected_components(G))))
        props_data.append((subgraphs_with_other_classes / nx.algorithms.components.number_connected_components(G)))
        print('Percentage subgraphs with other classes: {}'.format(subgraphs_with_other_classes_perc_per_graph))
        props_data.append(subgraphs_with_other_classes_perc_per_graph)
        print('Percentage other classes nodes: {}'.format(other_nodes_total / G.number_of_nodes()))
        props_data.append(other_nodes_total / G.number_of_nodes())
        print('Percentage main class nodes: {}'.format(
            (G.number_of_nodes() - other_nodes_total) / G.number_of_nodes()))
        props_data.append((G.number_of_nodes() - other_nodes_total) / G.number_of_nodes())
        # how much time does tissue needs to destory the current one by exploring MST
        # it is like tree with injections
        print('Destory power: {}'.format(destory_power))
        props_data.append(destory_power)
        props_data.append(subgraphs_with_other_classes_perc_per_subgraph)

    if advance:
        props_data_advance.append(eccentricity_averg)
        props_data_advance.append(diameter_averg)
        props_data_advance.append(radius_averg)
        props_data_advance.append(degree_centrality_averg)
        props_data_advance.append(betweenness_centrality_averg)
        props_data_advance.append(edge_betweenness_centrality)
        props_data_advance.append(closeness_centrality_averg)
    # number of cluster per unit
    # til cluster size
    # svm of the scores
    # correlation

    # plot  clustering coefficient and average path length  using watts strogatz or Barabási–Albert
    # they can be applied to each cluster to find out how can we connect them (small world)
    # chain #-1 ,1 depends on the label # the strenght of the attack

    if advance:
        return props_data, props_data_advance

    return props_data

    # plot  clustering coefficient and average path length  using watts strogatz or Barabási–Albert
    # they can be applied to each cluster to find out how can we connect them (small world)
    # chain #-1 ,1 depends on the label # the strenght of the attack
    # how much time does tissue needs to destory the current one by exploring MST
    # it is link tree with injections

    # Pairs = nx.all_pairs_shortest_path_length(G)
    # for pair in Pairs:
    #     print(pair)
    # print('all shortest path length= {}'.format(Pairs[1]))
    # print('SD of all shortest path length= {}'.format(np.std(Pairs[1])))
    # print('Mean of all shortest path length= {}'.format(np.mean(Pairs[1])))

    # #A graph is commonly classified as small-world if sigma>1.
    # print('Sigma: {}'.format(
    #     nx.algorithms.smallworld.sigma(G.subgraph(value))))

    '''The small-world coefficient (omega) ranges between -1 and 1.
    Values close to 0 means the G features small-world characteristics.
    Values close to -1 means G has a lattice shape whereas values close to 1 means G is a random graph.'''
    #
    # print('Omega: {}'.format(
    #     nx.algorithms.smallworld.omega(G.subgraph(value))))

    # print('S_metric: {}'.format(
    #     nx.algorithms.smetric.s_metric(G.subgraph(value))))


def calculate_props_all(param, classes=[], file_name='', advance=False, pickle_file=False):
    with open(param.graph_info + file_name + '.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        header = []
        header.append("ID")
        header.append("number_connected_components")
        header.append("number_of_nodes")
        header.append("number_of_edges")

        header.append("tissue_frag")
        header.append("tissue_sheet")
        header.append("tissue_TB")
        header.append("edges_mst_num_local")
        header.append("edges_mst_num_whole")
        header.append("MST_average")
        header.append("MST_average_Whole")
        header.append("MST_weights_sum_avrg")
        header.append("MST_weights_sum_perc")
        header.append("node_connectivity_sum")
        header.append("average_clustering_coefficient")
        header.append("average_cut_value")
        header.append("average_shortest_path_length_sum")
        header.append("shared_edges_perc")
        header.append("other_nodes_perc")
        header.append("wiener_index_sum")
        header.append("wiener_index_aver_sum")

        header.append("voterank_nodes_number")
        header.append("voterank_nodes_number_other")
        header.append("voterank_nodes_number_perc")

        if len(classes) > 1:
            header.append("subgraphs with others")
            header.append("perc subgraphs with others")
            header.append("perc others nodes")
            header.append("perc main class nodes")
            header.append("destory power")
            header.append("subgraphs_with_other_classes_perc_per_subgraph")

        if advance:
            header.append("eccentricity")
            header.append("diameter ")
            header.append("radius ")
            header.append("eccentricity_averg")
            header.append("diameter_averg")
            header.append("radius_averg")
            header.append("degree_centrality_averg")
            header.append("betweenness_centrality_averg")
            header.append("edge_betweenness_centrality")
            header.append("closeness_centrality_averg")

        writer.writerow(header)

    for wsi_name_ in sorted(glob.glob(param.wsi_path + "/*" + param.wsi_ext)):

        # try:
        with open(param.graph_info + file_name + '.csv', 'a+', newline='') as csvfile:
            writer = csv.writer(csvfile)
            base = os.path.basename(wsi_name_)
            wsi_name = os.path.splitext(base)[0]
            row = []
            row.append(wsi_name)
            if pickle_file:
                X, S, C, W = gc.read_graph_from_np(param, wsi_name)
                G = nx.from_numpy_matrix(W)
            else:
                N = gc.read_graph_from_pickle(param, wsi_name)
                H = gp.extract_subgraphs_based_on_classes(N, classes)
                C = H.copy()

                for u in H.nodes():
                    if H.degree[u] > 2:
                        pass
                    else:
                        C.remove_node(u)
                G = C.copy()
                for u, v, a in C.edges(data=True):
                    if G.nodes[u]['class_name'] == G.nodes[v]['class_name']:
                        G.remove_edge(u, v)

            if advance:
                props_data, props_data_advance = calculate_props_whole_graph_connected_components(G, classes,
                                                                                                  advance=advance)
                for item in props_data:
                    row.append(item)

                for item in props_data_advance:
                    row.append(item)
            else:
                props_data = calculate_props_whole_graph_connected_components(G, classes, advance=advance)
                for item in props_data:
                    row.append(item)

            writer.writerow(row)

        # except:
        #     print('Problem with the slide:%s!' %(wsi_name))
        #     pass


def Shannon_Wiener_Index_two_classes(G, classes, node_size_flag=False):
    print('Start Shannon Wiener Index calculating props!')
    '''species diversity
    H ln of number of each class with their size and prediction in each components
    Hmax = ln of number of species (components)
    '''

    H = 0
    related_comps = 0  # to calculate only components which contain class 1
    if node_size_flag:
        for component in nx.connected_components(G):
            comp = G.subgraph(list(component))

            nodes_types = list(nx.get_node_attributes(comp, 'class_name').values())
            nodes_sizes = list(nx.get_node_attributes(comp, 'node_size').values())
            nodes_predictions = list(nx.get_node_attributes(comp, 'predictions').values())

            class_count_one = sum(map(lambda x: x == classes[1], nodes_types))

            if class_count_one > 0:
                related_comps += 1
                for item in classes:
                    class_count_G = sum(map(lambda x: x == item, nodes_types))
                    if class_count_G > 0:
                        class_indeces = list(map(lambda x: x == item, nodes_types))
                        print(class_indeces)
                        class_sizes_sum = sum(
                            [value for i, value in zip(range(0, len(class_indeces)), nodes_sizes) if class_indeces[i]])
                        class_pred_sum = sum(
                            [value[[item]] for i, value in zip(range(0, len(class_indeces)), nodes_predictions) if
                             class_indeces[i]])

                        H += (class_count_G / class_sizes_sum) * (np.log(class_count_G / class_sizes_sum)) * (
                                    class_pred_sum / class_count_G)
    else:
        for component in nx.connected_components(G):
            comp = G.subgraph(list(component))

            nodes_types = list(nx.get_node_attributes(comp, 'class_name').values())

            nodes_predictions = list(nx.get_node_attributes(comp, 'predictions').values())

            class_count_one = sum(map(lambda x: x == classes[1], nodes_types))

            if class_count_one > 0:
                related_comps += 1
                for item in classes:
                    class_count_G = sum(map(lambda x: x == item, nodes_types))
                    if class_count_G > 0:
                        class_indeces = list(map(lambda x: x == item, nodes_types))

                        class_pred_sum = sum(
                            [value[item] for i, value in zip(range(0, len(class_indeces)), nodes_predictions)
                             if class_indeces[i]])

                        H += (class_count_G / len(nodes_types)) * (np.log(class_count_G / len(nodes_types))) * (
                                class_pred_sum / class_count_G)

    Hmax = np.log(related_comps)

    if Hmax != 0:
        E = np.abs(H) / Hmax
    else:
        E = 0

    return np.abs(H), np.abs(E)


def all_calculate_shanon_wiener_index(param, file_name, classes=[], indeces=[]):
    with open(param.graph_info + file_name + '.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        header = []
        header.append("ID")
        for distance in [100, 500, 1500, 2000, 2500, 3000]:
            header.append("shanon_wiener_index_" + str(distance))
            header.append("evenness_" + str(distance))
        writer.writerow(header)
    for wsi_name_ in sorted(glob.glob(param.wsi_path + "/*" + param.wsi_ext)):
        with open(param.graph_info + file_name + '.csv', 'a+', newline='') as csvfile:
            writer = csv.writer(csvfile)
            base = os.path.basename(wsi_name_)
            wsi_name = os.path.splitext(base)[0]
            print(wsi_name)
            row = []
            row.append(wsi_name)
            N = gc.read_graph_from_pickle(param, wsi_name)
            H = gp.extract_subgraphs_based_on_classes(N, classes)
            for distance in [100, 500, 1500, 2000, 2500, 3000]:
                print(distance)
                G = H.copy()
                for u, v, a in H.edges(data=True):
                    if a['weight'] > distance:
                        if (u, v) in G.edges():
                            G.remove_edge(u, v)

                # props_data = Shannon_Wiener_Index(G, classes,indeces)

                # for item in props_data:
                #     row.append(item)
                print(row)
            writer.writerow(row)


def wiener_index_two_classes(G, classes):
    print('Start wiener_index calculating props!')
    wiener_index_sum_nodes = 0
    wiener_index_sum_edges = 0
    props_data = []

    for component in nx.connected_components(G):
        if len(list(component)) > 1:
            comp = G.subgraph(list(component))

            nodes_types = list(nx.get_node_attributes(comp, 'class_name').values())

            class_one_count = sum(map(lambda x: x == classes[0], nodes_types))
            class_two_count = sum(map(lambda x: x == classes[1], nodes_types))

            if ((class_one_count / (class_two_count + class_one_count)) * 100) > 0:
                # The Wiener index of a graph is the sum of the shortest-path distances between each pair of reachable nodes
                W_main = nx.wiener_index(comp)  # , weight='weight')
                if class_one_count > 0 and class_two_count > 0:
                    W = W_main * (1 / class_one_count) * (1 / (class_two_count))
                elif class_one_count == 0 and class_two_count > 0:
                    W = W_main * (1 / (class_two_count))

                elif class_two_count == 0 and class_one_count > 0:
                    W = W_main * (1 / (class_one_count))
                else:
                    W = 0
                wiener_index_sum_nodes += W

                wiener_index_sum_edges += W_main * (1 / (comp.number_of_edges()))

    props_data.append(wiener_index_sum_nodes)
    props_data.append(wiener_index_sum_edges)
    return props_data


def mst_two_class(G, classes, node_size_flag=False):
    print('Start MST calculating props!')
    total_weights_sum = 0
    total_weights_sum_pure = 0
    props_data = []
    average_clustering_coefficient = 0

    for component in nx.connected_components(G):
        comp = G.subgraph(list(component))
        # average_clustering_coefficient += nx.algorithms.cluster.average_clustering(comp)
        edges_mst = minimum_spanning_edges_function(comp)

        # print('Edges mst length: {}'.format(len(edges_mst)))
        weights_sum = 0
        weights_sum_pure = 0
        if node_size_flag:

            nodes_types = list(nx.get_node_attributes(comp, 'class_name').values())

            nodes_sizes = list(nx.get_node_attributes(comp, 'node_size').values())

            class_one_indeces = list(map(lambda x: x == classes[0], nodes_types))

            class_one_sizes = [value for i, value in zip(range(0, len(class_one_indeces)), nodes_sizes) if
                               class_one_indeces[i]]

            nodes_sizes_sum = sum(nodes_sizes)
            nodes_one_sum = sum(class_one_sizes)

            # have to consider the nodes size
            if ((nodes_one_sum / nodes_sizes_sum) * 100) > 10:
                if len(edges_mst) > 0:

                    for edge in edges_mst:
                        u_predictions = list(G.nodes[edge[0]]['predictions'].values())
                        u_pred = u_predictions[u_predictions.index(max(u_predictions))]

                        v_predictions = list(G.nodes[edge[1]]['predictions'].values())
                        v_pred = v_predictions[v_predictions.index(max(v_predictions))]

                        u_size = G.nodes[edge[0]]['node_size']
                        v_size = G.nodes[edge[1]]['node_size']

                        # give the prediction less prediction if its from class two + node_size
                        if G.nodes[edge[0]]['class_name'] == classes[1]:
                            u_pred = u_pred / 2
                            u_size = u_size / 2

                        if G.nodes[edge[1]]['class_name'] == classes[1]:
                            v_pred = v_pred / 2
                            v_size = v_size / 2

                        if G.edges[edge]['weight'] > 0:
                            weights_sum += (u_size * u_pred) * (v_size * v_pred) * (1 / G.edges[edge]['weight'])

        else:
            nodes_types = list(nx.get_node_attributes(comp, 'class_name').values())

            class_one_count = sum(map(lambda x: x == classes[0], nodes_types))
            class_two_count = sum(map(lambda x: x == classes[1], nodes_types))

            if ((class_one_count / (class_two_count + class_one_count)) * 100) > 10:

                if len(edges_mst) > 0:

                    for edge in edges_mst:
                        weights_sum_pure += G.edges[edge]['weight']

                        u_predictions = list(G.nodes[edge[0]]['predictions'].values())
                        u_pred = u_predictions[u_predictions.index(max(u_predictions))]

                        v_predictions = list(G.nodes[edge[1]]['predictions'].values())
                        v_pred = v_predictions[v_predictions.index(max(v_predictions))]

                        # give the prediction less prediction if its from class two + node_size
                        if G.nodes[edge[0]]['class_name'] == classes[1]:
                            u_pred = u_pred / 2

                        if G.nodes[edge[1]]['class_name'] == classes[1]:
                            v_pred = v_pred / 2

                        if G.edges[edge]['weight'] > 0:
                            weights_sum += (u_pred) * (v_pred) * (1 / G.edges[edge]['weight'])

        total_weights_sum += weights_sum
        total_weights_sum_pure += weights_sum_pure

    props_data.append(total_weights_sum)
    props_data.append(total_weights_sum_pure)
    # props_data.append(average_clustering_coefficient)
    return props_data


def Shannon_Wiener_Index_one_class(G, class_name):
    '''species diversity
    H ln of number of each class with their size and prediction in each components
    Hmax = ln of number of species (components)
    '''
    print('Start Shannon Wiener Index calculating props!')
    H = 0
    H_ = 0
    _H_ = 0
    E = 0
    E_ = 0
    _E_ = 0
    props_data = []

    class_count_G = len(G)
    nodes_predictions_all = list(nx.get_node_attributes(G, 'predictions').values())

    class_pred_sum_all = sum([value[class_name] for value in nodes_predictions_all])

    comp_count = 0
    for component in nx.connected_components(G):
        comp = G.subgraph(list(component))

        comp_count += 1

        nodes_predictions = list(nx.get_node_attributes(comp, 'predictions').values())

        class_pred_sum = sum([value[class_name] for value in nodes_predictions])

        class_count = len(comp)

        H += class_count * (np.log(class_count) * (class_pred_sum / class_count))

        H_ += (class_count / class_count_G) * (np.log((class_count / class_count_G)) * (class_pred_sum / class_count))

        _H_ += (class_count / class_count_G) * (
                    np.log((class_count / class_count_G)) * (class_pred_sum / class_pred_sum_all))

    Hmax = np.log(comp_count)

    if Hmax != 0:
        E = np.abs(H) / Hmax
        E_ = np.abs(H_) / Hmax
        _E_ = np.abs(_H_) / Hmax
    else:
        E = 0
        E_ = 0
        _E_ = 0
    props_data.append(H)
    props_data.append(H_)
    props_data.append(_H_)
    props_data.append(E)
    props_data.append(E_)
    props_data.append(_E_)
    return props_data


def wiener_index_one_classes(G):
    print('Start wiener_index calculating props!')
    wiener_index_sum_nodes = 0
    wiener_index_sum_edges = 0
    props_data = []
    solo_count = 0

    for component in nx.connected_components(G):
        if len(list(component)) > 1:
            comp = G.subgraph(list(component))
            class_one_count = comp.number_of_nodes()
            W = nx.wiener_index(comp)
            W_count = W * (1 / (class_one_count))
            wiener_index_sum_nodes += W_count
            wiener_index_sum_edges += W * (1 / (comp.number_of_edges()))

        else:
            solo_count += 1

    props_data.append(wiener_index_sum_nodes)
    props_data.append(wiener_index_sum_edges)
    props_data.append(solo_count)
    return props_data


def mst_one_class(G, node_size_flag=False):
    print('Start MST calculating props!')
    total_weights_sum_mst = 0
    total_mst_prediction_mult = 0
    props_data = []
    average_clustering_coefficient = 0

    for component in nx.connected_components(G):
        comp = G.subgraph(list(component))
        # average_clustering_coefficient += nx.algorithms.cluster.average_clustering(comp)
        edges_mst = minimum_spanning_edges_function(comp)
        weights_sum = 0
        mst_prediction_mult = 0
        if node_size_flag:

            if len(edges_mst) > 0:

                for edge in edges_mst:

                    u_predictions = list(G.nodes[edge[0]]['predictions'].values())
                    u_pred = u_predictions[u_predictions.index(max(u_predictions))]

                    v_predictions = list(G.nodes[edge[1]]['predictions'].values())
                    v_pred = v_predictions[v_predictions.index(max(v_predictions))]

                    u_size = G.nodes[edge[0]]['node_size']
                    v_size = G.nodes[edge[1]]['node_size']

                    mst_prediction_mult += (u_size * u_pred) * (v_size * v_pred)

                    if G.edges[edge]['weight'] > 0:
                        weights_sum += (u_size * u_pred) * (v_size * v_pred) * (1 / G.edges[edge]['weight'])

        else:
            if len(edges_mst) > 0:
                for edge in edges_mst:
                    u_predictions = list(G.nodes[edge[0]]['predictions'].values())
                    u_pred = u_predictions[u_predictions.index(max(u_predictions))]

                    v_predictions = list(G.nodes[edge[1]]['predictions'].values())
                    v_pred = v_predictions[v_predictions.index(max(v_predictions))]

                    mst_prediction_mult += (u_pred) * (v_pred)
                    if G.edges[edge]['weight'] > 0:
                        weights_sum += (u_pred) * (v_pred) * (1 / G.edges[edge]['weight'])

        total_weights_sum_mst += weights_sum
        total_mst_prediction_mult += mst_prediction_mult

    props_data.append(total_weights_sum_mst)
    props_data.append(total_mst_prediction_mult)
    # props_data.append(average_clustering_coefficient)
    return props_data


def show_one_class_slide_histogram(G, wsi_name, param, classes):
    print('Calculating the graph histogram!')
    import collections
    import matplotlib.pyplot as plt

    degree_sequence = sorted([d for n, d in G.degree()], reverse=True)  # degree sequence
    degreeCount = collections.Counter(degree_sequence)
    try:
        deg, cnt = zip(*degreeCount.items())

        # fig, ax = plt.subplots()
        # plt.bar(deg, cnt, width=0.80, color="b")
        #
        # plt.title("Degree Histogram")
        # plt.ylabel("Count")
        # plt.xlabel("Degree")
        # ax.set_xticks([d + 0.4 for d in deg])
        # ax.set_xticklabels(deg)
        # fig.savefig(os.path.join(param.graph_images, wsi_name + '_'+classes+'_histogram_graph' + '.png'),
        #             bbox_inches='tight')
        # plt.close()
        return deg, cnt
    except:
        print('No %s!' % classes)
        pass


def count_classes(G, classes, depths):
    row = []
    for i in depths:
        for c in classes[1:]:
            H = gp.extract_around_class_nodes(G, depth=i, classes=[classes[0], c])
            class_two_graph = gp.extract_subgraphs_based_on_classes(H, [c])
            if len(class_two_graph.nodes()) > 0:
                row.append(len(class_two_graph.nodes()))
            else:
                row.append(0)

    return row


def all_mixed_count_info(param, file_name, classes, depths):
    with open(param.graph_info + file_name + '.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        header = []
        header.append("ID")
        for i in depths:
            for item in classes[1:]:
                header.append(item + "__" + str(i))
        writer.writerow(header)

    # import pandas as pd
    # acquired_data = pd.read_csv(os.path.join(param.survival_data_path, file_name_two+'.csv'),error_bad_lines=False)
    # acquired_data_two = pd.read_csv(os.path.join(param.survival_data_path, 'more_tumour' + '.csv'), error_bad_lines=False)

    for wsi_name_ in sorted(glob.glob(param.graph_info + "/*.gpickle")):
        # try:
        with open(param.graph_info + file_name + '.csv', 'a+', newline='') as csvfile:
            writer = csv.writer(csvfile)
            base = os.path.basename(wsi_name_)
            wsi_name = os.path.splitext(base)[0]
            print(wsi_name)

            row = []
            row.append(wsi_name)
            N = gc.read_graph_from_pickle(param, wsi_name)
            H = gp.extract_subgraphs_based_on_classes(N, [classes[0]])
            # props_data = count_classes(N,classes,depths)
            class_one_area = calculate_area_of_class(H)
            class_one_diameter = 0

            for component in nx.connected_components(H):
                comp = H.subgraph(list(component))
                if (len(comp) > 0):
                    d = nx.algorithms.distance_measures.diameter(comp)
                    class_one_diameter += d

            row.append(class_one_diameter)

            print([class_one_diameter, class_one_area])

            writer.writerow(row)


def all_calculate_choosen_props_one_class(param, file_name, classes):
    if not os.path.exists(param.graph_info + file_name + '.csv'):
        with open(param.graph_info + file_name + '.csv', 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            header = []
            header.append("ID")
            # header.append("Nodes_number")
            # header.append("Nodes_edges")
            header.append("degrees")
            header.append("counts")
            header.append("radius")
            header.append("diameter")

            # header.append("wiener_index_sum_nodes")
            # header.append("wiener_index_sum_edges")
            # header.append("solo_count")
            # header.append("shanon_index_sum_H")
            # header.append("shanon_index_sum_H_")
            # header.append("shanon_index_sum_H__")
            # header.append("shanon_index_sum_e")
            # header.append("shanon_index_sum_e_")
            # header.append("shanon_index_sum_e__")
            # header.append("total_mst_weights_sum")
            # header.append("total_mst_prediction_mult")
            # header.append('average_clustering_coefficient')

            writer.writerow(header)
    import pandas as pd
    acquired_data = pd.read_csv(os.path.join(param.graph_info + file_name + '.csv'), error_bad_lines=False)
    for wsi_name_ in sorted(glob.glob(param.graph_info + "/*.gpickle")):
        # try:
        with open(param.graph_info + file_name + '.csv', 'a+', newline='') as csvfile:
            writer = csv.writer(csvfile)
            base = os.path.basename(wsi_name_)
            wsi_name = os.path.splitext(base)[0]
            print(wsi_name)
            data = acquired_data.loc[acquired_data['ID'] == wsi_name]
            print(data)
            if data.empty:
                row = []
                row.append(wsi_name)
                N = gc.read_graph_from_pickle(param, wsi_name)
                H = gp.extract_subgraphs_based_on_classes(N, classes)

                # row.append(len(H.nodes()))
                # row.append(len(H.edges()))

                if len(H.nodes()) > 0:
                    deg, cnt = show_one_class_slide_histogram(H, wsi_name, param, classes[0])

                    row.append(list(deg))
                    row.append(list(cnt))

                    radius = 0
                    diameter = 0

                    for component in nx.connected_components(H):
                        comp = H.subgraph(list(component))
                        if (len(comp) > 0):
                            d = nx.algorithms.distance_measures.diameter(comp)
                            diameter += d
                            r = nx.algorithms.distance_measures.radius(comp)
                            radius += r

                    row.append(radius)
                    row.append(diameter)

                    # props_data = wiener_index_one_classes(H)
                    # for item in props_data:
                    #     row.append(item)
                    #
                    # props_data = Shannon_Wiener_Index_one_class(H, class_name = classes[0])
                    # for item in props_data:
                    #     row.append(item)
                    #
                    # props_data = mst_one_class(H)
                    # for item in props_data:
                    #     row.append(item)
                print(row)

                writer.writerow(row)

            else:
                pass
    # except:
    #     print('No %s found!' %classes)

def calculate_area_of_class(H):
    # print("Calculating graph area!")
    area_sum = 0
    for c in nx.enumerate_all_cliques(H):
        if len(c) <= 3:
            if len(c) == 3:
                a = float(H.edges[c[0], c[1]]['weight'])
                b = float(H.edges[c[1], c[2]]['weight'])
                c = float(H.edges[c[2], c[0]]['weight'])

                s = (a + b + c) / 2

                # calculate the area
                area_sum += (s * (s - a) * (s - b) * (s - c)) ** 0.5

    return area_sum


def calculate_props_all_coms(param, coms_name='louvain', classes=[], advance=False):
    with open(param.graph_info + 'graph_props_coms_' + coms_name + '.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        header = []
        header.append("slide")
        header.append("number_connected_components")
        header.append("number_of_nodes")
        header.append("number_of_edges")

        header.append("tissue_frag")
        header.append("tissue_sheet")
        header.append("tissue_TB")
        header.append("edges_mst_num_local")
        header.append("edges_mst_num_whole")
        header.append("MST_average")
        header.append("MST_average_Whole")
        header.append("MST_weights_sum_avrg")
        header.append("MST_weights_sum_perc")
        header.append("node_connectivity_sum")
        header.append("average_clustering_coefficient")
        header.append("average_cut_value")
        header.append("average_shortest_path_length_sum")
        header.append("shared_edges_perc")
        header.append("other_nodes_perc")
        header.append("wiener_index_aver_sum")

        header.append("voterank_nodes_number")
        header.append("voterank_nodes_number_other")
        header.append("voterank_nodes_number_perc")

        if len(classes) > 1:
            header.append("subgraphs with others")
            header.append("perc subgraphs with others")
            header.append("perc others nodes")
            header.append("perc main class nodes")
            header.append("destory power")
            header.append("subgraphs_with_other_classes_perc_per_subgraph")

        if advance:
            header.append("eccentricity")
            header.append("diameter ")
            header.append("radius ")
            header.append("eccentricity_averg")
            header.append("diameter_averg")
            header.append("radius_averg")
            header.append("degree_centrality_averg")
            header.append("betweenness_centrality_averg")
            header.append("edge_betweenness_centrality")
            header.append("closeness_centrality_averg")

        writer.writerow(header)

    for wsi_name_ in sorted(glob.glob(param.wsi_path + "/*" + param.wsi_ext)):
        # try:
        with open(param.graph_info + 'graph_props_coms_' + coms_name + '.csv', 'a+', newline='') as csvfile:
            writer = csv.writer(csvfile)
            base = os.path.basename(wsi_name_)
            wsi_name = os.path.splitext(base)[0]
            row = []
            row.append(wsi_name)
            N = gc.read_graph_from_pickle(param, wsi_name)
            G = gp.extract_subgraphs_based_on_classes(N, classes)

            if coms_name == 'louvain':
                coms = gcom.louvain_function(G)
            elif coms_name == 'walktrap':
                coms = gcom.walktrap_function(G)
            elif coms_name == 'girvan_newman':
                coms = gcom.girvan_newman_function(G)
            elif coms_name == 'greedy_modularity':
                coms = gcom.greedy_modularity_function(G)
            elif coms_name == 'label_propagation':
                coms = gcom.label_propagation_function(G)
            else:
                coms = gcom.louvain_function(G)

            if advance:
                props_data, props_data_advance = calculate_props_of_communities(coms, G, classes, advance=advance)
                for item in props_data:
                    row.append(item)

                for item in props_data_advance:
                    row.append(item)
            else:
                props_data = calculate_props_of_communities(coms, G, classes, advance=advance)
                for item in props_data:
                    row.append(item)

            writer.writerow(row)

    # except:
    #     print('Problem with the slide:%s!' %(wsi_name))
    #     pass


def nodes_number_all_classes_with_depths(param, file_name, studied_class='', depths=None):
    if depths is None:
        depths = [1, 2, 3, 4]
    with open(os.path.join(param.survival_data_path, file_name + '.csv'), 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        header = []
        header.append("ID")

        # for depth in depths:
        #     header.append( 'stroma_' + str(depth) + "_number_of_nodes")
        #         # header.append(class_name + '_' + str(depth) + "_number_of_edges")

        for class_name in param.class_names:
            if class_name != 'lobules' and class_name != 'others':
                header.append(class_name + "_number_connected_components")
                header.append(class_name + "_number_of_nodes")
                # header.append(class_name+"_number_of_edges")
                if class_name != studied_class and class_name != 'lobules' and class_name != 'others':
                    for depth in depths:
                        header.append(studied_class + '_' + class_name + '_' + str(depth) + "_number_of_nodes")
                        # header.append(class_name + '_' + str(depth) + "_number_of_edges")

        writer.writerow(header)
    print(header)
    import pandas as pd
    acquired_data = pd.read_csv(os.path.join(param.survival_data_path, file_name + '.csv'))
    for wsi_name_ in sorted(glob.glob(param.graph_info + "*.gpickle"), reverse=False):
        with open(os.path.join(param.survival_data_path, file_name + '.csv'), 'a+', newline='') as csvfile:
            writer = csv.writer(csvfile)
            base = os.path.basename(wsi_name_)
            wsi_name = os.path.splitext(base)[0]

            data = acquired_data.loc[acquired_data['ID'] == wsi_name]
            print(data)

            if data.empty:  # and os.path.exists(os.path.join(param.wsi_path, wsi_name + param.wsi_ext)):

                row = []

                row.append(wsi_name)

                print(wsi_name)
                Y = gc.read_graph_from_pickle(param, wsi_name)

                # N = gp.select_ROIs(param, Y,node_degree=node_degree, radius=radius, num_ROI=ROIs)
                # G = gp.extract_subgraphs_based_on_classes(Y, ['tumour','stroma','lym','necrosis'])
                # if len(list(G)) > 0:
                #     largest_cc = max(nx.connected_components(G), key=len)
                # else:
                #     largest_cc = 0
                # N = Y.copy()
                # for u in Y.nodes():
                #     if u in largest_cc:
                #         pass
                #     else:
                #         N.remove_node(u)
                # N = N.subgraph(largest_cc).copy()

                # tumour_count = 0
                # G = gp.extract_subgraphs_based_on_classes(Y, ['tumour'])
                #
                # for component in list(nx.connected_components(G)):
                #     if len(component) > 3:
                #         tumour_count +=1
                #
                # row.append(tumour_count)
                # row.append(G.number_of_nodes())
                # row.append(G.number_of_edges())
                # for depth in depths:
                #     H = gp.extract_around_class_nodes(N, depth=depth, classes=['lym', 'stroma'])
                #     D = gp.extract_subgraphs_based_on_classes(N, classes=['lym'])
                #
                #     row.append(H.number_of_nodes()-D.number_of_nodes())

                for class_name in param.class_names:
                    if class_name != 'lobules' and class_name != 'others':

                        G = gp.extract_subgraphs_based_on_classes(Y, [class_name])
                        row.append(nx.algorithms.components.number_connected_components(G))
                        row.append(G.number_of_nodes())

                        # print([class_name,G.number_of_nodes()])
                        # row.append(G.number_of_edges())
                        if class_name != studied_class and class_name != 'lobules' and class_name != 'others':
                            for depth in depths:
                                T = gp.extract_around_class_nodes(Y, depth=depth, classes=[studied_class, class_name])
                                # from graph_analysis import graph_visualization as gv
                                # gv.display_save_graph_cv(param, wsi_name, Y, graph_name=wsi_name, display=True)
                                # D = gp.extract_subgraphs_based_on_classes(Y, classes=['tumour'])
                                Z = gp.extract_subgraphs_based_on_classes(T, classes=[class_name])

                                row.append(Z.number_of_nodes())
                # for depth in depths:
                #     T = gp.extract_around_class_nodes(Y, depth=depth, classes=['tumour', 'stroma'])
                #     # from graph_analysis import graph_visualization as gv
                #     # gv.display_save_graph_cv(param, wsi_name, Y, graph_name=wsi_name, display=True)
                #     # D = gp.extract_subgraphs_based_on_classes(Y, classes=['tumour'])
                #     Z = gp.extract_subgraphs_based_on_classes(T, classes=['stroma'])
                #
                #     row.append(Z.number_of_nodes())

                print(row)
                writer.writerow(row)


def extract_tumour_bud_all_slides(param, file_name):
    with open(os.path.join(param.survival_data_path, file_name + '.csv'), 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        header = []
        header.append("ID")

        # header.append("tumourBud3")
        # header.append("tumourBud2")
        # header.append("tumourBud1")

        header.append("inner_stroma_cc")
        header.append("inner_stroma_nodes")

        header.append("outer_stroma_cc")
        header.append("outer_stroma_nodes")

        header.append("inner_lym_cc")
        header.append("inner_lym_nodes")

        header.append("outer_lym_cc")
        header.append("outer_lym_nodes")

        writer.writerow(header)

    import pandas as pd
    acquired_data = pd.read_csv(os.path.join(param.survival_data_path, file_name + '.csv'))
    for wsi_name_ in sorted(glob.glob(param.graph_info + "*.gpickle"), reverse=False):
        with open(os.path.join(param.survival_data_path, file_name + '.csv'), 'a+', newline='') as csvfile:

            writer = csv.writer(csvfile)
            base = os.path.basename(wsi_name_)
            wsi_name = os.path.splitext(base)[0]

            data = acquired_data.loc[acquired_data['ID'] == wsi_name]
            print(data)

            if data.empty and os.path.exists(os.path.join(param.wsi_path, wsi_name + param.wsi_ext)):
                row = []

                row.append(wsi_name)

                print(wsi_name)
                Y = gc.read_graph_from_pickle(param, wsi_name)
                # budsList = gp.extract_tumour_buds(Y)
                # row.append(budsList[0])
                # row.append(budsList[1])
                # row.append(budsList[2])

                inner_subgraphs, outer_subgraphs = gp.extract_inner_outer_tissue(Y, class_one='tumour',
                                                                                 class_two='stroma')

                row.append(nx.algorithms.components.number_connected_components(inner_subgraphs))
                row.append(inner_subgraphs.number_of_nodes())

                row.append(nx.algorithms.components.number_connected_components(outer_subgraphs))
                row.append(outer_subgraphs.number_of_nodes())

                inner_subgraphs, outer_subgraphs = gp.extract_inner_outer_tissue(Y, class_one='tumour',
                                                                                 class_two='lym')

                row.append(nx.algorithms.components.number_connected_components(inner_subgraphs))
                row.append(inner_subgraphs.number_of_nodes())

                row.append(nx.algorithms.components.number_connected_components(outer_subgraphs))
                row.append(outer_subgraphs.number_of_nodes())

                print(row)
                writer.writerow(row)
