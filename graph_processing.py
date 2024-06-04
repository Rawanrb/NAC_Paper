
import networkx as nx

scale_factor = 16

def extract_tumour_buds(G,classes=['tumour']):
    tumourBudCount3 = 0
    tumourBudCount2 = 0
    tumourBudCount1 = 0
    # [G.remove_node(node) for node, degree in dict(G.degree()).items() if degree < 5]

    N = extract_subgraphs_based_on_classes(G, classes)

    for component in list(nx.connected_components(N)):
        if len(component) < 4:
            tumourBudCount3 +=1
            if len(component) < 3:
                tumourBudCount2 += 1
                if len(component) < 2:
                    tumourBudCount1 += 1

    return [tumourBudCount3,tumourBudCount2,tumourBudCount1]




def extract_subgraphs_based_on_condition(H,classes=[],condition=''):
    sub_graphs = nx.Graph()

    for c in nx.connected_components(H):
        if nx.graph_number_of_cliques(H.subgraph(c)) > 100:
            # nx.is_chordal(H.subgraph(c))
            # len(H.subgraph(c).nodes())>10 and len(nx.cycle_basis(H.subgraph(c)))>200 and (len(H.subgraph(c).nodes())/len(nx.cycle_basis(H.subgraph(c))))<1:
            sub_graphs.add_edges_from(H.subgraph(c).edges(data=True))
            sub_graphs.add_nodes_from(H.subgraph(c).nodes(data=True))
            #print("radius: %d" % nx.radius(H.subgraph(c)))

            for path in nx.all_simple_paths(H.subgraph(c)):
                print(path)

            # print("diameter: %d" % nx.diameter(H.subgraph(c)))
            # #print("eccentricity: %s" % nx.eccentricity(H.subgraph(c)))
            # print("center: %s" % nx.center(H.subgraph(c)))
            # print("periphery: %s" % nx.periphery(H.subgraph(c)))
            # print("density: %s" % nx.density(H.subgraph(c)))

            # large_clique_size
            # graph_number_of_cliques

            # for n in nx.center(H.subgraph(c)):
            #     for m in nx.periphery(H.subgraph(c)):
            #         print(nx.node_connectivity(H.subgraph(c),n,m))
            # print(nx.dijkstra_path(H.subgraph(c),n,m))
            # print(len(H.subgraph(c).nodes()))
            # print("node cycles: %s" % len(nx.cycle_basis(H.subgraph(c))))
            # print("Number of triangles: %s" % nx.triangles(H.subgraph(c)))
            # print("Number of triangles tran: %s" % nx.transitivity(H.subgraph(c)))
            # print("Number of constraint tran: %s" % nx.constraint(H.subgraph(c)))
            # print("Clique:" + str(nx.graph_number_of_cliques(H.subgraph(c))))

            # print(nx.all_pairs_node_connectivity(H.subgraph(c)))
        # print(["Nodes:", H.subgraph(c).nodes(data=True)])
        # print(["Edges:", H.subgraph(c).edges(data=True)])

        # print(G.edges(data=True))

    return sub_graphs

def extract_simple_association_based_whole_graph(G, classes=[], min_connection_percentage=10):

    class_one_subgraphs = extract_subgraphs_based_on_classes(G, classes=classes[0])
    class_two_subgraphs = extract_subgraphs_based_on_classes(G, classes=classes[1])

    for c in nx.connected_components(class_two_subgraphs):
        number_of_outer_edges_per_node = 0
        number_of_related_outer_edges_per_node = 0
        class_two_subgraphs_modif = class_two_subgraphs.subgraph(c).copy()
        for node in class_two_subgraphs.subgraph(c).nodes():
            for u,v in G.edges(node):
                if G.nodes[u]['class_name'] != G.nodes[v]['class_name']:
                    number_of_outer_edges_per_node += 1
                if G.nodes[u]['class_name']==classes[0] and G.nodes[v]['class_name']==classes[1]:
                    class_two_subgraphs_modif.add_edge(u, v, **G.get_edge_data(u, v))

                    attrs = {}
                    attrs[u] = G.nodes[u]
                    attrs[v] = G.nodes[v]
                    nx.set_node_attributes(class_two_subgraphs_modif, attrs)

                    number_of_related_outer_edges_per_node += 1
                if G.nodes[u]['class_name']==classes[1] and G.nodes[v]['class_name']==classes[0]:
                    class_two_subgraphs_modif.add_edge(u, v, **G.get_edge_data(u, v))
                    attrs = {}
                    attrs[u] = G.nodes[u]
                    attrs[v] = G.nodes[v]
                    nx.set_node_attributes(class_two_subgraphs_modif, attrs)

                    number_of_related_outer_edges_per_node += 1


        if (number_of_related_outer_edges_per_node/number_of_outer_edges_per_node)*100 > min_connection_percentage:
            class_one_subgraphs.add_edges_from(class_two_subgraphs_modif.edges(data=True))
            class_one_subgraphs.add_nodes_from(class_two_subgraphs_modif.nodes(data=True))

    return class_one_subgraphs

def extract_simple_association_based_two_subgraphs(G,subgraphs_one,subgraphs_two,min_connection_percentage=10):

    for c in nx.connected_components(subgraphs_two):
        number_of_nodes_per_subgraph = subgraphs_two.subgraph(c).nodes()
        number_of_related_outer_edges_per_subgraph = 0
        for node in subgraphs_two.subgraph(c).nodes():
            for u,v in G.edges(node):
                if u in subgraphs_one.nodes():
                    number_of_related_outer_edges_per_subgraph +=1
                if v in subgraphs_one.nodes():
                    number_of_related_outer_edges_per_subgraph += 1

        if (number_of_related_outer_edges_per_subgraph/number_of_nodes_per_subgraph)*100 > min_connection_percentage:
            subgraphs_one.add_edges_from(subgraphs_two.subgraph(c).edges(data=True))
            subgraphs_one.add_nodes_from(subgraphs_two.subgraph(c).nodes(data=True))

    return subgraphs_one

def extract_sheet_frag_of_tissue(G,classes='tumour',clique_threshold=1000):
    class_graph = extract_subgraphs_based_on_classes(G,classes)

    sheets_subgraphs = nx.Graph()
    frag_subgraphs = nx.Graph()

    for c in nx.connected_components(class_graph):
        if nx.graph_number_of_cliques(class_graph.subgraph(c)) > clique_threshold:

            node_removal_flag = True

            processed_graph_with_degree = class_graph.subgraph(c).copy()

            modified_graph = class_graph.subgraph(c).copy()

            while node_removal_flag:
                node_removal_flag = False
                for key in processed_graph_with_degree.nodes().keys():
                    if processed_graph_with_degree.degree(key) < 3:
                        node_removal_flag = True
                        modified_graph.remove_node(key)

                processed_graph_with_degree = modified_graph.copy()



            for cc in nx.connected_components(processed_graph_with_degree):
                if nx.graph_number_of_cliques(processed_graph_with_degree.subgraph(cc))> clique_threshold:
                    sheets_subgraphs.add_edges_from(processed_graph_with_degree.subgraph(cc).edges(data=True))
                    sheets_subgraphs.add_nodes_from(processed_graph_with_degree.subgraph(cc).nodes(data=True))
                else:
                    frag_subgraphs.add_edges_from(processed_graph_with_degree.subgraph(cc).edges(data=True))
                    frag_subgraphs.add_nodes_from(processed_graph_with_degree.subgraph(cc).nodes(data=True))

        else:
            frag_subgraphs.add_edges_from(class_graph.subgraph(c).edges(data=True))
            frag_subgraphs.add_nodes_from(class_graph.subgraph(c).nodes(data=True))


    return sheets_subgraphs,frag_subgraphs

def extract_inner_outer_tissue(G,class_one='tumour',class_two='lym',max_connection_percentage= 70,min_connection_percentage=10):
    class_graph = extract_subgraphs_based_on_classes(G, classes=[class_one,class_two])


    class_one_subgraphs = extract_subgraphs_based_on_classes(class_graph, classes=[class_one])
    class_two_subgraphs = extract_subgraphs_based_on_classes(class_graph, classes=[class_two])

    inner_subgraphs = class_one_subgraphs.copy()
    outer_subgraphs = class_one_subgraphs.copy()

    #sewing process
    for c in nx.connected_components(class_two_subgraphs):
        number_of_outer_edges_per_node = 0
        number_of_related_outer_edges_per_node = 0
        for node in class_two_subgraphs.subgraph(c).nodes():
            for u,v in G.edges(node):
                if G.nodes[u]['class_name'] != G.nodes[v]['class_name']:
                    number_of_outer_edges_per_node += 1
                if G.nodes[u]['class_name']==class_one and G.nodes[v]['class_name']==class_two:
                    number_of_related_outer_edges_per_node += 1
                if G.nodes[u]['class_name']==class_two and G.nodes[v]['class_name']==class_one:
                    number_of_related_outer_edges_per_node += 1
        #in case the node is solo no connection with any
        if number_of_outer_edges_per_node>0:
            if (number_of_related_outer_edges_per_node/number_of_outer_edges_per_node)*100 > max_connection_percentage:
                inner_subgraphs.add_edges_from(class_two_subgraphs.subgraph(c).edges(data=True))
                inner_subgraphs.add_nodes_from(class_two_subgraphs.subgraph(c).nodes(data=True))

            elif (number_of_related_outer_edges_per_node/number_of_outer_edges_per_node)*100 >min_connection_percentage:
                outer_subgraphs.add_edges_from(class_two_subgraphs.subgraph(c).edges(data=True))
                outer_subgraphs.add_nodes_from(class_two_subgraphs.subgraph(c).nodes(data=True))

    inner = extract_subgraphs_based_on_classes(inner_subgraphs, classes=[class_two])
    outer = extract_subgraphs_based_on_classes(outer_subgraphs, classes=[class_two])
    return inner,outer


def extract_associate_sheet_frag_tissue_with_diff_tissue(G,class_one='tumour',class_two='lym',connection_min_clique_threshold= 300,min_connection_percentage=10):

    class_graph = extract_subgraphs_based_on_classes(G, classes=[class_one, class_two])

    class_two_subgraphs = extract_subgraphs_based_on_classes(class_graph, classes=[class_two])

    sheet_subgraphs, frag_subgraphs = extract_sheet_frag_of_tissue(G, classes=class_one, clique_threshold=connection_min_clique_threshold)



    for c in nx.connected_components(class_two_subgraphs):
        number_of_outer_edges_per_node = 0
        number_of_related_outer_edges_per_node = 0
        related_nodes_tissue_one = []

        class_two_subgraphs_modif = class_two_subgraphs.subgraph(c).copy()

        for node in class_two_subgraphs.subgraph(c).nodes():
            for u, v in G.edges(node):
                if G.nodes[u]['class_name'] != G.nodes[v]['class_name']:
                    number_of_outer_edges_per_node += 1

                if G.nodes[u]['class_name'] == class_one and G.nodes[v]['class_name'] == class_two:
                    class_two_subgraphs_modif.add_edge(u, v, **G.get_edge_data(u, v))

                    attrs = {}
                    attrs[u] = G.nodes[u]
                    attrs[v] = G.nodes[v]
                    nx.set_node_attributes(class_two_subgraphs_modif, attrs)

                    number_of_related_outer_edges_per_node += 1
                    related_nodes_tissue_one.append(u)

                if G.nodes[u]['class_name'] == class_two and G.nodes[v]['class_name'] == class_one:


                    class_two_subgraphs_modif.add_edge(u, v,**G.get_edge_data(u,v))

                    attrs = {}
                    attrs[u] = G.nodes[u]
                    attrs[v] = G.nodes[v]
                    nx.set_node_attributes(class_two_subgraphs_modif, attrs)

                    number_of_related_outer_edges_per_node += 1
                    related_nodes_tissue_one.append(v)


        if number_of_outer_edges_per_node >0:
            if (number_of_related_outer_edges_per_node / number_of_outer_edges_per_node) * 100 > min_connection_percentage:

                #print('The connection is more than %d%s !'%(min_connection_percentage,'%'))


                if any(item in related_nodes_tissue_one for item in frag_subgraphs.nodes()):
                    #print("Found in frag!")
                    frag_subgraphs.add_edges_from(class_two_subgraphs_modif.edges(data=True))
                    frag_subgraphs.add_nodes_from(class_two_subgraphs_modif.nodes(data=True))


                elif any(item in related_nodes_tissue_one for item in sheet_subgraphs.nodes()):
                   # print("Found in sheet!")
                    sheet_subgraphs.add_edges_from(class_two_subgraphs_modif.edges(data=True))
                    sheet_subgraphs.add_nodes_from(class_two_subgraphs_modif.nodes(data=True))


                # else:
                #     print("Found in neigther!")



    return sheet_subgraphs,frag_subgraphs

def extract_associate_sheet_frag_tissue_with_diff_tissue_with_perc(G, class_one='tumour', class_two='lym',connection_min_clique_threshold= 300,max_connection_percentage=70,min_connection_percentage=10):

    sheet_subgraphs_association, frag_subgraphs_association = \
        extract_associate_sheet_frag_tissue_with_diff_tissue(G,class_one=class_one,class_two=class_two,connection_min_clique_threshold= connection_min_clique_threshold,min_connection_percentage=min_connection_percentage)

    sheet_subgraphs_association_huge, sheet_subgraphs_association_little = nx.Graph(),nx.Graph()

    frag_subgraphs_association_huge, frag_subgraphs_association_little = nx.Graph(),nx.Graph()



    for c in nx.connected_components(sheet_subgraphs_association):
        tissue_one_nodes = 0
        tissue_two_nodes = 0
        for n,a in sheet_subgraphs_association.subgraph(c).nodes(data=True):
            if a['class_name'] == class_one:
                tissue_one_nodes += 1
            elif a['class_name'] == class_two:
                tissue_two_nodes += 1

        if (tissue_two_nodes/tissue_one_nodes)*100 > max_connection_percentage:
            sheet_subgraphs_association_huge.add_edges_from(sheet_subgraphs_association.subgraph(c).edges(data=True))
            sheet_subgraphs_association_huge.add_nodes_from(sheet_subgraphs_association.subgraph(c).nodes(data=True))
        else:
            sheet_subgraphs_association_little.add_edges_from(sheet_subgraphs_association.subgraph(c).edges(data=True))
            sheet_subgraphs_association_little.add_nodes_from(sheet_subgraphs_association.subgraph(c).nodes(data=True))



    for c in nx.connected_components(frag_subgraphs_association):
        tissue_one_nodes = 0
        tissue_two_nodes = 0
        for n, a in frag_subgraphs_association.subgraph(c).nodes(data=True):
            if a['class_name'] == class_one:
                tissue_one_nodes += 1
            elif a['class_name'] == class_two:
                tissue_two_nodes += 1

        if (tissue_two_nodes / tissue_one_nodes) * 100 > max_connection_percentage:
            frag_subgraphs_association_huge.add_edges_from(frag_subgraphs_association.subgraph(c).edges(data=True))
            frag_subgraphs_association_huge.add_nodes_from(frag_subgraphs_association.subgraph(c).nodes(data=True))
        else:
            frag_subgraphs_association_little.add_edges_from(frag_subgraphs_association.subgraph(c).edges(data=True))
            frag_subgraphs_association_little.add_nodes_from(frag_subgraphs_association.subgraph(c).nodes(data=True))


    return sheet_subgraphs_association_huge,sheet_subgraphs_association_little,frag_subgraphs_association_huge, frag_subgraphs_association_little


def extract_around_class_nodes(G,depth=1,classes=[]):
    '''The main class nodes extracted and then take its neighbor from other classes
    based on the required depth(how far away the nodes are from the main class)'''
    print('Processing the neighbors based on the depth %d!'%depth)

    N = G.copy()
    for u, a in G.nodes(data=True):
        if a['class_name'] in classes[0]:
            pass
        else:
            N.remove_node(u)

    H = N.copy()
    #extract the neighbor from the other class and set their weights
    attrs = {}
    nodes = N.nodes()
    for i in range(0,depth):
        nodes = N.nodes()
        for node in nodes:
            for n in G.neighbors(node):
                if G.nodes[n]['class_name'] in classes:#classes[1]:
                    H.add_node(n, **G.nodes[n])
                    H.add_edge(node, n)

                    if [node, n] in G.edges():
                        attrs[(node, n)] = {"weight": G.edges[node, n]['weight'],"label":G.edges[node, n]['label']}
                    elif [n, node] in G.edges():
                        attrs[(n, node)] = {"weight": G.edges[n, node]['weight'],"label":G.edges[node, n]['label']}
        N = H.copy()
        nx.set_edge_attributes(H, attrs)

    return H

def change_class_type(param,G,classes,threshold=0.1,replace_no_condition=False):
    class_one = classes[0]
    class_two = classes[1]
    for node in G.nodes():
        predictions = G.nodes[node]['predictions'].values()

        predictions_=[]
        for item in predictions:
            predictions_.append(float(item))

        G.nodes[node]['predictions'] = dict(zip(param.class_names, predictions_))

        predictions = G.nodes[node]['predictions'].values()

        if list(predictions).index(max(list(predictions))) == param.class_names.index(class_one):
            # if the prediction is between the two classes
            data = list(predictions)
            data[list(predictions).index(max(list(predictions)))] = 0

            #in case you want to replace all nodes no matter about their prediction uncomment the lines down

            if replace_no_condition:
                G.nodes[node]['predictions'][class_two] = G.nodes[node]['predictions'][class_one]
                if isinstance(predictions, str):
                    G.nodes[node]['predictions'][class_one] = '0'
                else:
                    G.nodes[node]['predictions'][class_one] = 0

            #in case you want the props to be greater than something uncomment the thresholding
            else:
                if data.index(max(data)) == param.class_names.index(class_two) and max(data) > threshold:
                    G.nodes[node]['predictions'][class_two] = G.nodes[node]['predictions'][class_one]
                    if isinstance(predictions, str):
                        G.nodes[node]['predictions'][class_one] = '0'
                    else:
                        G.nodes[node]['predictions'][class_one] = 0
    return G



def select_ROIs(param,G, node_degree=3,radius=600, num_ROI=3):
    #import math
    import numpy as np
    from graph_analysis import graph_visualization as gv
    pos, graph_colors, cv_color = gv.extract_graph_pos_colors_for_overlay(param, G)
    graph_class = []

    for u, a in G.nodes(data=True):
        graph_class.append(a['class_name'])
    nodelist = list(G)
    xy = np.asarray([pos[v] for v in nodelist])
    #large_w, large_h = gv.slide_dimensions(param, wsi_name)

    import random



    H = extract_subgraphs_based_on_classes(G, ['tumour'])
    #[H.remove_node(node) for node, degree in dict(H.degree()).items() if (degree < node_degree)]
    # H = N.copy()
    # for component in list(nx.connected_components(N)):
    #     if len(component) < 10:
    #         for node in component:
    #             H.remove_node(node)
    percent = 50
    dec_perc = 0
    choosen_node_final = []
    ROI_counter = 0
    loop_counter = 0
    while ROI_counter < num_ROI:  # and counter1 < (10 + num_ROI):
        loop_counter +=1

        if loop_counter == num_ROI:
            dec_perc += 10
            loop_counter = 0

        non_other_num = 0
        choosen_num = 0
        tumour_counter = 0
        choosen_node = []

        # chosen_node = random.choice(tumour_nodes)
        # (center_x, center_y) = xy[chosen_node]
        sample = random.sample(list(H.nodes), 1)

        (center_x, center_y) = xy[sample[0]]


        # center_x = random.sample(range(radius, new_w - radius), 1)  # np.random.uniform(low=0, high=large_w, size=(3,))
        # center_y = random.sample(range(radius, new_h - radius), 1)  # np.random.uniform(low=0, high=large_h, size=(3,))

        for i, (x, y) in enumerate(xy):
            if (x - center_x) ** 2 + (y - center_y) ** 2 <= radius ** 2:
                if graph_class[i] != 'others':
                    non_other_num += 1
                    choosen_node.append(i)
                if graph_class[i] != 'tumour':
                    tumour_counter +=1

                choosen_num += 1


        if ((tumour_counter / choosen_num) * 100) > (percent-dec_perc):
            H.remove_node(sample[0])
            ROI_counter += 1
            for i in choosen_node:
                choosen_node_final.append(i)


        # if counter1 > 10:
        #     counter1 += 1
        #     counter += 1
        #     for i in choosen_node:
        #         choosen_node_final.append(i)
        #
        # elif choosen_num > 0 and ((non_other_num / choosen_num) * 100) > 50:
        #     counter += 1
        #     for i in choosen_node:
        #         choosen_node_final.append(i)

    for i, (x, y) in enumerate(xy):
        if i not in choosen_node_final:
            G.remove_node(i)

    return G

def extract_subgraphs_based_on_classes(G,classes=[]):
    N = G.copy()
    if len(classes)<3:
        for u, a in G.nodes(data=True):
            if a['class_name'] in classes:
                pass
            else:
                N.remove_node(u)
        return N

    elif len(classes)==3:
        for u, a in G.nodes(data=True):
            if a['class_name'] in classes:
                pass
            else:
                N.remove_node(u)

        for u, a in N.nodes(data=True):
            if a['class_name'] ==classes[2]:
                a['class_name'] = classes[1]
            else:
                pass

        return N

    elif len(classes)==4:
        for u, a in G.nodes(data=True):
            if a['class_name'] in classes:
                pass
            else:
                N.remove_node(u)

        for u, a in N.nodes(data=True):
            if a['class_name'] == classes[1]:
                a['class_name'] = classes[0]
            elif a['class_name'] == classes[3]:
                a['class_name'] = classes[2]
            else:
                pass



        return N



def extract_graph_connected_borders_based_on_classes(G,classes=[]):
    N = G.copy()
    for u, a in G.nodes(data=True):

        if a['class_name'] in classes:
            pass
        else:
            N.remove_node(u)
    H = N.copy()
    for u, v in N.edges():
        if N.nodes[u]['class_name'] == N.nodes[v]['class_name']:
            H.remove_edge(u, v)

    [H.remove_node(node) for node, degree in dict(H.degree()).items() if (degree < 1)]

    return H

def extract_subgraphs_based_on_classes_with_distance(G,distance,classes=[]):
    N = G.copy()
    if len(classes)<3:
        for u, a in G.nodes(data=True):
            if a['class_name'] in classes:
                pass
            else:
                N.remove_node(u)

    H = N.copy()
    for u, v, a in N.edges(data=True):
        if a['weight']>distance:
            H.remove_edge(u,v)

    [H.remove_node(node) for node, degree in dict(H.degree()).items() if
     (degree < 1)]

    return H

def extract_around_class_nodes_with_specific_subgraph(G, classes,flower_degree):
    H = G.copy()
    for node in G.nodes():
        comp = G.subgraph(list(G.neighbors(node)))
        nodes_types = list(nx.get_node_attributes(comp, 'class_name').values())

        class_count_one = sum(map(lambda x: x == classes[0], nodes_types))

        if class_count_one == len(list(G.neighbors(node))) and G.nodes[node]["class_name"]==classes[0]:
            H.remove_node(node)

    N = H.copy()
    count = 0
    for node in H.nodes():
        if G.nodes[node]["class_name"]==classes[0]:
            comp = G.subgraph(list(H.neighbors(node)))
            nodes_types = list(nx.get_node_attributes(comp, 'class_name').values())
            class_count_two = sum(map(lambda x: x == classes[1], nodes_types))
            if class_count_two > flower_degree:
                count +=1
            else:
                N.remove_node(node)


    return count



def extract_subgraphs_fayyaz(X, W,class_index = []):
    import numpy as np
    # S is node size,  X is props, C is positions, W is matrix
    G = nx.from_numpy_matrix(W)
    props = dict(zip(range(X.shape[0]), X))

    H= G.copy()

    for u in G.nodes():
        for u, v in G.edges(u):
            max_index_u = np.argmax(props[u])
            max_index_v = np.argmax(props[v])
            if max_index_u in class_index:
                if max_index_v in class_index:
                    if max_index_u == max_index_v and max_index_u == class_index[1]:
                        pass
                    elif max_index_u == max_index_v and max_index_u == class_index[0]:
                        if u in H.nodes() and v in H.nodes() and (u,v) in H.edges():
                            H.remove_edge(u,v)
            else:
                if u in H.nodes():
                    H.remove_node(u)
                    if max_index_v not in class_index:
                        if v in H.nodes():
                            H.remove_node(v)


    [H.remove_node(node) for node, degree in dict(H.degree()).items() if (degree < 1) and (np.argmax(props[node])==class_index[0])]
    return H


def extract_patches_based_on_nodes_coords(param,neighboring_threshold=95,classes='tumour',patch_number=1):
    import glob, os

    from graph_analysis import graph_creating as gc

    for wsi_name_ in sorted(glob.glob(param.wsi_path + "/*.ndpi")):
        try:
            base = os.path.basename(wsi_name_)
            wsi_name = os.path.splitext(base)[0]
            G = gc.read_graph_from_pickle(param, wsi_name)
            from general import slide_no_prefix
            from general import util

            wsi_slide = slide_no_prefix.open_slide(wsi_name_)

            count = 0
            node_list = []
            for u, a in G.nodes(data=True):
                if count < patch_number:
                    if u not in node_list:
                        if a['class_name'] == classes:
                            neighbor = 0
                            total = 0
                            total += len(list(G.neighbors(u)))
                            for n in G.neighbors(u):
                                node_list.append(n)
                                if G.nodes[n]['class_name'] == classes:
                                    neighbor += 1
                                    total += len(list(G.neighbors(n)))
                                    for d in G.neighbors(n):
                                        node_list.append(d)
                                        if G.nodes[d]['class_name'] == classes:
                                            neighbor += 1
                                            total += len(list(G.neighbors(d)))
                                            for c in G.neighbors(d):
                                                if G.nodes[c]['class_name'] == classes:
                                                    neighbor += 1

                            if total > 0 and ((neighbor/ total)*100 > neighboring_threshold):
                                count += 1
                                #if count > 2000 and count < (patch_number + 2000):

                                x = int(((a['coords'][0])))
                                y = int(((a['coords'][1])))

                                tile_region = wsi_slide.read_region((x, y), 1, (224, 224))
                                # RGBA to RGB
                                pil_img = tile_region.convert("RGB")
                                np_img = util.pil_to_np_rgb(pil_img)

                                #patch_path = os.path.join(param.base_dir,'unsupervised_patches',wsi_name+'_'+str(count)+'.png')
                                patch_path = os.path.join(param.base_dir, 'model_patches',classes,
                                                          wsi_name+'_'+str(count)+'.png')

                                dir = os.path.dirname(patch_path)
                                if not os.path.exists(dir):
                                    os.makedirs(dir)
                                    print(dir + ' created!')

                                patch = util.np_to_pil(np_img)
                                patch.save(patch_path)
                                print(patch_path)
                else:
                    break
        except:
            pass