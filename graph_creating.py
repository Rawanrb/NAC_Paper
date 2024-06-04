import numpy as np
import networkx as nx
from graph_analysis import graph_visualization as gv
from scipy.spatial import Delaunay, KDTree
from collections import defaultdict
import pandas as pd
from scipy.spatial import distance
import glob
import os
import csv

from scipy.cluster.hierarchy import fcluster
from scipy.cluster import hierarchy
from sklearn.neighbors import KDTree as sKDTree
from sklearn.neighbors import kneighbors_graph, radius_neighbors_graph

scale_factor = 16


def connectClusters(Cc, data, dthresh=1000):
    tess = Delaunay(Cc)
    neighbors = defaultdict(set)
    for simplex in tess.simplices:
        for idx in simplex:
            other = set(simplex)
            other.remove(idx)
            neighbors[idx] = neighbors[idx].union(other)
    nx = neighbors

    W = np.zeros((Cc.shape[0], Cc.shape[0]))

    for n in nx:
        nx[n] = np.array(list(nx[n]))
        nx[n] = nx[n][KDTree(Cc[nx[n], :]).query_ball_point(Cc[n], r=dthresh)]

        point_one = [[int(data.iloc[n]['x']), int(data.iloc[n]['y'])]]
        for i in nx[n]:
            point_two = [[data.iloc[i]['x'], data.iloc[i]['y']]]
            W[n, i] = int(distance.euclidean(point_one, point_two))
            W[i, n] = int(distance.euclidean(point_one, point_two))

    return W  # neighbors of each cluster and an affinity matrix


def save_adj_matrix_based_on_csv(param, data, wsi_name, save_matrix=False):
    print('Saving adjacency matrix!')
    X = []

    for index, t in data.iterrows():
        x, y = int(t['x']), int(t['y'])
        X.append([x, y])


    W = connectClusters(np.array(X), data)
    if save_matrix:
        np.save(os.path.join(param.graph_info, wsi_name + '_adj_matrix.npy'), W)


def save_adj_matrix_based_on_csv_eculidean_based(param, data, wsi_name, save_matrix=False):
    print('Saving adjacency matrix!')
    X = []
    for index, t in data.iterrows():
        x, y = int(t['x']), int(t['y'])  # int(abs(t['y'] - large_h))
        X.append([x, y])

    X = np.array(X)
    W = distance.cdist(X, X, 'euclidean')

    for i in range(0, X.shape[0] - 1):
        for j in range(0, X.shape[0] - 1):
            if W[i, j] > 500:
                W[i, j] = 0
            else:
                W[i, j] = int(W[i, j])

    graph_path = os.path.join(param.graph_info, wsi_name)
    os.makedirs(graph_path, exist_ok=True)
    if save_matrix:
        np.save(os.path.join(graph_path, 'adj_matrix.npy'), W)
    print("The data has been saved on: %s adj_matrix.npy" % graph_path)


def read_adj_matrix_to_graph(param, wsi_name):
    W = np.load(os.path.join(param.graph_info, wsi_name + '_adj_matrix.npy'), allow_pickle=True)
    print('The data has been loaded from adjacency matrix!')

    return nx.from_numpy_matrix(W)


def set_graph_attributes_from_csv(param, data, G):
    attrs = {}
    for index, t in data.iterrows():
        x, y = int(t['x']), int(t['y'])  # int(abs(t['y'] - large_h))


        p = t['predictions']

        if isinstance(p, str):
            patch_probs = list(p.strip('][').split(' '))
            props_list = []
            for item in patch_probs:

                item = item.strip("\n")

                props_list.append(item)
            # to get rid of empty items
            test_list = [i for i in props_list if i]
            props_list = test_list
        else:
            props_list = p

        d = dict(zip(param.class_names, props_list))

        attrs[index] = {'coords': [x, y], 'predictions': d, 'class_name': t['class'],
                        'graph_color': param.graph_color[param.class_names.index(t['class'])],
                        'cv_color': param.color_code[param.class_names.index(t['class'])]}

    nx.set_node_attributes(G, attrs)
    print('Nodes data has been saved to the graph!')
    return G


def update_nodes_attributes(param, G):
    for node in G.nodes():
        predictions = G.nodes[node]['predictions'].values()
        predictions_ = []
        for item in predictions:
            predictions_.append(float(item))

        G.nodes[node]['predictions'] = dict(zip(param.class_names, predictions_))

        predictions = G.nodes[node]['predictions'].values()

        predictions_ = []
        for item in predictions:
            predictions_.append(float(item))

        G.nodes[node]['predictions'] = dict(zip(param.class_names, predictions_))

        predictions = G.nodes[node]['predictions'].values()

        max_index = list(predictions).index(max(list(predictions)))

        G.nodes[node]['class_name'] = param.class_names[max_index]
        G.nodes[node]['graph_color'] = param.graph_color[max_index]
        G.nodes[node]['cv_color'] = param.color_code[max_index]

    return G




def set_edges_labels_based_on_nodes_classes(param, G):
    for u, v, a in G.edges(data=True):
        predictions = G.nodes[u]['predictions'].values()
        if len(predictions) > 6:
            print(predictions)

        G.nodes[u]['class_name'] = param.class_names[np.argmax(list(predictions))]

        predictions = G.nodes[v]['predictions'].values()

        G.nodes[v]['class_name'] = param.class_names[np.argmax(list(predictions))]

        G.edges[u, v]['label'] = G.nodes[u]['class_name'] + '_' + G.nodes[v]['class_name']

    print('Edges data has been saved to the graph!')

    return G


def save_graph_to_pickle(param, G, wsi_name):
    if not (os.path.exists(os.path.join(param.graph_info))):
        os.makedirs(os.path.join(param.graph_info), exist_ok=True)
    nx.write_gpickle(G, os.path.join(param.graph_info, wsi_name + ".gpickle"))
    print('The graph data has been saved to pickle file!')
    print(os.path.join(param.graph_info, wsi_name + ".gpickle"))


def read_graph_from_pickle(param, wsi_name):
    print('The graph data has been loaded!')
    print(os.path.join(param.graph_info, wsi_name + ".gpickle"))
    return nx.read_gpickle(os.path.join(param.graph_info, wsi_name + ".gpickle"))


def read_graph_from_np(param, wsi_name):
    graph_path = os.path.join(param.graph_info)
    print('The graph data has been loaded!')
    F = np.load(os.path.join(graph_path, wsi_name + ".npz"))
    X, S, C, W = F['X'], F['S'], F['C'], F['W']
    return X, S, C, W


def read_graph_from_np_patch_based(param, wsi_name):
    graph_path = os.path.join(param.graph_info_patch_based, wsi_name)
    print('The graph data has been loaded!')
    F = np.load(os.path.join(graph_path, wsi_name + ".npz"))
    X, S, C, W = F['X'], F['S'], F['C'], F['W']
    return X, S, C, W


def read_csv_file(param, wsi_name, folder=''):
    print('Reading CSV data for slide:' + wsi_name)
    print(os.path.join(param.prediction_data, folder, wsi_name + '.csv'))
    data = pd.read_csv(os.path.join(param.prediction_data, folder, wsi_name + '.csv'), error_bad_lines=False)
    return data


def read_csv_file_patch_based(param, wsi_name):
    print('Reading CSV data for slide:' + wsi_name)
    data = pd.read_csv(os.path.join(param.patch_based_tile_info, wsi_name + '.csv'), error_bad_lines=False)
    return data


def create_graph(param, wsi_name):
    print('Start creating a graph for ' + wsi_name)
    data = read_csv_file(param, wsi_name)
    save_adj_matrix_based_on_csv(param, data, wsi_name, save_matrix=True)
    G = read_adj_matrix_to_graph(param, wsi_name)
    G = set_graph_attributes_from_csv(param, data, G)
    G = set_edges_labels_based_on_nodes_classes(param, G)
    save_graph_to_pickle(param, G, wsi_name)
    os.remove(os.path.join(param.graph_info, wsi_name + '_adj_matrix.npy'))


def create_graphs_for_all_wsis(param, save_graph=False, overwrite=False, neighboring=False):  # 04-11-071.701.EX.1A

    for wsi_path in sorted(glob.glob(param.wsi_path + "/*" + param.wsi_ext), reverse=True):
       #try:
            base = os.path.basename(wsi_path)
            wsi_name = os.path.splitext(base)[0]
            # wsi_with_folder = 'centroid/'+wsi_name
            wsi_with_folder = wsi_name

            wsi_graph_path = os.path.join(param.graph_info, wsi_with_folder + '.gpickle')

            if overwrite:
                print('The data will be overwrittien!')
                create_graph(param, wsi_with_folder)
            if not (os.path.exists(wsi_graph_path)):
                create_graph(param, wsi_with_folder)
            else:
                print('The graph for slide %s already processed!' % wsi_name)
            if neighboring:
                print('Change the nodes prediction based on the neighboring!')
                from graph_analysis import graph_prediction_modification as gm
                G = gm.wsi_prediction_change_based_on_neighbors(param, wsi_path, wsi_name, save_graph=False)
                save_graph_to_pickle(param, G, wsi_name)

            if save_graph:
                print('Save wsi overlay image!')
                G = read_graph_from_pickle(param, wsi_with_folder)

                gv.display_save_graph_cv(param, wsi_path, wsi_name, G, save_graph=True, graph_name=wsi_name,
                                         display=False)
        # except:
        #     pass


def wsicluster(C, F, lambda_d=3e-3, lambda_f=1e-3, lambda_h=0.85):
    TC = sKDTree(C)
    I, D = TC.query_radius(C, r=6 / lambda_d, return_distance=True, sort_results=True)
    DX = np.zeros(int(C.shape[0] * (C.shape[0] - 1) / 2))
    idx = 0
    for i in range(C.shape[0] - 1):
        f = np.exp(-lambda_f * np.linalg.norm(F[i] - F[I[i]], axis=1))
        d = np.exp(-lambda_d * D[i])
        df = 1 - f * d
        dfi = np.ones(C.shape[0])
        dfi[I[i]] = df
        dfi = dfi[i + 1:]
        DX[idx:idx + len(dfi)] = dfi
        idx = idx + len(dfi)

    Z = hierarchy.linkage(DX, method='average')
    clusters = fcluster(Z, lambda_h, criterion='distance') - 1

    return clusters, groupByCluster(clusters, C, F)


def groupByCluster(clusters, coords, x):
    C, X, S, I = [], [], [], []
    for c in set(clusters):
        idx = clusters == c
        I.append(np.nonzero(idx)[0])
        xm = np.mean(x[idx], axis=0)
        cm = np.mean(coords[idx], axis=0)
        C.append(cm)
        X.append(xm)
        S.append(np.sum(idx))
    return np.array(C), np.array(X), np.array(S), np.array(I)


def conectMetaClusters(clusters, I, Wc):
    Nc = len(I)
    Wcc = np.zeros((Nc, Nc))
    for cid in range(Nc):
        nodes = I[cid]
        connections = Wc[nodes]
        union = np.nonzero(np.sum(connections, axis=0))[0]
        cluster_connections = np.unique(clusters[union])
        Wcc[cid, cluster_connections] = 1.0
        Wcc[cluster_connections, cid] = 1.0
    np.fill_diagonal(Wcc, 0)

    return Wcc

def changing_props_based_on_neighbors(param, wsi_name, class_name):  # 04-08-036.909.EX.1A
    G = read_graph_from_pickle(param, wsi_name)
    from graph_analysis import graph_processing as gp
    H = gp.extract_subgraphs_based_on_classes(G, class_name)
    for component in nx.connected_components(H):
        comp = G.subgraph(list(component))
        node_count = len(comp)
        count = 0
        '''change necrosis based in two conditions:
            - if it is solo that means it should be surrounded by tumour 
            - if it is less than 10 nodes the tumour should be more than 30 percent at least 
            - should be surrounded by all stroma or others 
        '''
        if class_name[0] == 'necrosis':
            for u, v, a in G.edges(comp.nodes(), data=True):
                # in case the node is not surrounded by all tumour then change it to the node class to the next max props
                if 'tumour' not in G.edges[u, v]['label']:
                    count += 1

            if node_count == 1:
                if (count / comp.number_of_edges()) * 100 < 90:
                    G.nodes[comp.nodes()]['predictions'][np.argmax(G.nodes[comp.nodes()]['predictions'])] = 0

            if node_count < 10:
                if (count / comp.number_of_edges()) * 100 < 30:
                    for node in comp.nodes():
                        G.nodes[node]['predictions'][np.argmax(G.nodes[node]['predictions'])] = 0

        '''in case of lobules:
            - if it is one node then it is tumour which most likely to be the second in prediction
            - if comes as many then should be never surrounded by others'''

        if class_name[0] == 'lobules':
            if node_count == 1:
                print(node_count)
                for node in comp.nodes():
                    G.nodes[node]['predictions']['lobules'] = 0

                    predictions = G.nodes[node]['predictions'].values()

                    G.nodes[node]['class_name'] = param.class_names[np.argmax(list(predictions))]


            elif node_count < 10:
                print(node_count)
                for node in comp.nodes():
                    G.nodes[node]['predictions']['lobules'] = 0

                    predictions = G.nodes[node]['predictions'].values()

                    G.nodes[node]['class_name'] = param.class_names[np.argmax(list(predictions))]


    save_graph_to_pickle(param, G, wsi_name)


def from_graph_to_csv(param):
    graph_folder = param.graph_info
    csv_folder = param.prediction_data
    for wsi_name_ in sorted(glob.glob(graph_folder + "*.gpickle"), reverse=True):
        base = os.path.basename(wsi_name_)

        wsi_name = os.path.splitext(base)[0]
        print(wsi_name)
        if not os.path.exists(csv_folder + wsi_name + '.csv'):
            try:
                G = nx.read_gpickle(os.path.join(graph_folder, wsi_name + ".gpickle"))

                with open(csv_folder + wsi_name + '.csv', 'w', newline='') as csvfile:
                    writer = csv.writer(csvfile)
                    header = []
                    header.append('x')
                    header.append('y')
                    header.append('predictions')
                    header.append("class")
                    writer.writerow(header)

                    for u, a in G.nodes(data=True):
                        row = []
                        x = a['coords'][0]
                        y = a['coords'][1]
                        predictions = G.nodes[u]['predictions'].values()

                        graph_class = param.class_names[np.argmax(list(predictions))]

                        prediction = list(predictions)

                        row.append(x)
                        row.append(y)
                        row.append(prediction)
                        row.append(graph_class)
                        writer.writerow(row)


            except:
                print("Error with data of " + wsi_name)


def update_graph_from_csv(param):
    graph_folder = param.graph_info

    for wsi_name_ in sorted(glob.glob(param.prediction_data + "*.csv"), reverse=True):
        # try:
        base = os.path.basename(wsi_name_)
        wsi_name = os.path.splitext(base)[0]
        print(wsi_name)
        data = read_csv_file(param, wsi_name)
        G = nx.read_gpickle(os.path.join(graph_folder, wsi_name + ".gpickle"))
        G = set_graph_attributes_from_csv(param, data, G)
        G = set_edges_labels_based_on_nodes_classes(param, G)
        save_graph_to_pickle(param, G, wsi_name)
    # except:
    #     pass


def _csv(param):
    for wsi_name_ in sorted(glob.glob(param.graph_info + "*.gpickle"), reverse=True):
        base = os.path.basename(wsi_name_)

        wsi_name = os.path.splitext(base)[0]
        print(wsi_name)
        data = pd.read_csv(os.path.join(param.prediction_data, wsi_name + '_.csv'), error_bad_lines=False)

        with open(param.prediction_data + wsi_name + '.csv', 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            header = []
            header.append('x')
            header.append('y')
            header.append('predictions')
            header.append("class")
            writer.writerow(header)

            for index, t in data.iterrows():
                row = []
                row.append(t['x'])
                row.append(t['y'])
                p = t['predictions']

                if isinstance(p, str):
                    patch_probs = list(p.strip('][').split(','))
                    props_list = []
                    for item in patch_probs:
                        # item_ = item.rstrip("/n")
                        # item_ = item.strip("'")
                        item1 = item.replace("'", '')
                        item2 = item1.replace("'", '')
                        item3 = item2.replace(' ', '')

                        props_list.append(float(item3))

                    # # to get rid of empty items
                    # test_list = [i for i in props_list if i]
                    # props_list = test_list

                row.append(props_list)
                row.append(t['class'])

                writer.writerow(row)