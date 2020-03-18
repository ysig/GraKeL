"""The base file for loading default datasets."""
# Python 2/3 cross-compatibility import
from __future__ import print_function

import os
import shutil
import zipfile
import ssl
try:
    # Python 2
    from urllib2 import HTTPError
    from urllib2 import urlopen
except ImportError:
    # Python 3+
    from urllib.error import HTTPError
    from urllib.request import urlopen

import numpy as np

from shutil import copyfileobj

from collections import Counter

from sklearn.utils import Bunch

from grakel.graph import Graph

global datasets_metadata, symmetric_dataset

dataset_metadata = {
    "AIDS": {"nl": True, "el": True, "na": True, "ea": False,
             "link": "https://ls11-www.cs.uni-dortmund.de/people/" +
             "morris/graphkerneldatasets/AIDS.zip"},
    "BZR": {"nl": True, "el": False, "na": True, "ea": False,
            "link": "https://ls11-www.cs.uni-dortmund.de/people/" +
            "morris/graphkerneldatasets/BZR.zip"},
    "BZR_MD": {"nl": True, "el": True, "na": False, "ea": True,
               "link": "https://ls11-www.cs.uni-dortmund.de/people/" +
               "morris/graphkerneldatasets/BZR_MD.zip"},
    "COIL-DEL": {"nl": False, "el": True, "na": True, "ea": False,
                 "link": "https://ls11-www.cs.uni-dortmund.de/people/morris/" +
                 "graphkerneldatasets/COIL-DEL.zip"},
    "COIL-RAG": {"nl": False, "el": False, "na": True, "ea": True,
                 "link": "https://ls11-www.cs.uni-dortmund.de/people/morris/" +
                 "graphkerneldatasets/COIL-RAG.zip"},
    "COLLAB": {"nl": False, "el": False, "na": False, "ea": False,
               "link": "https://ls11-www.cs.uni-dortmund.de/people/morris/" +
               "graphkerneldatasets/COLLAB.zip"},
    "COX2": {"nl": True, "el": False, "na": True, "ea": False,
             "link": "https://ls11-www.cs.uni-dortmund.de/people/morris/" +
             "graphkerneldatasets/COX2.zip"},
    "COX2_MD": {"nl": True, "el": True, "na": False, "ea": True,
                "link": "https://ls11-www.cs.uni-dortmund.de/people/morris/" +
                "graphkerneldatasets/COX2_MD.zip"},
    "DHFR": {"nl": True, "el": False, "na": True, "ea": False,
             "link": "https://ls11-www.cs.uni-dortmund.de/people/morris/" +
             "graphkerneldatasets/DHFR.zip"},
    "DHFR_MD": {"nl": True, "el": True, "na": False, "ea": True,
                "link": "https://ls11-www.cs.uni-dortmund.de/people/morris/" +
                "graphkerneldatasets/DHFR_MD.zip"},
    "ER_MD": {"nl": True, "el": True, "na": False, "ea": True,
              "link": "https://ls11-www.cs.uni-dortmund.de/people/morris/" +
              "graphkerneldatasets/ER_MD.zip"},
    "DD": {"nl": True, "el": False, "na": False, "ea": False,
           "link": "https://ls11-www.cs.uni-dortmund.de/people/morris/" +
           "graphkerneldatasets/DD.zip"},
    "ENZYMES": {"nl": True, "el": False, "na": True, "ea": False,
                "link": "https://ls11-www.cs.uni-dortmund.de/people/morris/" +
                "graphkerneldatasets/ENZYMES.zip"},
    "Cuneiform": {"nl": True, "el": True, "na": True, "ea": True,
                  "link": "https://ls11-www.cs.uni-dortmund.de/people/morris/"
                          "graphkerneldatasets/Cuneiform.zip"},
    "FINGERPRINT": {"nl": False, "el": False, "na": True, "ea": True,
                    "link": "https://ls11-www.cs.uni-dortmund.de/people/" +
                    "morris/graphkerneldatasets/Fingerprint.zip"},
    "FIRSTMM_DB": {"nl": True, "el": False, "na": True, "ea": True,
                   "link": "https://ls11-www.cs.uni-dortmund.de/people/" +
                   "morris/graphkerneldatasets/FIRSTMM_DB.zip"},
    "FRANKENSTEIN": {"nl": False, "el": False, "na": True, "ea": False,
                     "link": "https://ls11-www.cs.uni-dortmund.de/people/mor" +
                     "ris/graphkerneldatasets/FRANKENSTEIN.zip"},
    "IMDB-BINARY": {"nl": False, "el": False, "na": False, "ea": False,
                    "link": "https://ls11-www.cs.uni-dortmund.de/people/mor" +
                    "ris/graphkerneldatasets/IMDB-BINARY.zip"},
    "IMDB-MULTI": {"nl": False, "el": False, "na": False, "ea": False,
                   "link": "https://ls11-www.cs.uni-dortmund.de/people/" +
                   "morris/graphkerneldatasets/IMDB-MULTI.zip"},
    "Letter-high": {"nl": False, "el": False, "na": True, "ea": False,
                    "link": "https://ls11-www.cs.uni-dortmund.de/people/mor" +
                    "ris/graphkerneldatasets/Letter-high.zip"},
    "Letter-low": {"nl": False, "el": False, "na": True, "ea": False,
                   "link": "https://ls11-www.cs.uni-dortmund.de/people/mor" +
                   "ris/graphkerneldatasets/Letter-low.zip"},
    "Letter-med": {"nl": False, "el": False, "na": True, "ea": False,
                   "link": "https://ls11-www.cs.uni-dortmund.de/people/mor" +
                   "ris/graphkerneldatasets/Letter-med.zip"},
    "Mutagenicity": {"nl": True, "el": True, "na": False, "ea": False,
                     "link": "https://ls11-www.cs.uni-dortmund.de/peo" +
                     "ple/morris/graphkerneldatasets/Mutagenicity.zip"},
    "MSRC_9": {"nl": True, "el": False, "na": False, "ea": False,
               "link": "https://ls11-www.cs.uni-dortmund.de/people/morris/" +
               "graphkerneldatasets/MSRC_9.zip"},
    "MSRC_21": {"nl": True, "el": False, "na": False, "ea": False,
                "link": "https://ls11-www.cs.uni-dortmund.de/people/mor" +
                "ris/graphkerneldatasets/MSRC_21.zip"},
    "MSRC_21C": {"nl": True, "el": False, "na": False, "ea": False,
                 "link": "https://ls11-www.cs.uni-dortmund.de/people/mor" +
                 "ris/graphkerneldatasets/MSRC_21C.zip"},
    "MUTAG": {"nl": True, "el": True, "na": False, "ea": False,
              "link": "https://ls11-www.cs.uni-dortmund.de/people/mor" +
              "ris/graphkerneldatasets/MUTAG.zip"},
    "NCI1": {"nl": True, "el": False, "na": False, "ea": False,
             "link": "https://ls11-www.cs.uni-dortmund.de/people/mor" +
             "ris/graphkerneldatasets/NCI1.zip"},
    "NCI109": {"nl": True, "el": False, "na": False, "ea": False,
               "link": "https://ls11-www.cs.uni-dortmund.de/people/mor" +
               "ris/graphkerneldatasets/NCI109.zip"},
    "PTC_FM": {"nl": True, "el": True, "na": False, "ea": False,
               "link": "https://ls11-www.cs.uni-dortmund.de/people/mor" +
               "ris/graphkerneldatasets/PTC_FM.zip"},
    "PTC_FR": {"nl": True, "el": True, "na": False, "ea": False,
               "link": "https://ls11-www.cs.uni-dortmund.de/people/mor" +
               "ris/graphkerneldatasets/PTC_FR.zip"},
    "PTC_MM": {"nl": True, "el": True, "na": False, "ea": False,
               "link": "https://ls11-www.cs.uni-dortmund.de/people/mor" +
               "ris/graphkerneldatasets/PTC_MM.zip"},
    "PTC_MR": {"nl": True, "el": True, "na": False, "ea": False,
               "link": "https://ls11-www.cs.uni-dortmund.de/people/mor" +
               "ris/graphkerneldatasets/PTC_MR.zip"},
    "PROTEINS": {"nl": True, "el": False, "na": True, "ea": False,
                 "link": "https://ls11-www.cs.uni-dortmund.de/people/mor" +
                 "ris/graphkerneldatasets/PROTEINS.zip"},
    "PROTEINS_full": {"nl": True, "el": False, "na": True, "ea": False,
                      "link": "https://ls11-www.cs.uni-dortmund.de/people" +
                      "/morris/graphkerneldatasets/PROTEINS_full.zip"},
    "REDDIT-BINARY": {"nl": False, "el": False, "na": False, "ea": False,
                      "link": "https://ls11-www.cs.uni-dortmund.de/people" +
                      "/morris/graphkerneldatasets/REDDIT-BINARY.zip"},
    "REDDIT-MULTI-5K": {"nl": False, "el": False, "na": False, "ea": False,
                        "link": "https://ls11-www.cs.uni-dortmund.de/people/" +
                        "morris/graphkerneldatasets/REDDIT-MULTI-5K.zip"},
    "REDDIT-MULTI-12K": {"nl": False, "el": False, "na": False, "ea": False,
                         "link": "https://ls11-www.cs.uni-dortmund.de/peop" +
                         "le/morris/graphkerneldatasets/REDDIT-MULTI-12K.zip"},
    "SYNTHETIC": {"nl": False, "el": False, "na": True, "ea": False,
                  "link": "https://ls11-www.cs.uni-dortmund.de/people" +
                  "/morris/graphkerneldatasets/SYNTHETIC.zip"},
    "SYNTHETICnew": {"nl": False, "el": False, "na": True, "ea": False,
                     "link": "https://ls11-www.cs.uni-dortmund.de/people/" +
                     "morris/graphkerneldatasets/SYNTHETICnew.zip"},
    "Synthie": {"nl": False, "el": False, "na": True, "ea": False,
                "link": "https://ls11-www.cs.uni-dortmund.de/people/" +
                "morris/graphkerneldatasets/Synthie.zip"},
    "Tox21_AHR": {"nl": True, "el": True, "na": False, "ea": False,
                  "link": "https://ls11-www.cs.uni-dortmund.de/people/" +
                  "morris/graphkerneldatasets/Tox21_AHR.zip"},
    "Tox21_AR": {"nl": True, "el": True, "na": False, "ea": False,
                 "link": "https://ls11-www.cs.uni-dortmund.de/people/morris/" +
                 "graphkerneldatasets/COX2_MD.zip"},
    "Tox21_AR-LBD": {"nl": True, "el": True, "na": False, "ea": False,
                     "link": "https://ls11-www.cs.uni-dortmund.de/people/" +
                     "morris/graphkerneldatasets/Tox21_AR-LBD.zip"},
    "Tox21_ARE": {"nl": True, "el": True, "na": False, "ea": False,
                  "link": "https://ls11-www.cs.uni-dortmund.de/people/" +
                  "morris/graphkerneldatasets/Tox21_ARE.zip"},
    "Tox21_aromatase": {"nl": True, "el": True, "na": False, "ea": False,
                        "link": "https://ls11-www.cs.uni-dortmund.de/people/" +
                        "morris/graphkerneldatasets/Tox21_aromatase.zip"},
    "Tox21_ATAD5": {"nl": True, "el": True, "na": False, "ea": False,
                    "link": "https://ls11-www.cs.uni-dortmund.de/people/" +
                    "morris/graphkerneldatasets/Tox21_ATAD5.zip"},
    "Tox21_ER": {"nl": True, "el": True, "na": False, "ea": False,
                 "link": "https://ls11-www.cs.uni-dortmund.de/people/" +
                 "morris/graphkerneldatasets/Tox21_ER.zip"},
    "Tox21_ER_LBD": {"nl": True, "el": True, "na": False, "ea": False,
                     "link": "https://ls11-www.cs.uni-dortmund.de/people/" +
                     "morris/graphkerneldatasets/Tox21_ER_LBD.zipp"},
    "Tox21_HSE": {"nl": True, "el": True, "na": False, "ea": False,
                  "link": "https://ls11-www.cs.uni-dortmund.de/people/" +
                  "morris/graphkerneldatasets/Tox21_HSE.zip"},
    "Tox21_MMP": {"nl": True, "el": True, "na": False, "ea": False,
                  "link": "https://ls11-www.cs.uni-dortmund.de/people/" +
                  "morris/graphkerneldatasets/Tox21_MMP.zip"},
    "Tox21_p53": {"nl": True, "el": True, "na": False, "ea": False,
                  "link": "https://ls11-www.cs.uni-dortmund.de/people/" +
                  "morris/graphkerneldatasets/Tox21_p53.zip"},
    "Tox21_PPAR-gamma": {"nl": True, "el": True, "na": False, "ea": False,
                         "link": "https://ls11-www.cs.uni-dortmund.de/peop" +
                         "le/morris/graphkerneldatasets/Tox21_PPAR-gamma.zip"}
}

symmetric_dataset = False


def read_data(
        name,
        with_classes=True,
        prefer_attr_nodes=False,
        prefer_attr_edges=False,
        produce_labels_nodes=False,
        as_graphs=False,
        is_symmetric=symmetric_dataset):
    """Create a dataset iterable for GraphKernel.

    Parameters
    ----------
    name : str
        The dataset name.

    with_classes : bool, default=False
        Return an iterable of class labels based on the enumeration.

    produce_labels_nodes : bool, default=False
        Produce labels for nodes if not found.
        Currently this means labeling its node by its degree inside the Graph.
        This operation is applied only if node labels are non existent.

    prefer_attr_nodes : bool, default=False
        If a dataset has both *node* labels and *node* attributes
        set as labels for the graph object for *nodes* the attributes.

    prefer_attr_edges : bool, default=False
        If a dataset has both *edge* labels and *edge* attributes
        set as labels for the graph object for *edge* the attributes.

    as_graphs : bool, default=False
        Return data as a list of Graph Objects.

    is_symmetric : bool, default=False
        Defines if the graph data describe a symmetric graph.

    Returns
    -------
    Gs : iterable
        An iterable of graphs consisting of a dictionary, node
        labels and edge labels for each graph.

    classes : np.array, case_of_appearance=with_classes==True
        An one dimensional array of graph classes aligned with the lines
        of the `Gs` iterable. Useful for classification.

    """
    indicator_path = "./"+str(name)+"/"+str(name)+"_graph_indicator.txt"
    edges_path = "./" + str(name) + "/" + str(name) + "_A.txt"
    node_labels_path = "./" + str(name) + "/" + str(name) + "_node_labels.txt"
    node_attributes_path = "./"+str(name)+"/"+str(name)+"_node_attributes.txt"
    edge_labels_path = "./" + str(name) + "/" + str(name) + "_edge_labels.txt"
    edge_attributes_path = \
        "./" + str(name) + "/" + str(name) + "_edge_attributes.txt"
    graph_classes_path = \
        "./" + str(name) + "/" + str(name) + "_graph_labels.txt"

    # node graph correspondence
    ngc = dict()
    # edge line correspondence
    elc = dict()
    # dictionary that keeps sets of edges
    Graphs = dict()
    # dictionary of labels for nodes
    node_labels = dict()
    # dictionary of labels for edges
    edge_labels = dict()

    # Associate graphs nodes with indexes
    with open(indicator_path, "r") as f:
        for (i, line) in enumerate(f, 1):
            ngc[i] = int(line[:-1])
            if int(line[:-1]) not in Graphs:
                Graphs[int(line[:-1])] = set()
            if int(line[:-1]) not in node_labels:
                node_labels[int(line[:-1])] = dict()
            if int(line[:-1]) not in edge_labels:
                edge_labels[int(line[:-1])] = dict()

    # Extract graph edges
    with open(edges_path, "r") as f:
        for (i, line) in enumerate(f, 1):
            edge = line[:-1].replace(' ', '').split(",")
            elc[i] = (int(edge[0]), int(edge[1]))
            Graphs[ngc[int(edge[0])]].add((int(edge[0]), int(edge[1])))
            if is_symmetric:
                Graphs[ngc[int(edge[1])]].add((int(edge[1]), int(edge[0])))

    # Extract node attributes
    if (prefer_attr_nodes and
        dataset_metadata[name].get(
                "na",
                os.path.exists(node_attributes_path)
                )):
        with open(node_attributes_path, "r") as f:
            for (i, line) in enumerate(f, 1):
                node_labels[ngc[i]][i] = \
                    [float(num) for num in
                     line[:-1].replace(' ', '').split(",")]
    # Extract node labels
    elif dataset_metadata[name].get(
            "nl",
            os.path.exists(node_labels_path)
            ):
        with open(node_labels_path, "r") as f:
            for (i, line) in enumerate(f, 1):
                node_labels[ngc[i]][i] = int(line[:-1])
    elif produce_labels_nodes:
        for i in range(1, len(Graphs)+1):
            node_labels[i] = dict(Counter(s for (s, d) in Graphs[i] if s != d))

    # Extract edge attributes
    if (prefer_attr_edges and
        dataset_metadata[name].get(
            "ea",
            os.path.exists(edge_attributes_path)
            )):
        with open(edge_attributes_path, "r") as f:
            for (i, line) in enumerate(f, 1):
                attrs = [float(num)
                         for num in line[:-1].replace(' ', '').split(",")]
                edge_labels[ngc[elc[i][0]]][elc[i]] = attrs
                if is_symmetric:
                    edge_labels[ngc[elc[i][1]]][(elc[i][1], elc[i][0])] = attrs

    # Extract edge labels
    elif dataset_metadata[name].get(
            "el",
            os.path.exists(edge_labels_path)
            ):
        with open(edge_labels_path, "r") as f:
            for (i, line) in enumerate(f, 1):
                edge_labels[ngc[elc[i][0]]][elc[i]] = int(line[:-1])
                if is_symmetric:
                    edge_labels[ngc[elc[i][1]]][(elc[i][1], elc[i][0])] = \
                        int(line[:-1])

    Gs = list()
    if as_graphs:
        for i in range(1, len(Graphs)+1):
            Gs.append(Graph(Graphs[i], node_labels[i], edge_labels[i]))
    else:
        for i in range(1, len(Graphs)+1):
            Gs.append([Graphs[i], node_labels[i], edge_labels[i]])

    if with_classes:
        classes = []
        with open(graph_classes_path, "r") as f:
            for line in f:
                classes.append(int(line[:-1]))

        classes = np.array(classes, dtype=np.int)
        return Bunch(data=Gs, target=classes)
    else:
        return Bunch(data=Gs)


def _download_zip(url, output_name):
    """Download a file from a requested url and store locally.

    Parameters
    ----------
    url : str
        The url from where the file will be downloaded.

    output_name : str
        The name of the file in the local directory.

    Returns
    -------
    None.

    """
    ctx = ssl.create_default_context(ssl.Purpose.CLIENT_AUTH)
    filename = output_name + ".zip"
    try:
        data_url = urlopen(url, context=ctx)
    except HTTPError as e:
        if e.code == 404:
            e.msg = "Dataset '%s' not found on mldata.org." % output_name
        raise

    # Store Zip File
    try:
        with open(filename, 'w+b') as zip_file:
            copyfileobj(data_url, zip_file)
    except Exception:
        os.remove(filename)
        raise
    data_url.close()


def fetch_dataset(
        name,
        verbose=True,
        data_home=None,
        download_if_missing=True,
        with_classes=True,
        produce_labels_nodes=False,
        prefer_attr_nodes=False,
        prefer_attr_edges=False,
        as_graphs=False):
    """Access a large collection of benchmark datasets from TU Dortmund :cite:`KKMMN2016`.

    For more info visit: :xref:`gd`

    Parameters
    ----------
    name : str
        The name of the dataset (as found in :xref:`gd`).

    verbose : bool, default=True
        Print messages, throughout execution.

    data_home : string, default=None
        Specify another download and cache folder for the datasets.
        By default all grakel data is stored in '~/grakel_data' subfolders.

    download_if_missing : boolean, default=True
        If False, raise a IOError if the data is not locally available instead
        of trying to download the data from the source site.

    with_classes : bool, default=False
        Return an iterable of class labels based on the enumeration.

    produce_labels_nodes : bool, default=False
        Produce labels for nodes if not found.
        Currently this means labeling its node by its degree inside the Graph.
        This operation is applied only if node labels are non existent.

    prefer_attr_nodes : bool, default=False
        If a dataset has both *node* labels and *node* attributes
        set as labels for the graph object for *nodes* the attributes.

    prefer_attr_edges : bool, default=False
        If a dataset has both *edge* labels and *edge* attributes
        set as labels for the graph object for *edge* the attributes.

    as_graphs : bool, default=False
        Return data as a list of Graph Objects.

    Returns
    -------
    graphs : iterable
        Returns an iterable of the produced *valid-graph-format*
        and labels for each node.

    classes : list
        Returns a list of all the classes corresponding to each graph by
        order of input.

    """
    name = str(name)
    if name in dataset_metadata:
        if data_home is None:
            data_home = os.path.join(os.path.expanduser("~"), 'grakel_data')

        exists = os.path.isdir(data_home)
        missing = not os.path.exists(os.path.join(data_home, name + ".zip"))
        cwd = os.getcwd()

        if missing:
            if download_if_missing:
                if not exists:
                    if verbose:
                        print("Initializing folder at", str(data_home))
                    os.makedirs(data_home)
                os.chdir(data_home)
                if verbose:
                    print("Downloading dataset for", name + "..")
                _download_zip(dataset_metadata[name]["link"], name)
            else:
                raise IOError('Dataset ' + name +
                              ' was not found on ' + str(data_home))
        else:
            # move to the general data directory
            os.chdir(data_home)

        with zipfile.ZipFile(str(name) + '.zip', "r") as zip_ref:
            if verbose:
                print("Extracting dataset ", str(name) + "..")
            zip_ref.extractall()

        if verbose:
            print("Parsing dataset ", str(name) + "..")

        data = read_data(name,
                         with_classes=with_classes,
                         prefer_attr_nodes=prefer_attr_nodes,
                         prefer_attr_edges=prefer_attr_edges,
                         produce_labels_nodes=produce_labels_nodes,
                         is_symmetric=symmetric_dataset,
                         as_graphs=as_graphs)
        if verbose:
            print("Parse was succesful..")

        if verbose:
            print("Deleting unzipped dataset files..")
        shutil.rmtree(str(name))

        if verbose:
            print("Going back to the original directory..")
        os.chdir(cwd)

        return data
    else:
        raise ValueError('Dataset: "'+str(name)+'" is currently unsupported.' +
                         '\nSupported datasets come from '
                         'https://ls11-www.cs.tu-dortmund.de/staff/morris/' +
                         'graphkerneldatasets. If your dataset name appears' +
                         ' them send us a pm, to explain you either why we ' +
                         'don\'t support it, or to update our dataset ' +
                         'database.')


def get_dataset_info(dataset_name, default=None):
    """Return the info concerning the existence of a certain dataset.

    Parameters
    ----------
    dataset_name : str
        The name of the dataset.

    default : Object, default=None
        The default return value if the dataset is not found.

    Returns
    -------
    dictionary_get : default or dictionary_entry
        The info of a dataset as a dictionary with
        fields:

            - **nl** : A boolean flag indicating if the
              dataset has node labels.

            - **el** : A boolean flag indicating if the
              dataset has edge labels.

            - **na** : A boolean flag indicating if the
              dataset has node attributes.

            - **ea** : A boolean flag indicating if the
              dataset has edge attributes.

            - **link** : A str corresponding to the
              download link of the dataset.

    """
    return dataset_metadata.get(dataset_name, default)
