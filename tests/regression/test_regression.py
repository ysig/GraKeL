"""Regression checks against frozen kernel matrices."""
import ast
import json
from pathlib import Path

import grakel
import numpy as np
import pytest


REGRESSION_DIR = Path(__file__).resolve().parent
DATA_DIR = REGRESSION_DIR / "data"
EXPECTED_DIR = REGRESSION_DIR / "expected"


def _read_json(path):
    with path.open("r", encoding="utf-8") as stream:
        return json.load(stream)


def _parse_edge_key(key):
    value = ast.literal_eval(key)
    if not isinstance(value, tuple) or len(value) != 2:
        raise ValueError("invalid edge key: {!r}".format(key))
    return tuple(int(node) for node in value)


def _load_graphs():
    graphs = {}
    for record in _read_json(DATA_DIR / "graphs.json")["graphs"]:
        edge_dictionary = {
            int(node): [int(neighbor) for neighbor in neighbors]
            for node, neighbors in record["edge_dictionary"].items()
        }
        node_labels = {
            int(node): label for node, label in record["node_labels"].items()
        }
        edge_labels = {
            _parse_edge_key(edge): label
            for edge, label in record["edge_labels"].items()
        }
        graphs[record["name"]] = grakel.Graph(
            edge_dictionary,
            node_labels=node_labels,
            edge_labels=edge_labels,
            graph_format="dictionary",
        )
    return graphs


def _build_kernel(case):
    specification = case["kernel"]
    params = dict(specification.get("params", {}))
    base = params.get("base_graph_kernel")
    if isinstance(base, dict):
        params["base_graph_kernel"] = (
            getattr(grakel, base["class"]), dict(base.get("params", {})))
    return getattr(grakel, specification["class"])(**params)


def _run_case(case, graphs, split):
    train = [graphs[name] for name in split["train"]]
    test = [graphs[name] for name in split["test"]]

    kernel = _build_kernel(case)
    K_train = kernel.fit_transform(train)

    kernel = _build_kernel(case)
    kernel.fit(train)
    K_test = kernel.transform(test)

    kernel = _build_kernel(case)
    K_repeat = kernel.fit_transform(train)
    return K_train, K_test, np.allclose(K_train, K_repeat, rtol=1e-9, atol=1e-9)


specifications = _read_json(DATA_DIR / "specifications.json")
graphs = _load_graphs()
expected_manifest = _read_json(EXPECTED_DIR / "manifest.json")


@pytest.mark.parametrize("case", specifications["cases"],
                         ids=[case["id"] for case in specifications["cases"]])
def test_kernel_matches_frozen_output(case):
    choice = expected_manifest["cases"][case["id"]]
    K_train, K_test, repeat_match = _run_case(case, graphs, specifications["split"])
    expected = _read_json(EXPECTED_DIR / (case["id"] + ".json"))

    np.testing.assert_allclose(K_train, expected["K_train"],
                               rtol=choice["rtol"], atol=choice["atol"])
    np.testing.assert_allclose(K_test, expected["K_test"],
                               rtol=choice["rtol"], atol=choice["atol"])
    assert repeat_match
