from astropy.table import Table, vstack
import logging
import re
import sys, yaml
from lsst.daf.butler import Butler
import lsst.geom as geom

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("GetDataFromButler")

operators_list = {
    "<": lambda x, y: x < y,
    ">": lambda x, y: x > y,
    "<=": lambda x, y: x <= y,
    ">=": lambda x, y: x >= y,
    "==": lambda x, y: x == y,
    "=": lambda x, y: x == y,
    "!=": lambda x, y: x != y,
}


def parse_cut(expr):
    # Match operators of 1 or 2 characters at the start of the string
    match = re.match(r"(<=|>=|!=|==|=|<|>)(.*)", expr.strip())
    if match:
        op, value = match.groups()
        return op, value.strip()
    else:
        raise ValueError(f"Cannot parse operator in expression: {expr}")


class GetDataFromButler:
    def __init__(self, repo, collections):
        self.butler = Butler(repo, collections=collections)

    def get_skyMap(self, collections=None, dataId={"skymap": "lsst_cells_v1"}):
        if collections is not None:
            return self.butler.get("skyMap", collections=collections, dataId=dataId)
        else:
            return self.butler.get(
                "skyMap", collections=self.collections, dataId=dataId
            )

    def find_tract_patch(self, ra, dec, skyMap):
        point = geom.SpherePoint(ra * geom.degrees, dec * geom.degrees)
        tp = skyMap.findTractPatchList([point])
        return tp  # returns (tractInfo, patchInfo)

    def get_datasets(self, datasetType, tract):  # condition is a config files
        registry = self.butler.registry
        position_query = "tract={}".format(tract)
        logger.info("Querying datasets with: %s", position_query)
        datasets = list(registry.queryDatasets(datasetType, where=position_query))
        if len(datasets) == 0:
            logger.warning("No datasets found for the given selection criteria.")
        else:
            logger.info("Found %d datasets.", len(datasets))
        return datasets

    def get_data(self, datasets, source_selection_cut=None, quantities=None):
        table_list = []
        for dataset in datasets:
            logger.info("Getting data for dataset: %s", dataset)
            table = self.butler.get(dataset)
            if not isinstance(table, Table):
                table = table.asAstropy()
            if source_selection_cut is not None:
                for key, expr in source_selection_cut.items():
                    if key not in table.colnames:
                        raise KeyError(f"Column '{key}' not found in the data table.")
                    op, value = parse_cut(expr)
                    logger.info("Applying cut: %s %s %s", key, op, value)
                    col_dtype = table[key].dtype
                    if col_dtype.kind in "if":  # integer or float
                        value = float(value)
                    elif col_dtype.kind == "b":  # boolean
                        value = value.lower() in ("true", "1")
                    if op not in operators_list:
                        raise ValueError(f"Unsupported operator: {op} not implemented.")
                    table = table[operators_list[op](table[key], value)]
            if quantities is not None:
                table = table[quantities]
            table_list.append(table)
        return vstack(table_list)


if __name__ == "__main__":
    with open(sys.argv[1]) as fstream:
        cfg = yaml.safe_load(fstream)

    repo = cfg["butler"]["repo"]
    collections = cfg["butler"]["collections"]
    datasetType = cfg["butler"]["datasetType"]
    source_selection_cuts = cfg["query_parameters"]["source_selection_cuts"]
    tracts = cfg["query_parameters"]["tract"]
    quantities = cfg["query_parameters"]["quantities"]
    data_fetcher = GetDataFromButler(repo, collections)
    datasets = []
    for tract in tracts:
        logger.info("Processing tract: %s", tract)
        datasets.extend(data_fetcher.get_datasets(datasetType, tract))
    data = data_fetcher.get_data(
        datasets, source_selection_cut=source_selection_cuts, quantities=quantities
    )
    logger.info(
        "Retrieved data with %d rows and %d columns", len(data), len(data.colnames)
    )
