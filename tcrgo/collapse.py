import numpy as np
import pandas as pd
import igraph as ig
import os
from typing import Optional, Tuple, List, Any
from .bio import hamming_distance
from .log import Log

Array = np.ndarray
DataFrame = pd.core.frame.DataFrame
Graph = ig.Graph
Vertex = ig.Vertex
VertexSeq = ig.VertexSeq
log = Log("root")

def preprocess_dataframe(df: DataFrame) -> DataFrame:
	# Ugly. topVJ should have been two separate reference columns.
	df["topVJ_region"] = df["topVJ_region"].apply(
		lambda x: tuple(x.replace("('", "").replace("')", "").split("', '"))
	)
	df["BCUMI"] = df["BC"] + df["UMI"]
	df = df.dropna(subset=["CDR3_nt"]) # TODO: Why are some CDR3s blank?
	return df
	
def get_neighbors(collapser_bcumi: str, collapser_cdr3: str, 
	collapsees_bcumi: Array, collapsees_cdr3: Array) -> Optional[object]:
	"""
	Find Collapsee BCUMIs with HD <= 1 to collapser BCUMI, return edges and
	CDR3 distances in an object.
	If I could do this over again, I wouldn't have made this an apply method
	and probably would have gone for a 2D collapser-collapsee matrix.
	Still wondering if there's a better approach to dropping None.
	"""
	hd_vectorized = np.vectorize(hamming_distance)
	# Compare collapser BCUMI to each collapsee BCUMI, drop those with HD > 1
	collapser_bcumi_repeated = np.repeat(collapser_bcumi, len(collapsees_bcumi))
	distances_bcumi = hd_vectorized(collapser_bcumi_repeated, collapsees_bcumi)
	neighbor_bcumis, neighbor_cdr3s = \
		np.where(distances_bcumi == 1 , (collapsees_bcumi, collapsees_cdr3), None)
	neighbor_bcumis = neighbor_bcumis[neighbor_bcumis != None]
	neighbor_cdr3s = neighbor_cdr3s[neighbor_cdr3s != None]

	edge_info = None
	if neighbor_bcumis.size > 0:
		# Retain BCUMIs which share same CDR3 length with collapser's CDR3
		neighbor_cdr3s, neighbor_bcumis = np.where(
			np.frompyfunc(len, 1, 1)(neighbor_cdr3s) == len(collapser_cdr3), 
			(neighbor_cdr3s, neighbor_bcumis), None
		)
		neighbor_bcumis = neighbor_bcumis[neighbor_bcumis != None]
		neighbor_cdr3s = neighbor_cdr3s[neighbor_cdr3s != None]
		
		if neighbor_bcumis.size > 0:
			# Compare collapser CDR3 to each collapsee CDR3, drop those with HD > 1
			collapser_cdr3_repeated = np.repeat(collapser_cdr3, len(neighbor_cdr3s))
			distances_cdr3 = hd_vectorized(collapser_cdr3_repeated, neighbor_cdr3s)
			neighbor_bcumis = np.where(distances_cdr3 <= 1, neighbor_bcumis, None)
			neighbor_bcumis = neighbor_bcumis[neighbor_bcumis != None]

			# Of those retained, return edges and CDR3 distances (in obj to avoid error)
			neighbors_distances = distances_cdr3[distances_cdr3 <= 1]
			collapser_bcumi_repeated = collapser_bcumi_repeated[:len(neighbor_bcumis)]
			edges = np.column_stack((neighbor_bcumis, collapser_bcumi_repeated))
			edge_info = {"edges": edges, "distances": neighbors_distances.astype(int)}
	log.vsep("-")
	log.verbose(f"{edge_info}")
	return edge_info

def populate_vertices(df_ref: DataFrame):
	graph = ig.Graph(directed=True)
	graph.add_vertices(df_ref["BCUMI"].to_list())
	for col in df_ref.columns:
		graph.vs[col] = df_ref[col].to_list()
	graph.vs["index"] = df_ref.index.to_list()
	return graph

def connect_vertices(graph: Graph, collapsers_edge_info: List[object]):
	edges = np.empty((0,2), dtype='<U32')
	distances= np.array([])
	for edge_info in collapsers_edge_info:
		edges = np.append(edges, edge_info["edges"], axis=0)
		distances = np.append(distances, edge_info["distances"])
	graph.add_edges(edges)
	graph.es["distances"] = distances
	return graph

def plot(graph: Graph, collapsers: Array, filepath: str):
	color = [
		"salmon" if v["BCUMI"] in collapsers else "gold" \
			for v in graph.vs
	]
	ig.plot(
		graph,
		layout=graph.layout("auto"),
		target=filepath,
		vertex_label=graph.vs["nReads"],
		vertex_label_size=8,
		edge_label=graph.es["distances"],
		edge_label_size=8,
		vertex_size=10 * np.log10(graph.vs["nReads"]) + 10,
		vertex_color=color
	)

def identify_vertices(cluster: Graph) -> Tuple[Vertex, VertexSeq]:
	collapser = cluster.vs.select(nReads = max(cluster.vs["nReads"]))
	if len(collapser) >= 2:
		collapser = collapser.select(_degree = collapser.maxdegree())
	collapser = collapser[0] # Arbitrarily select first if stil ties
	collapsees = cluster.vs.select(index_ne = collapser["index"])
	return collapser, collapsees

def collapse_information(collapser: Vertex, collapsees: VertexSeq) -> Vertex:
	# I didn't bother with the TR(A|B)(V|J|C) nReads columns.
	for attr in ("nReads", "topVJ_nReads", "CDR3_nReads"):
		collapser[attr] += sum(collapsees[attr])				
	collapser["topVJ_freq"] = collapser["topVJ_nReads"] / collapser["nReads"]
	collapser["CDR3_freq"] = collapser["CDR3_nReads"] / collapser["nReads"]
	return collapser

def build_row(vertex: Vertex, collapsees: List[str]=None) -> List[Any]:
	return [*vertex.attributes().values()] + [collapsees]

def collapse(df: DataFrame, threshold: int=10, do_plot: 
	bool=True, outpath: str='') -> DataFrame:
	""" """
	df = preprocess_dataframe(df)
	entries = list()
	for ref, df_ref in df.groupby("topVJ_region"):
		log.verbose(f"REF: {ref}")
		graph = populate_vertices(df_ref)
		df_collapsers = df_ref.loc[df_ref.nReads >= threshold]
		df_collapsees = df_ref.loc[df_ref.nReads < threshold]
		if len(df_collapsees.index) > 0:
			# index BCUMI=-1, index CDR3_nt=8
			collapsers_edge_info = df_collapsers.apply(
				lambda row: get_neighbors(
					row[-1], row[8], 
					df_collapsees["BCUMI"].to_numpy(), 
					df_collapsees["CDR3_nt"].to_numpy()),
				axis=1, raw=True).dropna()
			if len(collapsers_edge_info.index) > 0:
				graph = connect_vertices(graph, collapsers_edge_info)
				if graph.es:
					vertices_isolated = graph.vs.select(_degree=0)
					for vertex in vertices_isolated:
						entries.append([*vertex.attributes().values()] + [None])
					graph.delete_vertices(vertices_isolated)
					if do_plot:
						filepath = os.path.join(outpath, f"{ref[0]}_{ref[1]}.png")
						plot(graph, df_collapsers["BCUMI"].values, filepath)
					for cluster in graph.clusters(mode="WEAK").subgraphs():
						collapser, collapsees = identify_vertices(cluster)
						collapser = collapse_information(collapser, collapsees)
						entries.append([*collapser.attributes().values()] + [collapsees["BCUMI"]])
					continue # Everything has been written, skip to next. 
		for vertex in graph.vs:
			entries.append(build_row(vertex))
		log.vsep('+')
	columns = ["name"] + list(df.columns) + ["index", "Collapsed"]
	df = pd.DataFrame(entries, columns=columns)
	df = df.drop(["name", "index", "BCUMI"], axis=1)
	return df