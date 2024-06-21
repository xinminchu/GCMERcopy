#include <Rcpp.h>

#include "Header/coloring_algorithm.hpp"
#include "Header/dsatur.hpp"
#include "Header/mcs.hpp"
#include "Header/lmxrlf.hpp"
#include "Header/hybrid_dsatur.hpp"
#include "Header/hybrid_lmxrlf.hpp"

using namespace Rcpp;

using std::vector;
using std::map;
using std::string;

using GraphColoring::Dsatur;
using GraphColoring::Mcs;
using GraphColoring::Lmxrlf;
using GraphColoring::HybridDsatur;
using GraphColoring::HybridLmxrlf;
using GraphColoring::GraphColor;

map<string, vector<string> > as_input_graph(ListOf<IntegerVector> adj_list) {
  map<string, vector<string> > input_graph;

  for(ListOf<IntegerVector>::iterator it = adj_list.begin(); it != adj_list.end(); ++it) {
    IntegerVector neighbors = as<IntegerVector>(*it);

    string node = std::to_string(it.index() + 1);
    input_graph[node] = as<vector<string> >(as<CharacterVector>(neighbors));
  }

  return input_graph;
}

IntegerVector as_coloring(GraphColor *graph, int n) {
  map<string,int> coloring = graph->color();

  IntegerVector output(n);
  for(int i = 0; i < n; ++i) {
    string node = std::to_string(i + 1);
    output(i) = coloring[node] + 1;
  }

  return output;
}


// [[Rcpp::export]]
IntegerVector graph_coloring_dsatur(ListOf<IntegerVector> adj_list) {
  GraphColor *graph = new Dsatur(as_input_graph(adj_list));
  return as_coloring(graph, adj_list.size());
}

// [[Rcpp::export]]
IntegerVector graph_coloring_msc(ListOf<IntegerVector> adj_list) {
  GraphColor *graph = new Mcs(as_input_graph(adj_list));
  return as_coloring(graph, adj_list.size());
}

// [[Rcpp::export]]
IntegerVector graph_coloring_lmxrlf(ListOf<IntegerVector> adj_list) {
  GraphColor *graph = new Lmxrlf(as_input_graph(adj_list));
  return as_coloring(graph, adj_list.size());
}

// [[Rcpp::export]]
IntegerVector graph_coloring_hybrid_dsatur_tabucol(ListOf<IntegerVector> adj_list) {
  GraphColor *graph = new HybridDsatur(as_input_graph(adj_list));
  return as_coloring(graph, adj_list.size());
}

// [[Rcpp::export]]
IntegerVector graph_coloring_hybrid_lmxrlf_tabucol(ListOf<IntegerVector> adj_list) {
  GraphColor *graph = new HybridLmxrlf(as_input_graph(adj_list));
  return as_coloring(graph, adj_list.size());
}

// [[Rcpp::export]]
IntegerVector graph_coloring_tabucol(ListOf<IntegerVector> adj_list, int k, int tabu_size = 25, int rep = 100, int nbmax = 1000) {
  GraphColor *graph = new Tabucol(as_input_graph(adj_list), k, tabu_size, rep, nbmax);
  IntegerVector coloring = as_coloring(graph, adj_list.size());

  if(!graph->is_valid()) {
    stop("Graph cannot be colored with " + std::to_string(k) + " colors!");
  }

  return coloring;
}