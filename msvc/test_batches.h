#pragma once
using graphlib::edge_t;
using graphlib::graph_t;

#include <fixed/fixed.h>
#include <pde_solvers/pde_solvers.h>


TEST(Batches, develop)
{
    vector<edge_t> edges{ edge_t(0, 1), edge_t(1, 2), edge_t(1, 3) };

    graph_t g(edges);

    auto [V, E] = g.topological_sort();



}