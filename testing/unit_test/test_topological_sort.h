﻿#pragma once

using graphlib::edge_t;
using graphlib::graph_t;


TEST(TopologicalSort, test)
{
    vector<edge_t> edges{ edge_t(0, 1), edge_t(1, 2), edge_t(1, 3) };

    graph_t g(edges);

    auto [V, E] = g.topological_sort();
}

TEST(TopologicalSort, pump)
{
    vector<edge_t> edges{ edge_t(0, 1), edge_t(1, 3), edge_t(2, 1) };

    graph_t g(edges);

    auto [V, E] = g.topological_sort();
}

TEST(TopologicalSort, pump_two_level)
{
    vector<edge_t> edges{ edge_t(0, 1), edge_t(1, 3), edge_t(2, 1), edge_t(4, 2) };

    graph_t g(edges);

    auto [V, E] = g.topological_sort();
}

TEST(TopologicalSort, random_index)
{
    vector<edge_t> edges{ edge_t(0, 100), edge_t(100, 2), edge_t(100, 3) };

    graph_t g(edges);

    auto [V, E] = g.topological_sort();
}

TEST(TopologicalSort, cycle)
{
    vector<edge_t> edges{ edge_t(0, 1), edge_t(1, 2), edge_t(1, 3), edge_t(2, 3) };

    graph_t g(edges);

    auto [V, E] = g.topological_sort();
}

TEST(TopologicalSort, random_first_edge)
{
    vector<edge_t> edges{ edge_t(0, 1), edge_t(1, 2), edge_t(1, 3)};

    graph_t g(edges);

    auto [V, E] = g.topological_sort();
}

TEST(TopologicalSort, corman)
{
    vector<edge_t> edges{ edge_t(1, 2), edge_t(2, 3), edge_t(3, 4), edge_t(1, 4)};// , edge_t(4, 2)};// , edge_t(5, 3), edge_t(5, 6)};

    graph_t g(edges);

    auto [V, E] = g.topological_sort();
}
