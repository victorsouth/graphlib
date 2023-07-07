

#include <graphlib/graph.h>

using graphlib::edge_t;
using graphlib::graph_t;


int main()
{
    vector<edge_t> edges{ edge_t(100, 1), edge_t(1, 2), edge_t(1, 3) };


    graph_t g(edges);

    auto [V, E] = g.topological_sort({100});


}