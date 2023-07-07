

//#include <graphlib/graph.h>
//
//using graphlib::edge_t;
//using graphlib::graph_t;
//
//
//int main()
//{
//    vector<edge_t> edges{ edge_t(100, 1), edge_t(1, 2), edge_t(1, 3) };
//
//
//    graph_t g(edges);
//
//    auto [V, E] = g.topological_sort({100});
//
//
//}
//
//

#define GTEST_BREAK_ON_FAILURE 1
#define GTEST_CATCH_EXCEPTIONS 0
#define GTEST_HAS_SEH 0
#define _VARIADIC_MAX 10 /* for gtest */
#include "gtest/gtest.h"

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
#ifdef _WIN32
    std::wcout.imbue(std::locale("rus_rus.866"));
#endif
    int res = RUN_ALL_TESTS();
    return res;
}
