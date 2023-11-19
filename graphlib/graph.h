﻿#pragma once

#include <queue>
#include <vector>
#include <unordered_set>
#include <unordered_map>

using std::vector;
using std::unordered_set;
using std::deque;
using std::unordered_map;

namespace graphlib {
;

/// @brief Ребро графа. Описывает из какой вершины выходит, в какую приходит
struct edge_t {
    /// @brief Вершина, откуда ребро выходит
    std::size_t from;
    /// @brief Вершина, куда ребро приходит
    std::size_t to; 
    /// @brief Очевидный конструктор по двум вершинам
    edge_t(std::size_t from, std::size_t to) 
        : from(from)
        , to(to)
    {

    }
    edge_t() = default;
};

/// @brief Структуруа представляющая вершину графа,
/// списки смежности входящих и выходящих ребер
struct vertex_t {
    /// @brief Список индексов ребер, входящих в вершину
    vector<size_t> in;
    /// @brief Список индексов ребер, выходящих из вершины
    vector<size_t> out;
};

/// @brief Граф, заданные списком смежности ребер
class graph_t {
    /// @brief Список смежности ребер графа
    const vector<edge_t> edges;
public:
    graph_t(const vector<edge_t>& edges) // Конструктор с параметром.
        : edges(edges)
    {

    }
    /// @brief Возвращает список инцидентностей для вершин 
    /// с учетом ориентации ребер
    /// Без мемоизации!
    unordered_map<size_t, vertex_t> get_vertices() const {
        unordered_map<size_t, vertex_t> result;
        for (size_t index = 0; index < edges.size(); ++index) {
            size_t from = edges[index].from;
            result[from].out.push_back(index);

            size_t to = edges[index].to;
            result[to].in.push_back(index);
        }

        return result;
    }
   
    /// @brief Топологическая сортировка
    /// @param _start_vertices Начальные вершины, откуда начинать обход
    /// @return (Список пройденных вершин, Список пройденных ребер)
    std::pair<vector<size_t>, vector<size_t>>
        topological_sort() //const vector<size_t>& _start_vertices
    {   
        vector<size_t> vertex_order; 
        vector<size_t> edge_order;
        
        unordered_map<size_t, vertex_t> vertices = get_vertices(); // Получаем вершины графа.

        //
        vector<size_t> _start_vertices_indep;
        for (const auto& vertex : vertices) {
            if (vertex.second.in.empty()) { // Если у вершины нет входящих ребер.
                _start_vertices_indep.push_back(vertex.first); // Добавляем ее индекс в вектор _start_vertices.
            }
        }
        unordered_set<size_t> start_vertices(_start_vertices_indep.begin(), _start_vertices_indep.end());  // Создаем множество стартовых вершин из входящего вектора.        
        std::deque<size_t> bfs_queue(start_vertices.begin(), start_vertices.end()); // Создаем очередь для обхода в глубину, начиная с стартовых вершин.
        //
        
        vector<bool> visited_edges(edges.size(), false); // Инициализируем вектор, чтобы отслеживать посещенные ребра
        // Выполняем обход в ширину
        while (!bfs_queue.empty()) {
            // Получаем вершину из начала очереди.
            size_t vertex = bfs_queue.front();
            bfs_queue.pop_front();          
            vertex_order.push_back(vertex); // Добавляем вершину в вектор vertex_order.           
            const vector<size_t>& out = vertices.at(vertex).out; // Получаем исходящие ребра текущей вершины.
            // Перебираем исходящие ребра.
            for (size_t edge : out) {
                visited_edges[edge] = true; // Помечаем ребро как посещенное.
                edge_order.push_back(edge);  // Добавляем ребро в конец списка edge_order 
                size_t to = edges[edge].to;  // Получаем вершину, в которую ведет данное ребро
                // Проверяем, все ли входящие ребра в целевую вершину были посещены.
                bool to_has_visited_ins = true;  
                for (size_t edge_in : vertices.at(to).in) { 
                    if (visited_edges[edge_in] == false) {
                        to_has_visited_ins = false;
                        break;
                    }
                }
                // Если вершина имеет все посещенные входящие ребра, то добавляем ее в конец очереди bfs_queue
                if (to_has_visited_ins) {
                    bfs_queue.push_back(to);
                }
            }
        }
        // Возвращаем векторы vertex_order и edge_order.
        return std::make_pair(std::move(vertex_order), std::move(edge_order));
    }

};

}

