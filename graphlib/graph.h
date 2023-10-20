#pragma once

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

// ���������, �������������� ����� �����.
struct edge_t {
    std::size_t from;
    std::size_t to; 
    edge_t(std::size_t from, std::size_t to) // ����������� � �����������.
        : from(from)
        , to(to)
    {

    }
    edge_t() = default;
};
// ���������, �������������� ������� �����.
struct vertex_t {
    vector<size_t> in;
    vector<size_t> out;
};

class graph_t 
{
    const vector<edge_t> edges; // ������ ����� �����.
public:
    graph_t(const vector<edge_t>& edges) // ����������� � ����������.
        : edges(edges)
    {

    }
    // ����� ��� ��������� ������ �����.
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
   
    /// @brief �������������� ����������
    /// @param _start_vertices ��������� �������, ������ �������� �����
    /// @return (������ ���������� ������, ������ ���������� �����)
    std::pair<vector<size_t>, vector<size_t>>
        topological_sort() //const vector<size_t>& _start_vertices
    {   

        //unordered_set<size_t> start_vertices(_start_vertices.begin(), _start_vertices.end());  // ������� ��������� ��������� ������ �� ��������� �������.        
        //std::deque<size_t> bfs_queue(start_vertices.begin(), start_vertices.end()); // ������� ������� ��� ������ � �������, ������� � ��������� ������.
        // �������������� ������� ��� ���������� ������� ������ � ����� � �������������� ����������.
        vector<size_t> vertex_order; 
        vector<size_t> edge_order;
        
        unordered_map<size_t, vertex_t> vertices = get_vertices(); // �������� ������� �����.

        //
        vector<size_t> _start_vertices_indep;
        for (const auto& vertex : vertices) {
            if (vertex.second.in.empty()) { // ���� � ������� ��� �������� �����.
                _start_vertices_indep.push_back(vertex.first); // ��������� �� ������ � ������ _start_vertices.
            }
        }
        unordered_set<size_t> start_vertices(_start_vertices_indep.begin(), _start_vertices_indep.end());  // ������� ��������� ��������� ������ �� ��������� �������.        
        std::deque<size_t> bfs_queue(start_vertices.begin(), start_vertices.end()); // ������� ������� ��� ������ � �������, ������� � ��������� ������.
        //
        
        vector<bool> visited_edges(edges.size(), false); // �������������� ������, ����� ����������� ���������� �����
        // ��������� ����� � ������
        while (!bfs_queue.empty()) {
            // �������� ������� �� ������ �������.
            size_t vertex = bfs_queue.front();
            bfs_queue.pop_front();          
            vertex_order.push_back(vertex); // ��������� ������� � ������ vertex_order.           
            const vector<size_t>& out = vertices.at(vertex).out; // �������� ��������� ����� ������� �������.
            // ���������� ��������� �����.
            for (size_t edge : out) {
                visited_edges[edge] = true; // �������� ����� ��� ����������.
                edge_order.push_back(edge);  // ��������� ����� � ����� ������ edge_order 
                size_t to = edges[edge].to;  // �������� �������, � ������� ����� ������ �����
                // ���������, ��� �� �������� ����� � ������� ������� ���� ��������.
                bool to_has_visited_ins = true;  
                for (size_t edge_in : vertices.at(to).in) { 
                    if (visited_edges[edge_in] == false) {
                        to_has_visited_ins = false;
                        break;
                    }
                }
                // ���� ������� ����� ��� ���������� �������� �����, �� ��������� �� � ����� ������� bfs_queue
                if (to_has_visited_ins) {
                    bfs_queue.push_back(to);
                }
            }
        }
        // ���������� ������� vertex_order � edge_order.
        return std::make_pair(std::move(vertex_order), std::move(edge_order));
    }

};

}

