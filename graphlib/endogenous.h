#pragma once

#include <fixed/fixed.h>
#include <pde_solvers/pde_solvers.h>

namespace graphlib {
;


/// @brief Хранит все структуры данных для расчета эндогенных свойств сети
/// Пока заточен под расчет только плотности, явно будет дорабатываться
/// Включает в себя топологию графа и данные по каждой трубе
/// Данные по трубе включат в себя параметры труб, их расходы, 
/// экземпляры модели транспортного уравнения, буферы, в которых считаются ДУЧП
struct task_network_propagation_data_t {
    // Профиль переменных
    typedef quickest_ultimate_fv_solver_traits<1>::var_layer_data target_var_t;
    typedef quickest_ultimate_fv_solver_traits<1>::specific_layer specific_data_t;

    // Слой: переменных Vars + сколько угодно служебных Specific
    typedef composite_layer_t<target_var_t, specific_data_t> layer_t;

    /// @brief Трубы
    vector<pipe_properties_t> pipes;
    /// @brief Профили расходов по каждой трубе
    vector<vector<double>> Q;
    /// @brief Модели адвекции
    vector<PipeQAdvection> models;
    /// @brief Буферы для расчета адвекции
    vector<ring_buffer_t<layer_t>> buffers;
    /// @brief Топология графа
    vector<edge_t> edges;
    
    /// @brief Инициализация буферов по начальным условиям
    /// @param rho_initial 
    void init_buffers(double rho_initial)
    {
        for (auto& buffer : buffers) {
            for (double& rho : buffer.previous().vars.cell_double[0]) {
                rho = rho_initial;
            }
        }
    }

    /// @brief Делает advance всем буферам всех труб
    void advance_buffers()
    {
        for (auto& buffer : buffers) {
            buffer.advance(+1);
        }
    }

    /// @brief Генерирует простую сеть с разделением
    /// В будущем надо унести отсюда
    static task_network_propagation_data_t default_data()
    {
        // Упрощенное задание трубы - 50км, с шагом разбиения для расчтной сетки 1км, диаметром 700мм
        simple_pipe_properties simple_pipe_1;
        simple_pipe_properties simple_pipe_2;
        simple_pipe_properties simple_pipe_3;
        simple_pipe_1.length = 100e3;
        simple_pipe_2.length = 168e3;
        simple_pipe_3.length = 263e3;
        //simple_pipe.length = 700e3; // тест трубы 700км
        simple_pipe_1.diameter = 1;
        simple_pipe_2.diameter = 1;
        simple_pipe_3.diameter = 1;
        //simple_pipe.diameter = 0.514; // тест трубы 700км
        simple_pipe_1.dx = 1000;
        simple_pipe_2.dx = 1000;
        simple_pipe_3.dx = 1000;
        //simple_pipe.dx = 100; // тест трубы 700км
        auto pipe_1 = pipe_properties_t::build_simple_pipe(simple_pipe_1);
        auto pipe_2 = pipe_properties_t::build_simple_pipe(simple_pipe_2);
        auto pipe_3 = pipe_properties_t::build_simple_pipe(simple_pipe_3);

        task_network_propagation_data_t result;

        result.pipes = vector <pipe_properties_t>{ pipe_1, pipe_2, pipe_3 };

        vector<double> vol_flows{ 0, 2, 2 };
        for (size_t index = 0; index < result.pipes.size(); ++index) {
            result.Q.emplace_back(
                vector<double>(result.pipes[index].profile.getPointCount(), vol_flows[index]));
        }

        for (size_t index = 0; index < result.pipes.size(); ++index) {
            result.models.emplace_back(result.pipes[index], result.Q[index]);
        }

        for (size_t index = 0; index < result.pipes.size(); ++index) {
            result.buffers.emplace_back(2, result.pipes[index].profile.getPointCount());
        }

        result.edges = vector<edge_t>{
            edge_t(1, 2), edge_t(0, 1), edge_t(1, 3) }; //!! Изменение порядка рёбер влияет на порядок начальных условий

        return result;
    }
};


/// @brief Солвер сети для движения партий
class net_quickest_ultimate_solver
{
protected:
    task_network_propagation_data_t& data;
public:
    net_quickest_ultimate_solver(task_network_propagation_data_t& data)
        : data(data)
    {

    }
    /// @brief Делает расчет шага по заданным граничным условиям сети
    /// В узлах смешения берет взвешенное среднее по расходам
    /// @param dt 
    /// @param boundaries Граничные условия сети. 
    /// Ключ - индекс вершины, значение - величина граничного условия
    /// @return Временные ряды эндогенных свойств в вершинах
    /// Ключ - индекс вершины, значение - временной ряд ?? Почему ряд?
    std::map<size_t, vector<double>> step(
        double dt, const std::map<size_t, double>& boundaries)
    {
        graph_t g(data.edges);
        auto [V, E] = g.topological_sort();
        auto vertices = g.get_vertices();

        std::map<size_t, vector<double>> vertices_density;

        for (size_t vertex : V) {

            const auto& E = vertices[vertex];

            double density_vertex;
            if (E.in.empty()) { // Если у текущей вершины нет входных рёбер
                density_vertex = boundaries.at(vertex); // то должны быть начальные условия
            }
            else {
                double Summ_Rho_Q = 0;
                double Summ_Q = 0;
                for (const auto& edg : E.in) { // если их 1 и больше, то ищем смесь
                    const auto& current = data.buffers[edg].current().vars.cell_double[0];
                    double rho_from_edge = current.back();
                    double q_edge = data.Q[edg].back();
                    Summ_Q += q_edge;
                    Summ_Rho_Q += q_edge * rho_from_edge;
                }
                density_vertex = Summ_Rho_Q / Summ_Q;
            }

            vertices_density[vertex].push_back(density_vertex);

            for (const auto& edge : E.out) {

                double density_vertex_out = -1; // фиктивное краевое условие на выходе ребра

                quickest_ultimate_fv_solver solver(data.models[edge], data.buffers[edge]);
                solver.step(dt, density_vertex, density_vertex_out);

            }
        }
        return vertices_density;
    }



};


}