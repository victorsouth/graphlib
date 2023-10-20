#pragma once
using graphlib::edge_t;
using graphlib::graph_t;

#include <fixed/fixed.h>
#include <pde_solvers/pde_solvers.h>

using namespace std;


/// @brief Тесты (ТУ) для солвера quickest_ultimate_fv_solver
class QUICKEST_ULTIMATE_TU : public ::testing::Test {
protected:
    // Профиль переменных
    typedef quickest_ultimate_fv_solver_traits<1>::var_layer_data target_var_t;
    typedef quickest_ultimate_fv_solver_traits<1>::specific_layer specific_data_t;

    // Слой: переменных Vars + сколько угодно служебных Specific
    typedef composite_layer_t<target_var_t, specific_data_t> layer_t;
protected:
    ///// @brief Параметры трубы
    //PipeProperties pipe_1;
    //PipeProperties pipe_2;
    //PipeProperties pipe_3;
    ///// @brief Профиль расхода
    //vector<double> Q_1;
    //vector<double> Q_2;
    //vector<double> Q_3;
    //std::unique_ptr<PipeQAdvection> advection_model_1;
    //std::unique_ptr<PipeQAdvection> advection_model_2;
    //std::unique_ptr<PipeQAdvection> advection_model_3;
    //std::unique_ptr<custom_buffer_t<layer_t>> buffer_1;
    //std::unique_ptr<custom_buffer_t<layer_t>> buffer_2;
    //std::unique_ptr<custom_buffer_t<layer_t>> buffer_3;
    vector <PipeProperties> pipes;
    vector<vector<double>> Q;
    vector <PipeQAdvection> models;
    vector <custom_buffer_t<layer_t>> buffers;
    vector <layer_t> prevs;


protected:

    /// @brief Подготовка к расчету для семейства тестов
    virtual void SetUp() override {
        // Упрощенное задание трубы - 50км, с шагом разбиения для расчтной сетки 1км, диаметром 700мм

        // Упрощенное задание трубы - 50км, с шагом разбиения для расчтной сетки 1км, диаметром 700мм
        simple_pipe_properties simple_pipe_1;
        simple_pipe_properties simple_pipe_2;
        simple_pipe_properties simple_pipe_3;
        simple_pipe_1.length = 168e3;
        simple_pipe_2.length = 100e3;
        simple_pipe_3.length = 263e3;
        //simple_pipe.length = 700e3; // тест трубы 700км
        simple_pipe_1.diameter = 0.1;
        simple_pipe_2.diameter = 0.1;
        simple_pipe_3.diameter = 0.1;
        //simple_pipe.diameter = 0.514; // тест трубы 700км
        simple_pipe_1.dx = 1000;
        simple_pipe_2.dx = 1000;
        simple_pipe_3.dx = 1000;
        //simple_pipe.dx = 100; // тест трубы 700км
        auto pipe_1 = PipeProperties::build_simple_pipe(simple_pipe_1);
        auto pipe_2 = PipeProperties::build_simple_pipe(simple_pipe_2);
        auto pipe_3 = PipeProperties::build_simple_pipe(simple_pipe_3);

        pipes = vector <PipeProperties>{ pipe_1, pipe_2, pipe_3 };

        vector<double> vol_flows{0.5, 0.2, 0.3};
        for (size_t index = 0; index < pipes.size(); ++index) {
            Q.emplace_back(vector<double>(pipes[index].profile.getPointCount(),
                vol_flows[index]));
        }

        for (size_t index = 0; index < pipes.size(); ++index) {
            models.emplace_back(pipes[index], Q[index]);
        }

        for (size_t index = 0; index < pipes.size(); ++index) {
            buffers.emplace_back(2, pipes[index].profile.getPointCount());
        }

        //for (size_t index = 0; index < pipes.size(); ++index) {
        //    prevs.emplace_back(buffers[index].previous());
        //}
        //       prev.vars.cell_double[0] = vector<double>(prev.vars.cell_double[0].size(), 850);
    }
};

/// @brief Проверка QUICKEST-ULTIMATE, проверка изменения плотности после тройника смешения
TEST_F(QUICKEST_ULTIMATE_TU, MixDensity) {

    string path = prepare_test_folder();

    std::map<size_t, double> boundaries {
        { 0, 870},
        { 2, -1 },
        { 3, -1 }
    };
    double rho_initial = 850;
    for (auto& buffer : buffers) {
        for (double& rho : buffer.previous().vars.cell_double[0]) {
            rho = rho_initial;
        }
    }

    double T = 300000; // период моделирования
    //double T = 800000; // период моделирования (тест трубы 700км)

    vector<edge_t> edges{ edge_t(1, 2), edge_t(0, 1), edge_t(1, 3) };
    graph_t g(edges);
    auto [V, E] = g.topological_sort();
    auto vertices = g.get_vertices();

    const auto& x = models[0].get_grid();
    double dx = x[1] - x[0];
    double v = models[0].getEquationsCoeffs(0, 0);
    double dt_ideal = abs(dx / v);
    double Cr = 1;

    //std::unique_ptr<custom_buffer_t<layer_t>> buffer;
    //std::unique_ptr<PipeQAdvection> advection_model;

    double t = 0; // текущее время
    //double dt = 60; // 1 минута
    double dt = Cr * dt_ideal; // время в долях от Куранта
    size_t N = static_cast<int>(T / dt);


    std::map<size_t, vector<double>> vertices_density;

    //for (size_t vertex = 0; vertex < n_vertex; ++vertex) {
    for (size_t index = 0; index < N; ++index) {
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
                    const auto& current = buffers[edg].current().vars.cell_double[0];
                    double rho_from_edge = current.back();
                    double q_edge = Q[edg].back(); 
                    Summ_Q += q_edge;
                    Summ_Rho_Q += q_edge * rho_from_edge;
                }
                density_vertex = Summ_Rho_Q / Summ_Q;
            }

            vertices_density[vertex].push_back(density_vertex);

            for (const auto& edge : E.out) {
                // для vertex определить ребра, которые в него входят, из них взять выходные плотности, смешать        
                // для vertex определить ребра, которые ИЗ него ВЫХОДЯТ, 
                // то что намешал использовать как граничное условие расчета для QUICKEST
                
                //std::stringstream filename;
                //filename << path << "output " << edge << " Cr=" << Cr << ".csv";
                //std::ofstream output(filename.str());

            
                //if (index == 0) {
                //    prevs[edge] = buffers[edge].previous();
                //    prevs[edge].vars.print(t, output);
                //}

                //for (const auto& vert : vertices) {
                //    if (vert.first == vertex) { // Поиск текущей вершины.
                //        if (vert.second.in.empty()) { // Если у текущей вершины нет входных рёбер
                //            Rho_smesi = rho_in_1; // то должны быть начальные условия
                //        }
                //        else {
                //            double Summ_Rho_Q = 0;
                //            double Summ_Q = 0;
                //            for (const auto& j : vert.second.in) { // Проход по всем входящим ребрам
                //                Summ_Rho_Q = Summ_Rho_Q + (nexts[j].vars.cell_double[0][nexts[j].vars.cell_double[0].size() - 1] * Q[j][Q[j].size() - 1]);
                //                Summ_Q = Summ_Q + (Q[j][Q[j].size() - 1]);
                //            }
                //            Rho_smesi = Summ_Rho_Q / Summ_Q;
                //        }
                //    }
                //}

                double density_vertex_out = -1; // фиктивное краевое условие на выходе ребра

                quickest_ultimate_fv_solver solver(models[edge], buffers[edge]);
                solver.step(dt, density_vertex, density_vertex_out);

                //output.flush();
                //output.close();

            }
        }
        t += dt;

        for (auto& buffer : buffers) {
            buffer.advance(+1);
        }
    }

}
