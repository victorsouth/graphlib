﻿#pragma once


using graphlib::edge_t;
using graphlib::graph_t;
using graphlib::task_network_propagation_data_t;
using graphlib::net_quickest_ultimate_solver;
using namespace std;




struct timeseries_data {
    vector<double> density_in;
    vector<double> density_out;
    vector<double> volflow_in;
    vector<double> volflow_out;
    timeseries_data(const char* filename = "sikn_c.csv")
    {
        std::ifstream file(filename); // Открываем файл для чтения
        std::string line;

        std::vector<std::vector<double>> data;

        while (getline(file, line)) {
            std::stringstream ss(line);
            std::string value;
            size_t k = 0;
            while (getline(ss, value, ',')) {
                if (k == data.size()) { // Если вектора для столбца ещё нет, создаём его
                    data.emplace_back();
                }
                data[k].push_back(std::stod(value)); // Записываем значение в столбец
                k++; //1825
            }
        }

        file.close();

        density_in = data[1];
        density_out = data[0];
        volflow_in = data[3];
        volflow_out = data[2];

        for (double& q : volflow_in) {
            q /= 3600;
        }
        for (double& q : volflow_out) {
            q /= 3600;
        }

    }
    /// @brief Возвращает Dt, гарантирующий, что по всем трубам системы 
    /// не будет Cr > 1
    /// Пока что очень частное решение на основе одной трубы
    /// @param data 
    /// @return 
    double calc_ideal_dt(vector <PipeQAdvection>& models) const {
        //double v = models[1].getEquationsCoeffs(0, 0); // !!! В РУЧНУЮ ДЛЯ 1-ого РЕБРА
        auto it = std::max_element(volflow_in.begin(), volflow_in.end(), 
            [](double q1, double q2) {
                return abs(q1) < abs(q2);
            });

        double Q_max = abs(*it);

        double v_max = Q_max / models[1].get_pipe().wall.getArea();
        const auto& x = models[1].get_grid();
        double dx = x[1] - x[0];

        double dt_ideal = dx / v_max;
        return dt_ideal;
    }
};


/// @brief Тесты (ТУ) для солвера quickest_ultimate_fv_solver
class QUICKEST_ULTIMATE_TU : public ::testing::Test {
protected:
    task_network_propagation_data_t net_data;
protected:
    /// @brief Подготовка к расчету для семейства тестов
    virtual void SetUp() override {
        net_data = task_network_propagation_data_t::default_data();

        double rho_initial = 850;
        net_data.init_buffers(rho_initial);
    }


};

/// @brief Тесты (ТУ) для солвера quickest_ultimate_fv_solver RandomEdges
class QUICKEST_ULTIMATE_TU_RE : public ::testing::Test {
protected:
    // Профиль переменных
    typedef quickest_ultimate_fv_solver_traits<1>::var_layer_data target_var_t;
    typedef quickest_ultimate_fv_solver_traits<1>::specific_layer specific_data_t;

    // Слой: переменных Vars + сколько угодно служебных Specific
    typedef composite_layer_t<target_var_t, specific_data_t> layer_t;
protected:
    ///// @brief Параметры трубы

    vector <pipe_properties_t> pipes;
    vector<vector<double>> Q;
    vector <PipeQAdvection> models;
    vector <ring_buffer_t<layer_t>> buffers;
    vector <layer_t> prevs;


protected:

    /// @brief Подготовка к расчету для семейства тестов
    virtual void SetUp() override {

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

        pipes = vector <pipe_properties_t>{ pipe_3, pipe_2, pipe_1 };

        vector<double> vol_flows{2, 4, 2};
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

    timeseries_data data;
    double dt_ideal = data.calc_ideal_dt(net_data.models);
    double dt = 300; // время по реальным данным
    double Cr = dt_ideal/dt;

    double T = 300000; // период моделирования
    //double T = 800000; // период моделирования (тест трубы 700км)
    size_t N = static_cast<int>(T / dt);
    double t = 0; // текущее время

    std::stringstream filename;
    filename << path << "Rho" << ".csv";
    std::ofstream output(filename.str());

    for (size_t index = 0; index < N; ++index) {        
        // Учет краевых условий и расходов на новом шаге
        for (size_t i = 1; i < net_data.pipes.size(); ++i) { // Костыль
            double q_pipe = i == 1
                ? data.volflow_in[index]
                : data.volflow_out[index];
            net_data.Q[i] = vector<double>(net_data.pipes[i].profile.getPointCount(), q_pipe);
        }

        auto& next = net_data.buffers[2].current();
        next.vars.print(t, output);

        boundaries[0] = data.density_in[index];

        net_quickest_ultimate_solver solver(net_data);
        std::map<size_t, vector<double>> vertices_density
            = solver.step(dt, boundaries);
        net_data.advance_buffers();
        t += dt;
    }
    output.flush();
    output.close();

}

/// @brief Проверка QUICKEST-ULTIMATE, проверка задания ребер в случайно порядке
TEST_F(QUICKEST_ULTIMATE_TU_RE, RandomEdges) {

    string path = prepare_test_folder();

    timeseries_data data;

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

    vector<edge_t> edges{ edge_t(1, 2), edge_t(0, 1), edge_t(1, 3) }; //!! Изменение порядка рёбер влияет на порядок начальных условий
    graph_t g(edges);
    auto [V, E] = g.topological_sort();
    auto vertices = g.get_vertices();
    double dt_ideal = data.calc_ideal_dt(models);

    double t = 0; // текущее время
    double dt = 300; // время по реальным данным
    size_t N = static_cast<int>(T / dt);

    std::stringstream filename;
    filename << path << "Rho_3" << ".csv";
    std::ofstream output(filename.str());

    std::map<size_t, vector<double>> vertices_density;

    for (size_t index = 0; index < N; ++index) {

        // Учет краевых условий и расходов на новом шаге
        for (size_t i = 1; i < pipes.size(); ++i) { // Костыль
            double q_pipe = i == 1
                ? data.volflow_in[index]
                : data.volflow_out[index];
            Q[i] = vector<double>(pipes[i].profile.getPointCount(), q_pipe);
        }
        boundaries[0] = data.density_in[index];

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

                double density_vertex_out = -1; // фиктивное краевое условие на выходе ребра

                quickest_ultimate_fv_solver solver(models[edge], buffers[edge]);
                solver.step(dt, density_vertex, density_vertex_out);

            }
        }
        t += dt;

        layer_t& next = buffers[0].current();
        next.vars.print(t, output);

        for (auto& buffer : buffers) {
            buffer.advance(+1);
        }

    }
    output.flush();
    output.close();

}
