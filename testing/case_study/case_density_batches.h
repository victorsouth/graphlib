#pragma once


using graphlib::edge_t;
using graphlib::graph_t;
using graphlib::task_network_propagation_data_t;
using graphlib::net_quickest_ultimate_solver;
using namespace std;

/// @brief Представление временных рядов по объемному расходу и плотности
/// Очень частное представление - один вход, один выход
/// Используется для схемы с разветвлением, поскольку на ответвлении все равно нет данных
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
class DensityPropagation : public ::testing::Test {
protected:
    /// @brief Сеть с движением партий для реального трубопровода с разветвлением
    task_network_propagation_data_t net_data;
protected:
    /// @brief Подготовка к расчету для семейства тестов
    virtual void SetUp() override {
        net_data = task_network_propagation_data_t::default_data();

        double rho_initial = 850;
        net_data.init_buffers(rho_initial);
    }


};

/// @brief Расчетный кейс для движения партий с разветвлением
TEST_F(DensityPropagation, MixDensity) {
    
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
