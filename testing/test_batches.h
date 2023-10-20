#pragma once
using graphlib::edge_t;
using graphlib::graph_t;
using namespace std;

#include <fixed/fixed.h>
#include <pde_solvers/pde_solvers.h>

/// @brief �������� ����� ������ ��� ������ �������� ������� �� ������ QUICKEST-ULTIMATE 
template <size_t Dimension>
struct quickest_ultimate_fv_solver_traits
{
    typedef profile_collection_t<0, Dimension/*���������� - ������*/, 0, 0, 0, 0> var_layer_data;
    typedef profile_collection_t<Dimension /*������ F*/, 0,
        0, 0,
        0, 0> specific_layer;
};

double quickest_ultimate_border_approximation(double U_L, double U_C, double U_R, double hi, double dx, double dt, double v)
{
    double DEL = U_R - U_L;
    double ADEL = abs(DEL);
    double ACURV = abs(U_L + U_R - 2 * U_C);
    if (ACURV >= ADEL) {
        return U_C;
    }
    double Cour = abs((v * dt) / dx);
    double REF = U_L + ((U_C - U_L) / Cour);
    double Ub_linear = (U_C + U_R) / 2;
    double Grad = (U_R - U_C) / dx;
    double Curv = (U_L + U_R - 2 * U_C) / (dx * dx);
    double Ub_correction_first_order = -(dx * Cour * Grad) / 2;
    double Ub_correction_second_order = (dx * dx * (hi - (1 - (Cour * Cour)) / 3) * Curv) / 2;
    double Uf = Ub_linear + Ub_correction_first_order + Ub_correction_second_order;
    if (DEL > 0) {
        if (Uf < U_C) {
            Uf = U_C;
        }
        if (Uf > std::min(REF, U_R)) {
            Uf = std::min(REF, U_R);
        }
    }
    else {
        if (Uf > U_C) {
            Uf = U_C;
        }
        if (Uf < std::max(REF, U_R)) {
            Uf = std::max(REF, U_R);
        }
    }
    return Uf;
}

/// @brief ������ �� ������ QUICKEST-ULTIMATE, ������ ��� ����������� 1!
/// [Leonard 1991]
class quickest_ultimate_fv_solver {
public:
    typedef typename quickest_ultimate_fv_solver_traits<1>::var_layer_data var_layer_data;
    typedef typename quickest_ultimate_fv_solver_traits<1>::specific_layer specific_layer;
    typedef typename fixed_system_types<1>::matrix_type matrix_type;
    typedef typename fixed_system_types<1>::var_type vector_type;
protected:
    /// @brief ����
    pde_t<1>& pde;
    /// @brief �����, ���������� �� ����
    const vector<double>& grid;
    /// @brief ���������� ����� �����
    const size_t n;
    /// @brief ���������� ���� ����������
    const var_layer_data& prev_vars;
    /// @brief ����� (��������������) ���� ����������
    var_layer_data& curr_vars;
    /// @brief ���������� ������������� ���� (������ �� �����! ����� �� � �������?)
    const specific_layer& prev_spec;
    /// @brief ������� ������������� ����
    specific_layer& curr_spec;
public:
    /// @brief ����������� ��� ������ � ������� ������� ����:
    /// ����� � ���� ������ ���� ���� ������� ���������� � ���� ���� ��������� ������
    /// �� ������ ������� current() � previous()
    /// @param pde ����
    /// @param buffer ����� �����
    quickest_ultimate_fv_solver(pde_t<1>& pde,
        custom_buffer_t<composite_layer_t<var_layer_data, specific_layer>>& buffer)
        : quickest_ultimate_fv_solver(pde, buffer.previous(), buffer.current())
    {}

    /// @brief ����������� ��� ������� ����� - 
    /// ����� � ���� ������ ���� ���� ������� ���������� � ���� ���� ��������� ������
    /// @param pde ����
    /// @param prev ���������� ���� (��� ������������)
    /// @param curr ��������� (�����), ��� �������� ��������� ������� ������
    quickest_ultimate_fv_solver(pde_t<1>& pde,
        const composite_layer_t<var_layer_data, specific_layer>& prev,
        composite_layer_t<var_layer_data, specific_layer>& curr)
        : pde(pde)
        , grid(pde.get_grid())
        , n(pde.get_grid().size())
        , prev_vars(prev.vars)
        , curr_vars(curr.vars)
        , prev_spec(std::get<0>(prev.specific))
        , curr_spec(std::get<0>(curr.specific))
    {

    }
    /// @brief ������ ����
    /// @param dt �������� ������ �������
    /// @param u_in ����� ��������� �������
    /// @param u_out ������ ��������� �������
    void step(double dt, double u_in, double u_out) {
        auto& F = curr_spec.point_double[0]; // ������ �� �������� �����
        const auto& U = prev_vars.cell_double[0];
        auto& U_new = curr_vars.cell_double[0];

        // ������ ������� �� ������� �� ������ ��������� �������
        double v_in = pde.getEquationsCoeffs(0, U[0]);
        double v_out = pde.getEquationsCoeffs(F.size() - 1, U[U.size() - 1]);
        if (v_in >= 0) {
            F.front() = v_in * u_in;
        }
        if (v_out <= 0) {
            F.back() = v_out * u_out;
        }

        double v_pipe = pde.getEquationsCoeffs(0, U[0]);//�� ������ ���������, �������� � ������ ������� �� �������� �� �� ����� �������
        // ������ ������� �� ������� �� ������� QUICK
        if (v_pipe >= 0) {
            for (size_t cell = 0; cell < U.size(); ++cell) {
                size_t right_border = cell + 1;
                double Vb = v_pipe; // ������������, ��� �������� �� ������� �� ���� ������ ����� ���� � �� ��
                double Ub;
                if (cell == 0) {
                    Ub = quickest_ultimate_border_approximation(U[cell], U[cell], U[cell + 1], 0, grid[cell + 1] - grid[cell], dt, v_pipe); // ������� U_L = U_C
                }
                else if (cell == U.size() - 1) {
                    Ub = quickest_ultimate_border_approximation(U[cell - 1], U[cell], U[cell], 0, grid[cell + 1] - grid[cell], dt, v_pipe); // ������� U_R = U_C
                }
                else {
                    Ub = quickest_ultimate_border_approximation(U[cell - 1], U[cell], U[cell + 1], 0, grid[cell + 1] - grid[cell], dt, v_pipe); // ������� ������
                }
                F[right_border] = Ub * Vb;
            }
        }
        else {
            for (size_t cell = 0; cell < U.size(); ++cell) {
                size_t left_border = cell;
                double Vb = v_pipe; // ������������, ��� �������� �� ������� �� ���� ������ ����� ���� � �� ��
                double Ub;
                if (cell == 0) {
                    Ub = quickest_ultimate_border_approximation(U[cell + 1], U[cell], U[cell], 0, grid[cell + 1] - grid[cell], dt, v_pipe); // ������� U_L = U_C
                }
                else if (cell == U.size() - 1) {
                    Ub = quickest_ultimate_border_approximation(U[cell], U[cell], U[cell - 1], 0, grid[cell + 1] - grid[cell], dt, v_pipe); // ������� U_R = U_C
                }
                else {
                    Ub = quickest_ultimate_border_approximation(U[cell + 1], U[cell], U[cell - 1], 0, grid[cell + 1] - grid[cell], dt, v_pipe); // ������� ������
                }
                F[left_border] = Ub * Vb;
            }

        }

        for (size_t cell = 0; cell < U.size(); ++cell) {
            double dx = grid[cell + 1] - grid[cell]; // ������ ������ ���������� �����, �� ���� ��..
            U_new[cell] = U[cell] + dt / dx * ((F[cell] - F[cell + 1]));
        }

    }
};

/// @brief ����� (��) ��� ������� quickest_ultimate_fv_solver
class QUICKEST_ULTIMATE_TU : public ::testing::Test {
protected:
    // ������� ����������
    typedef quickest_ultimate_fv_solver_traits<1>::var_layer_data target_var_t;
    typedef quickest_ultimate_fv_solver_traits<1>::specific_layer specific_data_t;

    // ����: ���������� Vars + ������� ������ ��������� Specific
    typedef composite_layer_t<target_var_t, specific_data_t> layer_t;
protected:
    ///// @brief ��������� �����
    //PipeProperties pipe_1;
    //PipeProperties pipe_2;
    //PipeProperties pipe_3;
    ///// @brief ������� �������
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

    /// @brief ���������� � ������� ��� ��������� ������
    virtual void SetUp() override {
        // ���������� ������� ����� - 50��, � ����� ��������� ��� �������� ����� 1��, ��������� 700��

        // ���������� ������� ����� - 50��, � ����� ��������� ��� �������� ����� 1��, ��������� 700��
        simple_pipe_properties simple_pipe_1;
        simple_pipe_properties simple_pipe_2;
        simple_pipe_properties simple_pipe_3;
        simple_pipe_1.length = 168e3;
        simple_pipe_2.length = 100e3;
        simple_pipe_3.length = 263e3;
        //simple_pipe.length = 700e3; // ���� ����� 700��
        simple_pipe_1.diameter = 0.1;
        simple_pipe_2.diameter = 0.1;
        simple_pipe_3.diameter = 0.1;
        //simple_pipe.diameter = 0.514; // ���� ����� 700��
        simple_pipe_1.dx = 1000;
        simple_pipe_2.dx = 1000;
        simple_pipe_3.dx = 1000;
        //simple_pipe.dx = 100; // ���� ����� 700��
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

/// @brief �������� QUICKEST-ULTIMATE, �������� ��������� ��������� ����� �������� ��������
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

    double T = 300000; // ������ �������������
    //double T = 800000; // ������ ������������� (���� ����� 700��)

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

    double t = 0; // ������� �����
    //double dt = 60; // 1 ������
    double dt = Cr * dt_ideal; // ����� � ����� �� �������
    size_t N = static_cast<int>(T / dt);


    std::map<size_t, vector<double>> vertices_density;

    //for (size_t vertex = 0; vertex < n_vertex; ++vertex) {
    for (size_t index = 0; index < N; ++index) {
        for (size_t vertex : V) {

            const auto& E = vertices[vertex];

            double density_vertex;
            if (E.in.empty()) { // ���� � ������� ������� ��� ������� ����
                density_vertex = boundaries.at(vertex); // �� ������ ���� ��������� �������
            }
            else {
                double Summ_Rho_Q = 0;
                double Summ_Q = 0;
                for (const auto& edg : E.in) { // ���� �� 1 � ������, �� ���� �����
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
                // ��� vertex ���������� �����, ������� � ���� ������, �� ��� ����� �������� ���������, �������        
                // ��� vertex ���������� �����, ������� �� ���� �������, 
                // �� ��� ������� ������������ ��� ��������� ������� ������� ��� QUICKEST
                
                //std::stringstream filename;
                //filename << path << "output " << edge << " Cr=" << Cr << ".csv";
                //std::ofstream output(filename.str());

            
                //if (index == 0) {
                //    prevs[edge] = buffers[edge].previous();
                //    prevs[edge].vars.print(t, output);
                //}

                //for (const auto& vert : vertices) {
                //    if (vert.first == vertex) { // ����� ������� �������.
                //        if (vert.second.in.empty()) { // ���� � ������� ������� ��� ������� ����
                //            Rho_smesi = rho_in_1; // �� ������ ���� ��������� �������
                //        }
                //        else {
                //            double Summ_Rho_Q = 0;
                //            double Summ_Q = 0;
                //            for (const auto& j : vert.second.in) { // ������ �� ���� �������� ������
                //                Summ_Rho_Q = Summ_Rho_Q + (nexts[j].vars.cell_double[0][nexts[j].vars.cell_double[0].size() - 1] * Q[j][Q[j].size() - 1]);
                //                Summ_Q = Summ_Q + (Q[j][Q[j].size() - 1]);
                //            }
                //            Rho_smesi = Summ_Rho_Q / Summ_Q;
                //        }
                //    }
                //}

                double density_vertex_out = -1; // ��������� ������� ������� �� ������ �����

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

TEST(Batches, develop)
{
    vector<edge_t> edges{ edge_t(0, 1), edge_t(1, 2), edge_t(1, 3) };

    graph_t g(edges);

    auto [V, E] = g.topological_sort();



}