/**
  Levenberg-Marquardt Quantum Circuit Optimizer (LMQCO): Matches parameterized,
  layer-based quantum circuits over a three-qubit linear qubit topology to a
  target unitary matrix using the Levenberg-Marquardt nonlinear least squares
  optimization method. See README.md for additional information.
  @file main.cpp
  @version 0.0.1 July 26th, 2021
  @author Max Aksel Bowman (mbowman@anl.gov)
*/

#include <iostream>
#include <random>
#include <ctime>
#include <stdio.h>
#include <math.h>
#include <complex>
#include <mpi.h>
#include <Eigen/Dense>
#include "ceres/ceres.h"
#include "glog/logging.h"

using Eigen::Matrix;
using Eigen::MatrixXd;
using ceres::Jet;
using ceres::AutoDiffCostFunction;
using ceres::CostFunction;
using ceres::Problem;
using ceres::Solve;
using ceres::Solver;

template <typename T>
struct ccomplex {
    T real;
    T imag;

    ccomplex() {
      real = (T)(0.0);
      imag = (T)(0.0);
    }

    ccomplex(int i) {
      real = (T)(i);
      imag = (T)(0.0);
    }

    ccomplex(T c_real, T c_imag) {
        real = c_real;
        imag = c_imag;
    }

    inline ccomplex<T> operator+(const ccomplex<T>& obj) const {
        return ccomplex<T>(real + obj.real, imag + obj.imag);
    }

    inline ccomplex<T> operator*(const T& scalar) const {
        return ccomplex<T>(scalar * real, scalar * imag);
    }

    inline ccomplex<T> operator-(const ccomplex<T>& obj) const {
        return *this + ccomplex<T>((T)(-1.0), (T)(0.0)) * obj;
    }

    inline ccomplex<T> operator*(const ccomplex<T>& obj) const {
        return ccomplex<T>(real * obj.real - imag * obj.imag,
                           real * obj.imag + imag * obj.real);
    }

    inline ccomplex<T> operator/(const T& scalar) const {
        return *this * ((T)(1.0) / scalar);
    }

    inline ccomplex<T>& operator+=(const ccomplex<T>& obj) {
      real = real + obj.real;
      imag = imag + obj.imag;
      return *this;
    }

    inline ccomplex<T>& operator-=(const ccomplex<T>& obj) {
      real = real - obj.real;
      imag = imag - obj.imag;
      return *this;
    }

    inline ccomplex<T>& operator*=(const ccomplex<T>& obj) {
      real = real * obj.real - imag * obj.imag;
      imag = real * obj.imag + imag * obj.real;
      return *this;
    }

    friend std::ostream& operator<<(std::ostream& os, const ccomplex<T>& obj) {
      os << obj.real << '+' << obj.imag << 'j';
      return os;
    }
};

namespace Eigen {
    template<typename T>
    struct NumTraits<ccomplex<T>>: NumTraits<double> {
        typedef T Real;
        typedef ccomplex<T> NonInteger;
        typedef ccomplex<T> Nested;

        enum {
            IsComplex = 0,
            IsInteger = 0,
            IsSigned = 1,
            RequireInitialization = 1,
            ReadCost = 1,
            AddCost = 3,
            MulCost = 3
        };
    };
}

template <typename T>
Matrix<ccomplex<T>, 8, 8> u3_gate(const T theta, const T phi, const T lambda, int qubit) {
  Matrix<ccomplex<T>, 8, 8> unitary;

  // Useful constants for computing unitary
  ccomplex<T> one = ccomplex<T>((T)(1.0), (T)(0.0));
  ccomplex<T> zero = ccomplex<T>((T)(0.0), (T)(0.0));
  T half_theta = (T)(0.5)*theta;

  // Create complex variables a, b, c, and d corresponding to 2x2 single-qubit unitary [a b; c d];
  ccomplex<T> a = ccomplex<T>(cos(half_theta), (T)(0.0));
  ccomplex<T> b = ccomplex<T>(-cos(lambda)*sin(half_theta), -sin(lambda)*sin(half_theta));
  ccomplex<T> c = ccomplex<T>(cos(phi)*sin(half_theta), sin(phi)*sin(half_theta));
  ccomplex<T> d = ccomplex<T>(cos(phi+lambda)*cos(half_theta), sin(phi+lambda)*cos(half_theta));

  // Kronecker products were taken by hand due to the Kronercker product operation technically being unsupported by Eigen
  if (qubit == 0) {
    unitary << a, b, zero, zero, zero, zero, zero, zero,
        c, d, zero, zero, zero, zero, zero, zero,
        zero, zero, a, b, zero, zero, zero, zero,
        zero, zero, c, d, zero, zero, zero, zero,
        zero, zero, zero, zero, a, b, zero, zero,
        zero, zero, zero, zero, c, d, zero, zero,
        zero, zero, zero, zero, zero, zero, a, b,
        zero, zero, zero, zero, zero, zero, c, d;
  }
  else if (qubit == 1) {
      unitary << a, zero, b, zero, zero, zero, zero, zero,
          zero, a, zero, b, zero, zero, zero, zero,
          c, zero, d, zero, zero, zero, zero, zero,
          zero, c, zero, d, zero, zero, zero, zero,
          zero, zero, zero, zero, a, zero, b, zero,
          zero, zero, zero, zero, zero, a, zero, b,
          zero, zero, zero, zero, c, zero, d, zero,
          zero, zero, zero, zero, zero, c, zero, d;
  }
  else if (qubit == 2) {
    unitary << a, zero, zero, zero, b, zero, zero, zero,
        zero, a, zero, zero, zero, b, zero, zero,
        zero, zero, a, zero, zero, zero, b, zero,
        zero, zero, zero, a, zero, zero, zero, b,
        c, zero, zero, zero, d, zero, zero, zero,
        zero, c, zero, zero, zero, d, zero, zero,
        zero, zero, c, zero, zero, zero, d, zero,
        zero, zero, zero, c, zero, zero, zero, d;
  }

  return unitary;
}

template <typename T>
Matrix<ccomplex<T>, 8, 8> linear_mq_gate(int gate_type) {
    Matrix<ccomplex<T>, 8, 8> unitary;
    ccomplex<T> one = ccomplex<T>((T)(1.0), (T)(0.0));
    ccomplex<T> zero = ccomplex<T>((T)(0.0), (T)(0.0));

    if (gate_type == 0) {
        // Apply a CNOT gate between qubits 0 and 1 on a three-qubit circuit
        unitary << one, zero, zero, zero, zero, zero, zero, zero,
            zero, zero, zero, one, zero, zero, zero, zero,
            zero, zero, one, zero, zero, zero, zero, zero,
            zero, one, zero, zero, zero, zero, zero, zero,
            zero, zero, zero, zero, one, zero, zero, zero,
            zero, zero, zero, zero, zero, zero, zero, one,
            zero, zero, zero, zero, zero, zero, one, zero,
            zero, zero, zero, zero, zero, one, zero, zero;
    }
    else if (gate_type == 1) {
        // Apply a CNOT gate between qubits 0 and 2 on a three-qubit circuit
        unitary << one, zero, zero, zero, zero, zero, zero, zero,
            zero, zero, zero, zero, zero, one, zero, zero,
            zero, zero, one, zero, zero, zero, zero, zero,
            zero, zero, zero, zero, zero, zero, zero, one,
            zero, zero, zero, zero, one, zero, zero, zero,
            zero, one, zero, zero, zero, zero, zero, zero,
            zero, zero, zero, zero, zero, zero, one, zero,
            zero, zero, zero, one, zero, zero, zero, zero;
    }
    else if (gate_type == 2) {
        // Apply a FAN-OUT gate with qubit 0 as the control qubit
        unitary << one, zero, zero, zero, zero, zero, zero, zero,
            zero, zero, zero, zero, zero, zero, zero, one,
            zero, zero, one, zero, zero, zero, zero, zero,
            zero, zero, zero, zero, zero, one, zero, zero,
            zero, zero, zero, zero, one, zero, zero, zero,
            zero, zero, zero, one, zero, zero, zero, zero,
            zero, zero, zero, zero, zero, zero, one, zero,
            zero, one, zero, zero, zero, zero, zero, zero;
    }

    return unitary;
}

template <typename T>
Matrix<ccomplex<T>, 8, 8> build_unitary(int circuit_number, int num_layers, const T* const params) {
    const int num_params = (num_layers + 1) * 9;
    ccomplex<T> one = ccomplex<T>((T)(1.0), (T)(0.0));
    ccomplex<T> zero = ccomplex<T>((T)(0.0), (T)(0.0));

    Matrix<ccomplex<T>, 8, 8> unitary;

    unitary << one, zero, zero, zero, zero, zero, zero, zero,
          zero, one, zero, zero, zero, zero, zero, zero,
          zero, zero, one, zero, zero, zero, zero, zero,
          zero, zero, zero, one, zero, zero, zero, zero,
          zero, zero, zero, zero, one, zero, zero, zero,
          zero, zero, zero, zero, zero, one, zero, zero,
          zero, zero, zero, zero, zero, zero, one, zero,
          zero, zero, zero, zero, zero, zero, zero, one; // initialize unitary to identity matrix

    unitary *= u3_gate<T>(params[0], params[1], params[2], 0);
    unitary *= u3_gate<T>(params[3], params[4], params[5], 1);
    unitary *= u3_gate<T>(params[6], params[7], params[8], 2);

    int max_exponent;
    if (circuit_number == 0) {
        max_exponent = 0;
    }
    else {
        max_exponent = (int)(log(circuit_number) / log(3));
    }
    max_exponent = std::max(max_exponent, num_layers - 1); // add "padding" for leading 0-type multi-qubit operations

    int instruction;
    int instruction_index = 0;
    for (int base_exponent = max_exponent; base_exponent >= 0; base_exponent--) {
        // Extract next multi-qubit instruction from the circuit number
        instruction = (int)(circuit_number / pow(3, base_exponent));
        circuit_number = circuit_number % int(pow(3, base_exponent));

        // Add next layer
        unitary *= linear_mq_gate<T>(instruction); // add multi-qubit operation
        for (int qubit = 0; qubit < 3; qubit++) {
            T theta = params[(instruction_index + 1) * 9 + qubit * 3];
            T phi = params[(instruction_index + 1) * 9 + qubit * 3 + 1];
            T lambda = params[(instruction_index + 1) * 9 + qubit * 3 + 2];
            unitary *= u3_gate<T>(theta, phi, lambda, qubit); // add single-qubit operations
        }
        instruction_index++;
    }

    return unitary;
}

struct ResidualFunctor {
    int circuit_number;
    int num_layers;

    ResidualFunctor(int cn, int nl) {
        circuit_number = cn;
        num_layers = nl;
    }

    template <typename T>
    bool operator()(const T* const params, T* residual) const {
        Matrix<ccomplex<T>, 8, 8> cunitary = build_unitary(circuit_number, num_layers, params);
        ccomplex<T> one = ccomplex<T>((T)(1.0), (T)(0.0));
        ccomplex<T> zero = ccomplex<T>((T)(0.0), (T)(0.0));

        for (int i_ind = 0; i_ind < 8; i_ind++) {
            for (int j_ind = 0; j_ind < 8; j_ind++) {
                if ((i_ind == j_ind && i_ind != 3 && i_ind != 7) || (i_ind == 3 && j_ind == 7) || (i_ind == 7 && j_ind == 3)) {
                  residual[2 * (i_ind + 8 * j_ind)] = (cunitary(i_ind, j_ind) - one).real; // even vector elements contain real residuals
                  residual[2 * (i_ind + 8 * j_ind) + 1] = (cunitary(i_ind, j_ind) - one).imag; // odd vector element contain imaginary residuals
                } else {
                  residual[2 * (i_ind + 8 * j_ind)] = cunitary(i_ind, j_ind).real; // even vector elements contain real residuals
                  residual[2 * (i_ind + 8 * j_ind) + 1] = cunitary(i_ind, j_ind).imag; // odd vector element contain imaginary residuals
                }
            }
        }
        return true;
    }
};

template <typename T>
void print_array(T* arr, int length) {
  for (int ind = 0; ind < length; ind++) {
    std::cout << arr[ind];
    if (ind < length - 1) {
      std::cout << ", ";
    }
  }
}

/**
  Returns a pointer to an array of doubles. Array contains random parameters
  used as a starting point for the Levenberg-Marquardt optimization algorithm.
  @param an integer number of parameters to generate (should be a multiple of 3).
  @return a pointer to an array of doubles. Random parameters are of form
          [theta_1, phi_1, lambda_1, ... theta_n, phi_n, lambda_n]
*/
double* generate_initial_parameters(int num) {
  double* random_array = new double[num];
  std::default_random_engine generator(unsigned(time(nullptr)));
  std::uniform_real_distribution<double> theta_distribution(0.0, 4*M_PI);
  std::uniform_real_distribution<double> phi_lam_distribution(0.0, 2*M_PI);
  for (int i = 0; i < num; i++) {
    if (i % 3 == 0) {
      random_array[i] = theta_distribution(generator);
    } else {
      random_array[i] = phi_lam_distribution(generator);
    }
  }
  return random_array;
}

int main(int argc, char** argv) {
  // Initialize the MPI environment
  MPI_Init(NULL, NULL);

  // Get the number, rank, and name of processes from MPI
  int size, rank, name_len;
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Get_processor_name(processor_name, &name_len);

  // Optimization hyperparameters
  const int num_params = 72; // 18 -> 27 -> 36 -> 45 -> 54 -> 63 -> 72 -> 81
  int num_layers = num_params / 9 - 1; // derived from number of parameters
  int num_initial_points = 1; // number of initial points to start from per subproblem (set to 1 if using multiple machines)
  const int circuit_number_min = pow(3, num_layers) * rank / size;
  const int circuit_number_max = pow(3, num_layers) * (rank + 1) / size;

  printf("Optimizing for %d layers (%d parameters) with %d initial guesses per optimization problem.\n", num_layers, num_params, num_initial_points);

  google::InitGoogleLogging(argv[0]);

  for (int circuit_number = circuit_number_min; circuit_number < circuit_number_max; circuit_number++) {
  // for (int circuit_number = 2460; circuit_number <= 2460; circuit_number++) {
    for (int run_index = 0; run_index < num_initial_points; run_index++) {
      printf("Optimizing circuit %d on rank %d on machine %s (run id is %d)\n", circuit_number, rank, processor_name, run_index);
      // Initialize random starting point and print it
      double* opt_params = generate_initial_parameters(num_params); // this array will be updated at each iteration of Levenberg-Marquardt
      double initial_params[num_params];
      memcpy(initial_params, opt_params, sizeof(initial_params));

      // Build nonlinear least squares problem
      Matrix<ccomplex<double>, 8, 8> cunitary = build_unitary(circuit_number, num_layers, opt_params);
      Problem problem;
      CostFunction* cost_function = new AutoDiffCostFunction<ResidualFunctor, 128, num_params>(new ResidualFunctor(circuit_number, num_layers));
      problem.AddResidualBlock(cost_function, nullptr, opt_params);

      // Configure and run Levenberg-Marquardt nonlinear least squares solver
      Solver::Options options;
      // options.minimizer_progress_to_stdout = true;
      options.max_num_iterations = 200;
      options.linear_solver_type = ceres::DENSE_QR;
      options.function_tolerance = 1e-320;
      options.max_num_consecutive_invalid_steps = 50;
      options.gradient_tolerance = 1e-320;
      options.parameter_tolerance = 1e-320;

      Solver::Summary summary;
      Solve(options, &problem, &summary);

      // Save (good) results to disk
      if (summary.final_cost <= 1e-10) {
        std::cout << "Success. Trying to output to a file..." << std::endl;
        char filename[255];
        FILE* file;
        sprintf(filename, "optout_%d_%s_%d.txt", circuit_number, processor_name, run_index);
        file = fopen(filename, "w");
        for (int i = 0; i < num_params - 1; i++) {
          fprintf(file, "%.50f, ", opt_params[i]);
        }
        fprintf(file, "%.50f\n", opt_params[num_params - 1]);
        for (int i = 0; i < num_params - 1; i++) {
          fprintf(file, "%.50f, ", initial_params[i]);
        }
        fprintf(file, "%.50f\n", initial_params[num_params - 1]);
        fprintf(file, "%.50f\n", summary.final_cost);
        fclose(file);
      }
    }
  }

  printf("Rank %d is finished.\n", rank);
  // Finalize the MPI environment
  MPI_Finalize();

  return 0;
}
