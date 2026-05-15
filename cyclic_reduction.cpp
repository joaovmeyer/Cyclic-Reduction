#include <vector>
#include <iostream>




// tridiagonal system a[i - 1] * x[i - 1] + b[i] * x[i] + c[i] * x[i + 1] = rhs[i]
std::vector<double> cyclic_reduction(const std::vector<double>& a, const std::vector<double>& b, const std::vector<double>& c, const std::vector<double>& rhs) {

    const std::size_t n = b.size();

    if (n == 1) {
        return { rhs[0] / b[0] };
    }

    const std::size_t new_n = (n + 1) / 2;
    std::vector<double> new_a(new_n - 1);
    std::vector<double> new_b(new_n);
    std::vector<double> new_c(new_n - 1);
    std::vector<double> new_rhs(new_n);

    // independent iterations
    for (std::size_t i = 0; i < n; i += 2) {
        // start from equation i: a[i - 1] * x[i - 1] + b[i] * x[i] + c[i] * x[i + 1] = rhs[i]
        // isolate x[i - 1] and x[i + 1] from equations i - 1 and i + 1:
        // x[i - 1] = (rhs[i - 1] - a[i - 2] * x[i - 2] - c[i - 1] * x[i]) / b[i - 1]
        // x[i + 1] = (rhs[i + 1] - a[i] * x[i] - c[i + 1] * x[i + 2]) / b[i + 1]

        // substitute it back into i-th equation
        // => a[i - 1] * (rhs[i - 1] - a[i - 2] * x[i - 2] - c[i - 1] * x[i]) / b[i - 1]
        //     + b[i] * x[i] + c[i] * (rhs[i + 1] - a[i] * x[i] - c[i + 1] * x[i + 2]) / b[i + 1]
        //     = rhs[i]

        // rearranging, we see that we're left with a new linear system half the size: 
        // => (-a[i - 1] * a[i - 2] / b[i - 1]) * x[i - 2]
        //     + (b[i] - a[i - 1] * c[i - 1] / b[i - 1] - c[i] * a[i] / b[i + 1]) * x[i]
        //     + (-c[i] * c[i + 1] / b[i + 1]) * x[i + 2]
        //     = rhs[i] - a[i - 1] * rhs[i - 1] / b[i - 1] - c[i] * rhs[i + 1] / b[i + 1]
        

        if (i > 1) new_a[i / 2 - 1] = -a[i - 1] * a[i - 2] / b[i - 1];
        new_b[i / 2] = b[i] - (i > 0 ? a[i - 1] * c[i - 1] / b[i - 1] : 0.0) - (i + 1 < n ? c[i] * a[i] / b[i + 1] : 0.0);
        new_rhs[i / 2] = rhs[i] - (i > 0 ? a[i - 1] * rhs[i - 1] / b[i - 1] : 0.0) - (i + 1 < n ? c[i] * rhs[i + 1] / b[i + 1] : 0.0);

        if (i + 2 < n) new_c[i / 2] = -c[i] * c[i + 1] / b[i + 1];
    }

    auto new_sol = cyclic_reduction(new_a, new_b, new_c, new_rhs);

    // this is also all independent
    std::vector<double> sol(n);
    for (std::size_t i = 0; i < new_n; ++i) {
        sol[i * 2] = new_sol[i];
    }
    for (std::size_t i = 1; i < n; i += 2) {
        sol[i] = (rhs[i] - a[i - 1] * new_sol[(i - 1) / 2] - (i + 1 < n ? c[i] * new_sol[(i + 1) / 2] : 0.0)) / b[i];
    }

    return sol;
}








int main() {

    std::vector a{ 1.0, 1.0, -3.0 };
    std::vector b{ 4.0, 5.0, 2.0, -5.0 };
    std::vector c{ -1.0, 1.0, -10.0 };

    std::vector x_exact{ 1.0, 2.0, 5.0, -20.0 };

    auto n = b.size();
    std::vector<double> rhs(n);
    rhs[0] = b[0] * x_exact[0] + c[0] * x_exact[1];
    for (std::size_t i = 1; i < n - 1; ++i) {
        rhs[i] = a[i - 1] * x_exact[i - 1] + b[i] * x_exact[i] + c[i] * x_exact[i + 1];
    }
    rhs[n - 1] = a[n - 2] * x_exact[n - 2] + b[n - 1] * x_exact[n - 1];

    auto sol = cyclic_reduction(a, b, c, rhs);

    std::cout << "Computed:\n";
    for (auto v : sol)
        std::cout << v << ", ";

    std::cout << "\n\nExpected:\n";
    for (auto v : x_exact)
        std::cout << v << ", ";

    return 0;
}
