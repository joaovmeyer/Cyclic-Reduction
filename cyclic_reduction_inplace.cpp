#include <vector>
#include <iostream>
#include <ranges>
#include <cmath>
#include <numeric>
#include <cstddef>
#include <omp.h>








std::vector<double> cyclic_reduction_inplace2(auto& a, auto& b, auto& c, auto& rhs, std::size_t stride = 1) {

    const std::size_t n = b.size();

    if (stride >= n) {
        return std::vector<double>(b.size(), rhs[0] / b[0]);
    }

    for (std::size_t i = 0; i < n; i += 2 * stride) {
        // a[i - 1] * x[i - stride] + b[i] * x[i] + c[i] * x[i + stride] = rhs[i]
        // x[i - stride] = (rhs[i - stride] - a[i - stride - 1] * x[i - 2 * stride] - c[i - stride] * x[i]) / b[i - stride]
        // x[i + stride] = (rhs[i + stride] - a[i + stride - 1] * x[i] - c[i + stride] * x[i + 2 * stride]) / b[i + stride]

        // => a[i - 1] * (rhs[i - stride] - a[i - stride - 1] * x[i - 2 * stride] - c[i - stride] * x[i]) / b[i - stride]
        //     + b[i] * x[i] + c[i] * (rhs[i + stride] - a[i + stride - 1] * x[i] - c[i + stride] * x[i + 2 * stride]) / b[i + stride]
        //     = rhs[i]

        // updated linear system:
        // => (-a[i - 1] * a[i - stride - 1] / b[i - stride]) * x[i - 2 * stride]
        //     + (b[i] - a[i - 1] * c[i - stride] / b[i - stride] - c[i] * a[i + stride - 1] / b[i + stride]) * x[i]
        //     + (-c[i] * c[i + stride] / b[i + stride]) * x[i + 2 * stride]
        //     = rhs[i] - a[i - 1] * rhs[i - stride] / b[i - stride] - c[i] * rhs[i + stride] / b[i + stride]

        b[i] = b[i] - (i >= stride ? a[i - 1] * c[i - stride] / b[i - stride] : 0.0) - (i + stride < n ? c[i] * a[i + stride - 1] / b[i + stride] : 0.0);
        rhs[i] = rhs[i] - (i >= stride ? a[i - 1] * rhs[i - stride] / b[i - stride] : 0.0) - (i + stride < n ? c[i] * rhs[i + stride] / b[i + stride] : 0.0);

        if (i > stride) a[i - 1] = -a[i - 1] * a[i - stride - 1] / b[i - stride];
        if (i + 2 * stride < n) c[i] = -c[i] * c[i + stride] / b[i + stride];
    }

    auto sol = cyclic_reduction_inplace2(a, b, c, rhs, stride * 2);

    for (std::size_t i = stride; i < n; i += 2 * stride) {
        sol[i] = (rhs[i] - a[i - 1] * sol[i - stride] - (i + stride < n ? c[i] * sol[i + stride] : 0.0)) / b[i];
    }

    // NRVO?
    return sol;
}




std::vector<double> cyclic_reduction_inplace(auto& a, auto& b, auto& c, auto& rhs, std::size_t stride = 1) {

    const std::size_t n = b.size();
    
    std::size_t start_idx = 2 * stride - 1;

    if (start_idx >= n) {
        return std::vector<double>(b.size(), rhs[stride - 1] / b[stride - 1]);
    }
	
	for (std::size_t i = start_idx; i < n; i += 2 * stride) {
		// a[i - 1] * x[i - stride] + b[i] * x[i] + c[i] * x[i + stride] = rhs[i]
		// x[i - stride] = (rhs[i - stride] - a[i - stride - 1] * x[i - 2 * stride] - c[i - stride] * x[i]) / b[i - stride]
		// x[i + stride] = (rhs[i + stride] - a[i + stride - 1] * x[i] - c[i + stride] * x[i + 2 * stride]) / b[i + stride]

		// => a[i - 1] * (rhs[i - stride] - a[i - stride - 1] * x[i - 2 * stride] - c[i - stride] * x[i]) / b[i - stride]
		//     + b[i] * x[i] + c[i] * (rhs[i + stride] - a[i + stride - 1] * x[i] - c[i + stride] * x[i + 2 * stride]) / b[i + stride]
		//     = rhs[i]

		// updated linear system:
		// => (-a[i - 1] * a[i - stride - 1] / b[i - stride]) * x[i - 2 * stride]
		//     + (b[i] - a[i - 1] * c[i - stride] / b[i - stride] - c[i] * a[i + stride - 1] / b[i + stride]) * x[i]
		//     + (-c[i] * c[i + stride] / b[i + stride]) * x[i + 2 * stride]
		//     = rhs[i] - a[i - 1] * rhs[i - stride] / b[i - stride] - c[i] * rhs[i + stride] / b[i + stride]

		b[i] = b[i] - (i >= stride ? a[i - 1] * c[i - stride] / b[i - stride] : 0.0) - (i + stride < n ? c[i] * a[i + stride - 1] / b[i + stride] : 0.0);
		rhs[i] = rhs[i] - (i >= stride ? a[i - 1] * rhs[i - stride] / b[i - stride] : 0.0) - (i + stride < n ? c[i] * rhs[i + stride] / b[i + stride] : 0.0);

		if (i >= 2 * stride) a[i - 1] = -a[i - 1] * a[i - stride - 1] / b[i - stride];
		if (i + 2 * stride < n) c[i] = -c[i] * c[i + stride] / b[i + stride];
	}
	
	auto sol = cyclic_reduction_inplace(a, b, c, rhs, stride * 2);

    for (std::size_t i = start_idx - stride; i < n; i += 2 * stride) {
        sol[i] = (rhs[i] - (i >= stride ? a[i - 1] * sol[i - stride] : 0.0) - (i + stride < n ? c[i] * sol[i + stride] : 0.0)) / b[i];
    }

    // NRVO?
    return sol;
}



int main() {

	std::size_t n = 11;
	std::vector a(n - 1, -1.0);
	std::vector b(n, 2.0);
	std::vector c(n - 1, -1.0);

	std::vector expected_sol(n, 0.0);
	std::iota(expected_sol.begin(), expected_sol.end(), 0.0);

	std::vector<double> rhs(b.size());
	for (std::size_t i = 0; i < rhs.size(); ++i) {
		rhs[i] = (i > 0 ? a[i - 1] * expected_sol[i - 1] : 0.0) + b[i] * expected_sol[i] + (i < rhs.size() - 1 ? c[i] * expected_sol[i + 1] : 0.0);
	}

	auto found_sol = cyclic_reduction_inplace(a, b, c, rhs);

	 for (const auto& [expected, found] : std::ranges::views::zip(expected_sol, found_sol)) {
	 	std::cout << expected - found << "\n";
	}

	return 0;
}







