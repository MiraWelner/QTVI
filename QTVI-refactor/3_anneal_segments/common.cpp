#include "common.hpp"

// MATLAB-compatible round: "round half to even" (banker's rounding).
// std::round() rounds half away from zero which differs at exactly *.5.
// MATLAB: round(0.5)=0, round(1.5)=2, round(2.5)=2, round(-0.5)=0
double matlab_round(double x) {
    double r = std::round(x);
    double diff = x - std::floor(x);
    if (std::abs(diff - 0.5) < 1e-12) {
        double f = std::floor(x);
        if (std::fmod(std::abs(f), 2.0) < 0.5)
            r = f;
        else
            r = f + 1.0;
    }
    return r;
}

// MATLAB closest_idx: finds index in uniform time vector nearest to target.
// Equivalent to round(target_time * sr) + 1 using MATLAB rounding.
// Returns 1-based index.
uint64_t closest_idx(double target_time, double sr) {
    double raw = target_time * sr;
    double rounded = matlab_round(raw);
    if (rounded < 0.0) rounded = 0.0;
    return static_cast<uint64_t>(rounded) + 1;
}