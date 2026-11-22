// Each surface's gradient against central differences of its function, at
// scattered points. Exit 1 if one differs.

#include "tpms.hpp"

#include <cstdio>

int main() {
    using namespace tpms;
    int bad = 0;
    for (const Surface& s : surfaces())
        for (int k = 0; k < 200; ++k) {
            const Vec3 p{0.37 * k, 1.3 + 0.61 * k, 2.1 + 0.83 * k};
            const double h = 1e-6;
            const Vec3 d{(s.f(p.x + h, p.y, p.z) - s.f(p.x - h, p.y, p.z)) / (2 * h),
                         (s.f(p.x, p.y + h, p.z) - s.f(p.x, p.y - h, p.z)) / (2 * h),
                         (s.f(p.x, p.y, p.z + h) - s.f(p.x, p.y, p.z - h)) / (2 * h)};
            const Vec3 g = s.grad(p.x, p.y, p.z);
            if (norm(g - d) > 1e-6 * (1 + norm(d))) {
                std::printf("%s: gradient off at (%g, %g, %g)\n", s.name.c_str(), p.x, p.y, p.z);
                ++bad;
                break;
            }
        }
    return bad ? 1 : 0;
}
