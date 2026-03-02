// TriCommon.hpp
#pragma once
#include <array>
#include <cstddef>
#include <ostream>

namespace cgshop {

#ifndef CGSHOP_POINT_DEFINED
#define CGSHOP_POINT_DEFINED
    struct Point { double x{0}, y{0};
        int lbl;
        friend std::ostream &operator<<(std::ostream &os, const Point &point) {
            os << "x: " << point.x << " y: " << point.y << " label: " << point.lbl;
            return os;
        }
    };
#endif

    struct Edge {int u, v;};

#ifndef CGSHOP_EDGEKEY_DEFINED
#define CGSHOP_EDGEKEY_DEFINED
    struct EdgeKey {
        int a, b; // normalized a<b
        EdgeKey() : a(-1), b(-1) {}
        EdgeKey(int u, int v) { if (u < v) { a = u; b = v; } else { a = v; b = u; } }
        bool operator==(const EdgeKey& o) const noexcept { return a==o.a && b==o.b; }

        friend std::ostream &operator<<(std::ostream &os, const EdgeKey &key) {
            os << "u: " << key.a << " v: " << key.b;
            return os;
        }
    };
    struct EdgeKeyHash {
        size_t operator()(const EdgeKey& e) const noexcept {
            return (static_cast<size_t>(e.a) << 32) ^ static_cast<size_t>(e.b);
        }
    };
#endif

    using Face = std::array<int,3>;

} // namespace cgshop
