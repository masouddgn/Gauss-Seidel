#ifndef TRIANGLE_H
#define TRIANGLE_H
#include <memory>

struct Vertex {
    double x, y;
};

struct Face {
    int v1, v2, v3;
};
#endif // DAHA_H