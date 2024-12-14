#include "FPToolkit.c"

int window_size;
double gap = 20;

void grid() {
    int i;
    for (i = 0; i < window_size; i += gap) {
        G_line(i, 0, i, window_size);
        G_line(0, i, window_size, i);
    }
}

int intersect_2_lines(double A[2], double B[2], double C[2], double D[2], double intersection[2]) {
    double x1 = A[0], y1 = A[1];
    double x2 = B[0], y2 = B[1];
    double x3 = C[0], y3 = C[1];
    double x4 = D[0], y4 = D[1];

    // find denominators
    double denom = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);

    // check for parallel lines (denominator = 0)
    if (denom == 0) {
        return 0; // lines don't intersect
    }

    // find intersection point
    intersection[0] = ((x1 * y2 - y1 * x2) * (x3 - x4) - (x1 - x2) * (x3 * y4 - y3 * x4)) / denom;
    intersection[1] = ((x1 * y2 - y1 * x2) * (y3 - y4) - (y1 - y2) * (x3 * y4 - y3 * x4)) / denom;

    // check for intersection point in lines
    if (intersection[0] < fmin(x1, x2) || intersection[0] > fmax(x1, x2) ||
        intersection[0] < fmin(x3, x4) || intersection[0] > fmax(x3, x4) ||
        intersection[1] < fmin(y1, y2) || intersection[1] > fmax(y1, y2) ||
        intersection[1] < fmin(y3, y4) || intersection[1] > fmax(y3, y4)) {
        return 0; // intersection outside of lines
    }

    return 1; // lines intersect at point
}

int main() {
    double a[2], b[2], c[2], d[2];
    double intersect[2];
    double signal, xi, yi;
    char q;

    printf("Click two points to determine a line segment,\n");
    printf("then click two more for another line segment.\n");

    window_size = 800;
    G_init_graphics(window_size, window_size);
    G_rgb(0, 0, 0);
    G_clear();
    G_rgb(0.5, 0.5, 0.5);
    grid();
    G_rgb(0, 1, 0);

    G_wait_click(a);
    a[0] = gap * floor((a[0] + 0.5 * gap) / gap);
    a[1] = gap * floor((a[1] + 0.5 * gap) / gap);
    G_fill_circle(a[0], a[1], 3);

    G_wait_click(b);
    b[0] = gap * floor((b[0] + 0.5 * gap) / gap);
    b[1] = gap * floor((b[1] + 0.5 * gap) / gap);
    G_fill_circle(b[0], b[1], 3);

    G_line(a[0], a[1], b[0], b[1]);

    G_wait_click(c);
    c[0] = gap * floor((c[0] + 0.5 * gap) / gap);
    c[1] = gap * floor((c[1] + 0.5 * gap) / gap);
    G_fill_circle(c[0], c[1], 3);

    G_wait_click(d);
    d[0] = gap * floor((d[0] + 0.5 * gap) / gap);
    d[1] = gap * floor((d[1] + 0.5 * gap) / gap);
    G_fill_circle(d[0], d[1], 3);

    G_line(c[0], c[1], d[0], d[1]);

    signal = intersect_2_lines(a, b, c, d, intersect);
    if (signal == 0) {
        G_rgb(1, 0, 0);
        G_draw_string("The two lines do NOT intersect in a unique point.", 20, 20);
    } else {
        xi = intersect[0];
        yi = intersect[1];

        G_rgb(1, 1, 0);
        G_fill_circle(xi, yi, 5); // draw intersection point
    }

    G_display_image();
    q = G_wait_key();
}

