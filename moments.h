typedef struct {
    float pos[3];
    float m;
    float x, y, z;
    float x2, xy, y2, xz, yz, z2;
    float x3, x2y, xy2, y3, x2z, xyz, y2z, xz2, yz2, z3;
    float x4, x3y, x2y2, xy3, y4, x3z, x2yz, xy2z, y3z, x2z2, xyz2, y2z2, xz3, yz3, z4;
    float x5, x4y, x3y2, x2y3, xy4, y5, x4z, x3yz, x2y2z, xy3z, y4z, x3z2, x2yz2, xy2z2, y3z2, x2z3, xyz3, y2z3, xz4, yz4, z5;
    float x6, x5y, x4y2, x3y3, x2y4, xy5, y6, x5z, x4yz, x3y2z, x2y3z, xy4z, y5z, x4z2, x3yz2, x2y2z2, xy3z2, y4z2, x3z3, x2yz3, xy2z3, y3z3, x2z4, xyz4, y2z4, xz5, yz5, z6;
    float x7, x6y, x5y2, x4y3, x3y4, x2y5, xy6, y7, x6z, x5yz, x4y2z, x3y3z, x2y4z, xy5z, y6z, x5z2, x4yz2, x3y2z2, x2y3z2, xy4z2, y5z2, x4z3, x3yz3, x2y2z3, xy3z3, y4z3, x3z4, x2yz4, xy2z4, y3z4, x2z5, xyz5, y2z5, xz6, yz6, z7;
    float x8, x7y, x6y2, x5y3, x4y4, x3y5, x2y6, xy7, y8, x7z, x6yz, x5y2z, x4y3z, x3y4z, x2y5z, xy6z, y7z, x6z2, x5yz2, x4y2z2, x3y3z2, x2y4z2, xy5z2, y6z2, x5z3, x4yz3, x3y2z3, x2y3z3, xy4z3, y5z3, x4z4, x3yz4, x2y2z4, xy3z4, y4z4, x3z5, x2yz5, xy2z5, y3z5, x2z6, xyz6, y2z6, xz7, yz7, z8;
} cartesian_moments;

typedef struct {
    double pos[3];
    double m;
    double x, y, z;
    double x2, xy, y2, xz, yz, z2;
    double x3, x2y, xy2, y3, x2z, xyz, y2z, xz2, yz2, z3;
    double x4, x3y, x2y2, xy3, y4, x3z, x2yz, xy2z, y3z, x2z2, xyz2, y2z2, xz3, yz3, z4;
    double x5, x4y, x3y2, x2y3, xy4, y5, x4z, x3yz, x2y2z, xy3z, y4z, x3z2, x2yz2, xy2z2, y3z2, x2z3, xyz3, y2z3, xz4, yz4, z5;
    double x6, x5y, x4y2, x3y3, x2y4, xy5, y6, x5z, x4yz, x3y2z, x2y3z, xy4z, y5z, x4z2, x3yz2, x2y2z2, xy3z2, y4z2, x3z3, x2yz3, xy2z3, y3z3, x2z4, xyz4, y2z4, xz5, yz5, z6;
    double x7, x6y, x5y2, x4y3, x3y4, x2y5, xy6, y7, x6z, x5yz, x4y2z, x3y3z, x2y4z, xy5z, y6z, x5z2, x4yz2, x3y2z2, x2y3z2, xy4z2, y5z2, x4z3, x3yz3, x2y2z3, xy3z3, y4z3, x3z4, x2yz4, xy2z4, y3z4, x2z5, xyz5, y2z5, xz6, yz6, z7;
    double x8, x7y, x6y2, x5y3, x4y4, x3y5, x2y6, xy7, y8, x7z, x6yz, x5y2z, x4y3z, x3y4z, x2y5z, xy6z, y7z, x6z2, x5yz2, x4y2z2, x3y3z2, x2y4z2, xy5z2, y6z2, x5z3, x4yz3, x3y2z3, x2y3z3, xy4z3, y5z3, x4z4, x3yz4, x2y2z4, xy3z4, y4z4, x3z5, x2yz5, xy2z5, y3z5, x2z6, xyz6, y2z6, xz7, yz7, z8;
} cartesian_moments_d;
