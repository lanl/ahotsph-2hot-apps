typedef struct {
    float pos[3];
    float m;
    float x, y, z;
    float x2, xy, y2, xz, yz, z2;
    float x3, x2y, xy2, y3, x2z, xyz, y2z, xz2, yz2, z3;
    float x4, x3y, x2y2, xy3, y4, x3z, x2yz, xy2z, y3z, x2z2, xyz2, y2z2, xz3, yz3, z4;
    float x5, x4y, x3y2, x2y3, xy4, y5, x4z, x3yz, x2y2z, xy3z, y4z, x3z2, x2yz2, xy2z2, y3z2, x2z3,
        xyz3, y2z3, xz4, yz4, z5;
    float x6, x5y, x4y2, x3y3, x2y4, xy5, y6, x5z, x4yz, x3y2z, x2y3z, xy4z, y5z, x4z2, x3yz2,
        x2y2z2, xy3z2, y4z2, x3z3, x2yz3, xy2z3, y3z3, x2z4, xyz4, y2z4, xz5, yz5, z6;
    float x7, x6y, x5y2, x4y3, x3y4, x2y5, xy6, y7, x6z, x5yz, x4y2z, x3y3z, x2y4z, xy5z, y6z, x5z2,
        x4yz2, x3y2z2, x2y3z2, xy4z2, y5z2, x4z3, x3yz3, x2y2z3, xy3z3, y4z3, x3z4, x2yz4, xy2z4,
        y3z4, x2z5, xyz5, y2z5, xz6, yz6, z7;
    float x8, x7y, x6y2, x5y3, x4y4, x3y5, x2y6, xy7, y8, x7z, x6yz, x5y2z, x4y3z, x3y4z, x2y5z,
        xy6z, y7z, x6z2, x5yz2, x4y2z2, x3y3z2, x2y4z2, xy5z2, y6z2, x5z3, x4yz3, x3y2z3, x2y3z3,
        xy4z3, y5z3, x4z4, x3yz4, x2y2z4, xy3z4, y4z4, x3z5, x2yz5, xy2z5, y3z5, x2z6, xyz6, y2z6,
        xz7, yz7, z8;
    float x9, x8y, x7y2, x6y3, x5y4, x4y5, x3y6, x2y7, xy8, y9, x8z, x7yz, x6y2z, x5y3z, x4y4z,
        x3y5z, x2y6z, xy7z, y8z, x7z2, x6yz2, x5y2z2, x4y3z2, x3y4z2, x2y5z2, xy6z2, y7z2, x6z3,
        x5yz3, x4y2z3, x3y3z3, x2y4z3, xy5z3, y6z3, x5z4, x4yz4, x3y2z4, x2y3z4, xy4z4, y5z4, x4z5,
        x3yz5, x2y2z5, xy3z5, y4z5, x3z6, x2yz6, xy2z6, y3z6, x2z7, xyz7, y2z7, xz8, yz8, z9;
    float x10, x9y, x8y2, x7y3, x6y4, x5y5, x4y6, x3y7, x2y8, xy9, y10, x9z, x8yz, x7y2z, x6y3z,
        x5y4z, x4y5z, x3y6z, x2y7z, xy8z, y9z, x8z2, x7yz2, x6y2z2, x5y3z2, x4y4z2, x3y5z2, x2y6z2,
        xy7z2, y8z2, x7z3, x6yz3, x5y2z3, x4y3z3, x3y4z3, x2y5z3, xy6z3, y7z3, x6z4, x5yz4, x4y2z4,
        x3y3z4, x2y4z4, xy5z4, y6z4, x5z5, x4yz5, x3y2z5, x2y3z5, xy4z5, y5z5, x4z6, x3yz6, x2y2z6,
        xy3z6, y4z6, x3z7, x2yz7, xy2z7, y3z7, x2z8, xyz8, y2z8, xz9, yz9, z10;
    float x11, x10y, x9y2, x8y3, x7y4, x6y5, x5y6, x4y7, x3y8, x2y9, xy10, y11, x10z, x9yz, x8y2z,
        x7y3z, x6y4z, x5y5z, x4y6z, x3y7z, x2y8z, xy9z, y10z, x9z2, x8yz2, x7y2z2, x6y3z2, x5y4z2,
        x4y5z2, x3y6z2, x2y7z2, xy8z2, y9z2, x8z3, x7yz3, x6y2z3, x5y3z3, x4y4z3, x3y5z3, x2y6z3,
        xy7z3, y8z3, x7z4, x6yz4, x5y2z4, x4y3z4, x3y4z4, x2y5z4, xy6z4, y7z4, x6z5, x5yz5, x4y2z5,
        x3y3z5, x2y4z5, xy5z5, y6z5, x5z6, x4yz6, x3y2z6, x2y3z6, xy4z6, y5z6, x4z7, x3yz7, x2y2z7,
        xy3z7, y4z7, x3z8, x2yz8, xy2z8, y3z8, x2z9, xyz9, y2z9, xz10, yz10, z11;
    float x12, x11y, x10y2, x9y3, x8y4, x7y5, x6y6, x5y7, x4y8, x3y9, x2y10, xy11, y12, x11z, x10yz,
        x9y2z, x8y3z, x7y4z, x6y5z, x5y6z, x4y7z, x3y8z, x2y9z, xy10z, y11z, x10z2, x9yz2, x8y2z2,
        x7y3z2, x6y4z2, x5y5z2, x4y6z2, x3y7z2, x2y8z2, xy9z2, y10z2, x9z3, x8yz3, x7y2z3, x6y3z3,
        x5y4z3, x4y5z3, x3y6z3, x2y7z3, xy8z3, y9z3, x8z4, x7yz4, x6y2z4, x5y3z4, x4y4z4, x3y5z4,
        x2y6z4, xy7z4, y8z4, x7z5, x6yz5, x5y2z5, x4y3z5, x3y4z5, x2y5z5, xy6z5, y7z5, x6z6, x5yz6,
        x4y2z6, x3y3z6, x2y4z6, xy5z6, y6z6, x5z7, x4yz7, x3y2z7, x2y3z7, xy4z7, y5z7, x4z8, x3yz8,
        x2y2z8, xy3z8, y4z8, x3z9, x2yz9, xy2z9, y3z9, x2z10, xyz10, y2z10, xz11, yz11, z12;
} cartesian_moments;

typedef struct {
    double pos[3];
    double m;
    double x, y, z;
    double x2, xy, y2, xz, yz, z2;
    double x3, x2y, xy2, y3, x2z, xyz, y2z, xz2, yz2, z3;
    double x4, x3y, x2y2, xy3, y4, x3z, x2yz, xy2z, y3z, x2z2, xyz2, y2z2, xz3, yz3, z4;
    double x5, x4y, x3y2, x2y3, xy4, y5, x4z, x3yz, x2y2z, xy3z, y4z, x3z2, x2yz2, xy2z2, y3z2,
        x2z3, xyz3, y2z3, xz4, yz4, z5;
    double x6, x5y, x4y2, x3y3, x2y4, xy5, y6, x5z, x4yz, x3y2z, x2y3z, xy4z, y5z, x4z2, x3yz2,
        x2y2z2, xy3z2, y4z2, x3z3, x2yz3, xy2z3, y3z3, x2z4, xyz4, y2z4, xz5, yz5, z6;
    double x7, x6y, x5y2, x4y3, x3y4, x2y5, xy6, y7, x6z, x5yz, x4y2z, x3y3z, x2y4z, xy5z, y6z,
        x5z2, x4yz2, x3y2z2, x2y3z2, xy4z2, y5z2, x4z3, x3yz3, x2y2z3, xy3z3, y4z3, x3z4, x2yz4,
        xy2z4, y3z4, x2z5, xyz5, y2z5, xz6, yz6, z7;
    double x8, x7y, x6y2, x5y3, x4y4, x3y5, x2y6, xy7, y8, x7z, x6yz, x5y2z, x4y3z, x3y4z, x2y5z,
        xy6z, y7z, x6z2, x5yz2, x4y2z2, x3y3z2, x2y4z2, xy5z2, y6z2, x5z3, x4yz3, x3y2z3, x2y3z3,
        xy4z3, y5z3, x4z4, x3yz4, x2y2z4, xy3z4, y4z4, x3z5, x2yz5, xy2z5, y3z5, x2z6, xyz6, y2z6,
        xz7, yz7, z8;
    double x9, x8y, x7y2, x6y3, x5y4, x4y5, x3y6, x2y7, xy8, y9, x8z, x7yz, x6y2z, x5y3z, x4y4z,
        x3y5z, x2y6z, xy7z, y8z, x7z2, x6yz2, x5y2z2, x4y3z2, x3y4z2, x2y5z2, xy6z2, y7z2, x6z3,
        x5yz3, x4y2z3, x3y3z3, x2y4z3, xy5z3, y6z3, x5z4, x4yz4, x3y2z4, x2y3z4, xy4z4, y5z4, x4z5,
        x3yz5, x2y2z5, xy3z5, y4z5, x3z6, x2yz6, xy2z6, y3z6, x2z7, xyz7, y2z7, xz8, yz8, z9;
    double x10, x9y, x8y2, x7y3, x6y4, x5y5, x4y6, x3y7, x2y8, xy9, y10, x9z, x8yz, x7y2z, x6y3z,
        x5y4z, x4y5z, x3y6z, x2y7z, xy8z, y9z, x8z2, x7yz2, x6y2z2, x5y3z2, x4y4z2, x3y5z2, x2y6z2,
        xy7z2, y8z2, x7z3, x6yz3, x5y2z3, x4y3z3, x3y4z3, x2y5z3, xy6z3, y7z3, x6z4, x5yz4, x4y2z4,
        x3y3z4, x2y4z4, xy5z4, y6z4, x5z5, x4yz5, x3y2z5, x2y3z5, xy4z5, y5z5, x4z6, x3yz6, x2y2z6,
        xy3z6, y4z6, x3z7, x2yz7, xy2z7, y3z7, x2z8, xyz8, y2z8, xz9, yz9, z10;
    double x11, x10y, x9y2, x8y3, x7y4, x6y5, x5y6, x4y7, x3y8, x2y9, xy10, y11, x10z, x9yz, x8y2z,
        x7y3z, x6y4z, x5y5z, x4y6z, x3y7z, x2y8z, xy9z, y10z, x9z2, x8yz2, x7y2z2, x6y3z2, x5y4z2,
        x4y5z2, x3y6z2, x2y7z2, xy8z2, y9z2, x8z3, x7yz3, x6y2z3, x5y3z3, x4y4z3, x3y5z3, x2y6z3,
        xy7z3, y8z3, x7z4, x6yz4, x5y2z4, x4y3z4, x3y4z4, x2y5z4, xy6z4, y7z4, x6z5, x5yz5, x4y2z5,
        x3y3z5, x2y4z5, xy5z5, y6z5, x5z6, x4yz6, x3y2z6, x2y3z6, xy4z6, y5z6, x4z7, x3yz7, x2y2z7,
        xy3z7, y4z7, x3z8, x2yz8, xy2z8, y3z8, x2z9, xyz9, y2z9, xz10, yz10, z11;
    double x12, x11y, x10y2, x9y3, x8y4, x7y5, x6y6, x5y7, x4y8, x3y9, x2y10, xy11, y12, x11z,
        x10yz, x9y2z, x8y3z, x7y4z, x6y5z, x5y6z, x4y7z, x3y8z, x2y9z, xy10z, y11z, x10z2, x9yz2,
        x8y2z2, x7y3z2, x6y4z2, x5y5z2, x4y6z2, x3y7z2, x2y8z2, xy9z2, y10z2, x9z3, x8yz3, x7y2z3,
        x6y3z3, x5y4z3, x4y5z3, x3y6z3, x2y7z3, xy8z3, y9z3, x8z4, x7yz4, x6y2z4, x5y3z4, x4y4z4,
        x3y5z4, x2y6z4, xy7z4, y8z4, x7z5, x6yz5, x5y2z5, x4y3z5, x3y4z5, x2y5z5, xy6z5, y7z5, x6z6,
        x5yz6, x4y2z6, x3y3z6, x2y4z6, xy5z6, y6z6, x5z7, x4yz7, x3y2z7, x2y3z7, xy4z7, y5z7, x4z8,
        x3yz8, x2y2z8, xy3z8, y4z8, x3z9, x2yz9, xy2z9, y3z9, x2z10, xyz10, y2z10, xz11, yz11, z12;
} cartesian_moments_d;
