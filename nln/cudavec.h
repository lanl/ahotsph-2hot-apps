/*
 * © 2026. Triad National Security, LLC. All rights reserved.
 * This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare. derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.
 */

/* UGLY, but you can use these for parentheses inside the VV macros! */
#define LPAREN (
#define RPAREN )
#define COMMA ,
#define Dot(a, b) (VVinfix(a, *b, +))
#define VVinfix(a, b, op) a[0] b[0] op a[1] b[1] op a[2] b[2]

#define LPAREN (
#define RPAREN )

#define VECWIDTH 4
#define NSSE 8

/* map float4 to acc */
#define Ax x
#define Ay y
#define Az z
#define Phi w

#if VECWIDTH == 1

#define VS(a, b) \
    do { a[0] b; } while (0)

#define VV(a, b) \
    do { a[0] b[0]; } while (0)

#define VVS(a, b, c) \
    do { a[0] b[0] c; } while (0)

#define VVV(a, b, c) \
    do { a[0] b[0] c[0]; } while (0)

#define VVVV(a, b, c, d) \
    do { a[0] b[0] c[0] d[0]; } while (0)

// clang-format off
#define Vdecl(a, s)                          \
    {                                        \
        a[(s) * (index * VECWIDTH + 0) + 0], \
        a[(s) * (index * VECWIDTH + 0) + 1], \
        a[(s) * (index * VECWIDTH + 0) + 2]  \
    }
// clang-format on

#define V4decl(a) \
    { a[index] }

#endif

#if VECWIDTH == 2

#define VS(a, b) \
    do {         \
        a[0] b;  \
        a[1] b;  \
    } while (0)

#define VV(a, b)   \
    do {           \
        a[0] b[0]; \
        a[1] b[1]; \
    } while (0)

#define VVS(a, b, c) \
    do {             \
        a[0] b[0] c; \
        a[1] b[1] c; \
    } while (0)

#define VVV(a, b, c)    \
    do {                \
        a[0] b[0] c[0]; \
        a[1] b[1] c[1]; \
    } while (0)

#define VVVV(a, b, c, d)     \
    do {                     \
        a[0] b[0] c[0] d[0]; \
        a[1] b[1] c[1] d[1]; \
    } while (0)

// clang-format off
#define Vdecl(a, s)                              \
    {                                            \
        {                                        \
            a[(s) * (index * VECWIDTH + 0) + 0], \
            a[(s) * (index * VECWIDTH + 0) + 1], \
            a[(s) * (index * VECWIDTH + 0) + 2]  \
        },                                       \
        {                                        \
            a[(s) * (index * VECWIDTH + 1) + 0], \
            a[(s) * (index * VECWIDTH + 1) + 1], \
            a[(s) * (index * VECWIDTH + 1) + 2]  \
        }                                        \
    }
// clang-format on

#define V4decl(a) \
    { a[index], a[index + 1] }

#endif

#if VECWIDTH == 3

#define VS(a, b) \
    do {         \
        a[0] b;  \
        a[1] b;  \
        a[2] b;  \
    } while (0)

#define VV(a, b)   \
    do {           \
        a[0] b[0]; \
        a[1] b[1]; \
        a[2] b[2]; \
    } while (0)

#define VVS(a, b, c) \
    do {             \
        a[0] b[0] c; \
        a[1] b[1] c; \
        a[2] b[2] c; \
    } while (0)

#define VVV(a, b, c)    \
    do {                \
        a[0] b[0] c[0]; \
        a[1] b[1] c[1]; \
        a[2] b[2] c[2]; \
    } while (0)

#define VVVV(a, b, c, d)     \
    do {                     \
        a[0] b[0] c[0] d[0]; \
        a[1] b[1] c[1] d[1]; \
        a[2] b[2] c[2] d[2]; \
    } while (0)

// clang-format off
#define Vdecl(a, s)                              \
    {                                            \
        {                                        \
            a[(s) * (index * VECWIDTH + 0) + 0], \
            a[(s) * (index * VECWIDTH + 0) + 1], \
            a[(s) * (index * VECWIDTH + 0) + 2]  \
        },                                       \
        {                                        \
            a[(s) * (index * VECWIDTH + 1) + 0], \
            a[(s) * (index * VECWIDTH + 1) + 1], \
            a[(s) * (index * VECWIDTH + 1) + 2]},\
        {                                        \
            a[(s) * (index * VECWIDTH + 2) + 0], \
            a[(s) * (index * VECWIDTH + 2) + 1], \
            a[(s) * (index * VECWIDTH + 2) + 2]  \
        }                                        \
    }
// clang-format on

#define V4decl(a) \
    { a[index], a[index + 1], a[index + 2] }

#endif

#if VECWIDTH == 4

#define VS(a, b) \
    do {         \
        a[0] b;  \
        a[1] b;  \
        a[2] b;  \
        a[3] b;  \
    } while (0)

#define VV(a, b)   \
    do {           \
        a[0] b[0]; \
        a[1] b[1]; \
        a[2] b[2]; \
        a[3] b[3]; \
    } while (0)

#define VVS(a, b, c) \
    do {             \
        a[0] b[0] c; \
        a[1] b[1] c; \
        a[2] b[2] c; \
        a[3] b[3] c; \
    } while (0)

#define VVV(a, b, c)    \
    do {                \
        a[0] b[0] c[0]; \
        a[1] b[1] c[1]; \
        a[2] b[2] c[2]; \
        a[3] b[3] c[3]; \
    } while (0)

#define VVVV(a, b, c, d)     \
    do {                     \
        a[0] b[0] c[0] d[0]; \
        a[1] b[1] c[1] d[1]; \
        a[2] b[2] c[2] d[2]; \
        a[3] b[3] c[3] d[3]; \
    } while (0)

// clang-format off
#define Vdecl(a, s)                              \
    {                                            \
        {                                        \
            a[(s) * (index * VECWIDTH + 0) + 0], \
            a[(s) * (index * VECWIDTH + 0) + 1], \
            a[(s) * (index * VECWIDTH + 0) + 2]  \
        },                                       \
        {                                        \
            a[(s) * (index * VECWIDTH + 1) + 0], \
            a[(s) * (index * VECWIDTH + 1) + 1], \
            a[(s) * (index * VECWIDTH + 1) + 2]  \
        },                                       \
        {                                        \
            a[(s) * (index * VECWIDTH + 2) + 0], \
            a[(s) * (index * VECWIDTH + 2) + 1], \
            a[(s) * (index * VECWIDTH + 2) + 2]  \
        },                                       \
        {                                        \
            a[(s) * (index * VECWIDTH + 3) + 0], \
            a[(s) * (index * VECWIDTH + 3) + 1], \
            a[(s) * (index * VECWIDTH + 3) + 2]  \
        }                                        \
    }
// clang-format on

#define V4decl(a) \
    { a[index], a[index + 1], a[index + 2], a[index + 3] }

#endif

#if VECWIDTH == 6

#define VS(a, b) \
    do {         \
        a[0] b;  \
        a[1] b;  \
        a[2] b;  \
        a[3] b;  \
        a[4] b;  \
        a[5] b;  \
    } while (0)

#define VV(a, b)   \
    do {           \
        a[0] b[0]; \
        a[1] b[1]; \
        a[2] b[2]; \
        a[3] b[3]; \
        a[4] b[4]; \
        a[5] b[5]; \
    } while (0)

#define VVS(a, b, c) \
    do {             \
        a[0] b[0] c; \
        a[1] b[1] c; \
        a[2] b[2] c; \
        a[3] b[3] c; \
        a[4] b[4] c; \
        a[5] b[5] c; \
    } while (0)

#define VVV(a, b, c)    \
    do {                \
        a[0] b[0] c[0]; \
        a[1] b[1] c[1]; \
        a[2] b[2] c[2]; \
        a[3] b[3] c[3]; \
        a[4] b[4] c[4]; \
        a[5] b[5] c[5]; \
    } while (0)

#define VVVV(a, b, c, d)     \
    do {                     \
        a[0] b[0] c[0] d[0]; \
        a[1] b[1] c[1] d[1]; \
        a[2] b[2] c[2] d[2]; \
        a[3] b[3] c[3] d[3]; \
        a[4] b[4] c[4] d[4]; \
        a[5] b[5] c[5] d[5]; \
    } while (0)

// clang-format off
#define Vdecl(a, s)                              \
    {                                            \
        {                                        \
            a[(s) * (index * VECWIDTH + 0) + 0], \
            a[(s) * (index * VECWIDTH + 0) + 1], \
            a[(s) * (index * VECWIDTH + 0) + 2]  \
        },                                       \
        {                                        \
            a[(s) * (index * VECWIDTH + 1) + 0], \
            a[(s) * (index * VECWIDTH + 1) + 1], \
            a[(s) * (index * VECWIDTH + 1) + 2]  \
        },                                       \
        {                                        \
            a[(s) * (index * VECWIDTH + 2) + 0], \
            a[(s) * (index * VECWIDTH + 2) + 1], \
            a[(s) * (index * VECWIDTH + 2) + 2]  \
        },                                       \
        {                                        \
            a[(s) * (index * VECWIDTH + 3) + 0], \
            a[(s) * (index * VECWIDTH + 3) + 1], \
            a[(s) * (index * VECWIDTH + 3) + 2]  \
        },                                       \
        {                                        \
            a[(s) * (index * VECWIDTH + 4) + 0], \
            a[(s) * (index * VECWIDTH + 4) + 1], \
            a[(s) * (index * VECWIDTH + 4) + 2]  \
        },                                       \
        {                                        \
            a[(s) * (index * VECWIDTH + 5) + 0], \
            a[(s) * (index * VECWIDTH + 5) + 1], \
            a[(s) * (index * VECWIDTH + 5) + 2]  \
        }                                        \
    }
// clang-format on

#define V4decl(a) \
    { a[index], a[index + 1], a[index + 2], a[index + 3], a[index + 4], a[index + 5] }

#endif

#if VECWIDTH == 8

#define VS(a, b) \
    do {         \
        a[0] b;  \
        a[1] b;  \
        a[2] b;  \
        a[3] b;  \
        a[4] b;  \
        a[5] b;  \
        a[6] b;  \
        a[7] b;  \
    } while (0)

#define VV(a, b)   \
    do {           \
        a[0] b[0]; \
        a[1] b[1]; \
        a[2] b[2]; \
        a[3] b[3]; \
        a[4] b[4]; \
        a[5] b[5]; \
        a[6] b[6]; \
        a[7] b[7]; \
    } while (0)

#define VVS(a, b, c) \
    do {             \
        a[0] b[0] c; \
        a[1] b[1] c; \
        a[2] b[2] c; \
        a[3] b[3] c; \
        a[4] b[4] c; \
        a[5] b[5] c; \
        a[6] b[6] c; \
        a[7] b[7] c; \
    } while (0)

#define VVV(a, b, c)    \
    do {                \
        a[0] b[0] c[0]; \
        a[1] b[1] c[1]; \
        a[2] b[2] c[2]; \
        a[3] b[3] c[3]; \
        a[4] b[4] c[4]; \
        a[5] b[5] c[5]; \
        a[6] b[6] c[6]; \
        a[7] b[7] c[7]; \
    } while (0)

#define VVVV(a, b, c, d)     \
    do {                     \
        a[0] b[0] c[0] d[0]; \
        a[1] b[1] c[1] d[1]; \
        a[2] b[2] c[2] d[2]; \
        a[3] b[3] c[3] d[3]; \
        a[4] b[4] c[4] d[4]; \
        a[5] b[5] c[5] d[5]; \
        a[6] b[6] c[6] d[6]; \
        a[7] b[7] c[7] d[7]; \
    } while (0)

// clang-format off
#define Vdecl(a, s)                              \
    {                                            \
        {                                        \
            a[(s) * (index * VECWIDTH + 0) + 0], \
            a[(s) * (index * VECWIDTH + 0) + 1], \
            a[(s) * (index * VECWIDTH + 0) + 2]  \
            },                                   \
        {                                        \
            a[(s) * (index * VECWIDTH + 1) + 0], \
            a[(s) * (index * VECWIDTH + 1) + 1], \
            a[(s) * (index * VECWIDTH + 1) + 2]  \
            },                                   \
        {                                        \
            a[(s) * (index * VECWIDTH + 2) + 0], \
            a[(s) * (index * VECWIDTH + 2) + 1], \
            a[(s) * (index * VECWIDTH + 2) + 2]  \
            },                                   \
        {                                        \
            a[(s) * (index * VECWIDTH + 3) + 0], \
            a[(s) * (index * VECWIDTH + 3) + 1], \
            a[(s) * (index * VECWIDTH + 3) + 2]  \
            },                                   \
        {                                        \
            a[(s) * (index * VECWIDTH + 4) + 0], \
            a[(s) * (index * VECWIDTH + 4) + 1], \
            a[(s) * (index * VECWIDTH + 4) + 2]  \
            },                                   \
        {                                        \
            a[(s) * (index * VECWIDTH + 5) + 0], \
            a[(s) * (index * VECWIDTH + 5) + 1], \
            a[(s) * (index * VECWIDTH + 5) + 2]  \
            },                                   \
        {                                        \
            a[(s) * (index * VECWIDTH + 6) + 0], \
            a[(s) * (index * VECWIDTH + 6) + 1], \
            a[(s) * (index * VECWIDTH + 6) + 2]  \
            },                                   \
        {                                        \
            a[(s) * (index * VECWIDTH + 7) + 0], \
            a[(s) * (index * VECWIDTH + 7) + 1], \
            a[(s) * (index * VECWIDTH + 7) + 2]  \
        }                                        \
    }
// clang-format on

#define V4decl(a)                                                                       \
    {                                                                                   \
        a[index], a[index + 1], a[index + 2], a[index + 3], a[index + 4], a[index + 5], \
            a[index + 6], a[index + 7]                                                  \
    }

#endif
