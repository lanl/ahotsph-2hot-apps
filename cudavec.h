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

/* map float4 to accp */
#define Ax .x
#define Ay .y
#define Az .z
#define Phi .w

#define Mass .x
#define X .y
#define Y .z
#define Z .w

#if VECWIDTH == 2

#define VS(a, b) do { \
      a[0] b; \
      a[1] b; \
    } while (0)

#define VV(a, b) do { \
      a[0] b[0]; \
      a[1] b[1]; \
    } while (0)

#define VVS(a, b, c) do { \
      a[0] b[0] c; \
      a[1] b[1] c; \
    } while (0)


#define VVV(a, b, c) do { \
      a[0] b[0] c[0]; \
      a[1] b[1] c[1]; \
    } while (0)

#define VVVV(a, b, c, d) do { \
      a[0] b[0] c[0] c[0];   \
      a[1] b[1] c[1] d[1];   \
    } while (0)

#define Vdecl(a) {a[4*VECWIDTH*index], a[4*VECWIDTH*index+VECWIDTH]}
#define V4decl(a) {a[index], a[index+1]}

#endif

#if VECWIDTH == 3

#define VS(a, b) do { \
      a[0] b; \
      a[1] b; \
      a[2] b; \
    } while (0)

#define VV(a, b) do { \
      a[0] b[0]; \
      a[1] b[1]; \
      a[2] b[2]; \
    } while (0)

#define VVS(a, b, c) do { \
      a[0] b[0] c; \
      a[1] b[1] c; \
      a[2] b[2] c; \
    } while (0)


#define VVV(a, b, c) do { \
      a[0] b[0] c[0]; \
      a[1] b[1] c[1]; \
      a[2] b[2] c[2]; \
    } while (0)

#define VVVV(a, b, c, d) do {	\
      a[0] b[0] c[0] d[0]; \
      a[1] b[1] c[1] d[1]; \
      a[2] b[2] c[2] d[2]; \
    } while (0)

#define Vdecl(a) {a[4*VECWIDTH*index], a[4*VECWIDTH*index+VECWIDTH], a[4*VECWIDTH*index+2*VECWIDTH]}
#define V4decl(a) {a[index], a[index+1], a[index+2]}

#endif


#if VECWIDTH == 4

#define VS(a, b) do { \
      a[0] b; \
      a[1] b; \
      a[2] b; \
      a[3] b; \
    } while (0)

#define VV(a, b) do { \
      a[0] b[0]; \
      a[1] b[1]; \
      a[2] b[2]; \
      a[3] b[3]; \
    } while (0)

#define VVS(a, b, c) do { \
      a[0] b[0] c; \
      a[1] b[1] c; \
      a[2] b[2] c; \
      a[3] b[3] c; \
    } while (0)


#define VVV(a, b, c) do { \
      a[0] b[0] c[0]; \
      a[1] b[1] c[1]; \
      a[2] b[2] c[2]; \
      a[3] b[3] c[3]; \
    } while (0)

#define VVVV(a, b, c, d) do {	\
      a[0] b[0] c[0] d[0]; \
      a[1] b[1] c[1] d[1]; \
      a[2] b[2] c[2] d[2]; \
      a[3] b[3] c[3] d[3]; \
    } while (0)

#define Vdecl(a, s) { \
    {a[(s)*(index*VECWIDTH+0)+0], a[(s)*(index*VECWIDTH+0)+1], a[(s)*(index*VECWIDTH+0)+2], a[(s)*(index*VECWIDTH+0)+3]}, \
    {a[(s)*(index*VECWIDTH+1)+0], a[(s)*(index*VECWIDTH+1)+1], a[(s)*(index*VECWIDTH+1)+2], a[(s)*(index*VECWIDTH+1)+3]}, \
    {a[(s)*(index*VECWIDTH+2)+0], a[(s)*(index*VECWIDTH+2)+1], a[(s)*(index*VECWIDTH+2)+2], a[(s)*(index*VECWIDTH+1)+3]}, \
    {a[(s)*(index*VECWIDTH+3)+0], a[(s)*(index*VECWIDTH+3)+1], a[(s)*(index*VECWIDTH+3)+2], a[(s)*(index*VECWIDTH+1)+3]}}

#define V4decl(a) {a[index], a[index+1], a[index+2], a[index+3]}

#endif
