typedef struct {
    float pos[3];          /* particle position */
    float strength[3];     /* particle strength = omega X vol */
    float vol;             /* particle vol=(volume/4*pi) */
    PAD_DECL
    ACCUM psi[3];          /* particle streamfunction = psi */
    ACCUM vel[3];          /* particle velocity = curl(psi) */
    ACCUM gradvel[3][3];   /* particle grad(vel)            */
    ACCUM dstr[3];         /* particle "stretching + viscous interaction" */ 
    ACCUM omegat[3];  	/* part of gsw's relaxation scheme*/
    float vel_old[3];      /* old particle velocity */
    float dstr_old[3];     /* old particle "stretching" + visc. interaction */ 
    float errsum;		/* sum of a posteriori error bounds */
    float errsum2;		/* sum of squares of a.p. error bounds */
    unsigned int ident;    /* particle number */
    Key_t key;			/* hash key */
    int nterms;	           /* number of terms in field eval */
} body, *bodyptr;

/* Rely on ANSI-style string-concatenation.  */
#define WHOLEBODYDESC \
"struct {\n\
    float x,y,z;           /* particle position */\n\
    float strx,stry,strz;  /* particle strength = omega X vol */\n\
    float vol;	        	/* volume/4pi */\n\
    " PAD_DECL_S "		/* padding for aligned doubles to follow */\n\
    " ACCUM_S " psix, psiy, psiz;     /* particle streamfunction = psi */\n\
    " ACCUM_S " vx, vy, vz;           /* particle velocity = curl(psi) */\n\
    " ACCUM_S " dvxdx, dvxdy, dvxdz,\n\
	  dvydx, dvydy, dvydz,\n\
	  dvzdx, dvzdy, dvzdz; /* grad vel */\n\
    " ACCUM_S " dstrx, dstry, dstrz; /* stretching + viscous interaction */\n\
    float vx_old, vy_old, vz_old; /* old particle velocity */\n\
    float dstrx_old, dstry_old, dstrz_old;  /* old value of dstr */\n\
    float errsum;        /* sum of a posteriori error bounds */\n\
    float errsum2;       /* sum of squares of " " */\n\
    unsigned int ident;         /* particle number */\n\
    unsigned long key[2];	/* hash key */\n\
    int nterms;			/* number of terms in field eval */\n\
}"

/* Use this when just dumping positions and strengths */
typedef struct{
    float pos[3];
    float strength[3];
    float vol;
    int ident;
} outbody;

#define OUTBODYDESC \
"struct {\n\
  float x, y, z;\n\
  float strx, stry, strz;\n\
  float vol;\n\
  int ident;\n\
}"

struct dmoment{           /* dipole terms of multipole expansion */
    float x[3];              
    float y[3];
    float z[3];
};

struct qmoment{           /* quadrupole terms of multipole expansion */
    float xx[3];             
    float yy[3];
    float zz[3];
    float xy[3];
    float xz[3];
    float yz[3];
};

typedef struct {
    float pos[3];             /* expansion center (= geometric cell center)*/
    float strength[3];        /* monopole term */
    struct dmoment dpole;     /* dipole terms */
    struct qmoment qpole;     /* quadripole terms */
    float bmax;	              /* length of half-diagonal of cell */
    float rcrit2;		/* closest acceptable distance */
    float b3, b4;		/* parameters of the error bound */
    int daughters;            /* how many particles in that cell */
} cell, *cellptr;

typedef struct {
    body *bp;
} Sink;

typedef struct{
    float strength[NDIM];
    float pos[NDIM];
    struct dmoment dpole;
    struct qmoment qpole;
    float b0;
    float b2;
    float bmax;
    int daughters;
} cofm_data;

