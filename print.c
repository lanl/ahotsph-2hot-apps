#include <stdio.h>

#include "physics.h"
#include "protos.h"
#include "vop.h"

/* Make this return a ptr to static data so we can finesse the */
/* problem of what kind of FILE *! */
char *PrintCellContents(const cell *p) {
    static char contents_string[256];

    if (p == NULL) {
        sprintf(contents_string, "\t(null)");
    } else {
        sprintf(contents_string,
                "\tmass: %.4g, nd:%ld, bmax:%.2g, rcrit:%.2g\n"
                "\t" Sinfix("%.4f", " "),
                p->mass,
                (long int)p->daughters,
                p->bmax,
                p->rcrit,
                Vinfix(p->pos, COMMA));
    }
    return contents_string;
}

char *PrintCellContents4(const hexacell *p) {
    static char contents_string[256];

    if (p == NULL) {
        sprintf(contents_string, "\t(null)");
    } else {
        sprintf(contents_string,
                "\tmass: %.4g, nd:%ld, bmax:%.2g, rcrit:%.4g %.4g %.4g\n"
                "\t" Sinfix("%.4f", " "),
                p->mass,
                (long int)p->daughters,
                p->bmax,
                p->rcrit_m,
                p->rcrit_q,
                p->rcrit_h,
                Vinfix(p->pos, COMMA));
    }
    return contents_string;
}


char *PrintBodyContents(const body *p) {
    static char contents_string[256];

    if (p == NULL) {
        sprintf(contents_string, "\t(null)");
    } else {
        sprintf(contents_string,
                "\tid:%ld, mass:%.4g\n"
                "\t" Sinfix("%.4f", " "),
                p->ident,
                p->mass,
                Vinfix(p->pos, COMMA));
    }
    return contents_string;
}

/* For out of bits confirmation */
char *PrintBodyContentsLong(const body *p) {
    static char contents_string[256];

    if (p == NULL) {
        sprintf(contents_string, "\t(null)");
    } else {
        sprintf(contents_string,
                "\t@%#lx %7.4g\n"
                "\t" Sinfix("%14.10f", " ") " " Sinfix("%x", " "),
                (unsigned long)p,
                p->mass,
                Vinfix(p->pos, COMMA),
                Vinfix(*(int *)&p->pos, COMMA));
    }
    return contents_string;
}

char *PrintBranch(const cofmdata *cmp) {
    static char ret[512];
    sprintf(ret,
            "Br: mass: %.3g, ndaughters=%ld, pos=(%.3f %.3f %.3f), bmax:%.2g\n",
            cmp->m,
            (long int)cmp->ndaughters,
            cmp->center[0],
            cmp->center[1],
            cmp->center[2],
            cmp->bmax);
    return ret;
}
