#ifndef _SDFWriteDOTh
#define _SDFWRiteDOTh

#ifdef __cplusplus
extern "C" {
#endif
void SDFwrite(const char *filename, int gnobj, int nobj, const void *btab,
	    int bsize, const char *bodydesc, 
	      /* const char *name, SDF_type_enum type, <type> val */ ...);

/* A trivial special case, just don't write any body data */
void SDFwritehdr(const char *filename, const char *bodydesc, 
		 /* const char *name, SDF_type_enum type, <type> val */ ...);
#ifdef __cplusplus
}
#endif

#endif
