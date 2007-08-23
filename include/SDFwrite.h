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

/* An extension of SDFwrite to write wind sources and SPH particles to a
   single SDF file */
void SDFwritewind(const char *filename, int gnobj, int nobj, const void *btab,
		  int windnobj, const void *windbtab, 
		  int bsize, int wsize, const char *winddesc, 
		  const char *bodydesc,
		  /* const char *name, SDF_type_enum type, <type> val */ ...);
#ifdef __cplusplus
}
#endif

#endif
