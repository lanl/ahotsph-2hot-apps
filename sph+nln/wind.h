
void
SDFwritewind(const char *filename, int gnobj, int nobj, 
	     const void *btab, int windnobj, const void *windbtab, 
	     int bsize, int wsize, const char *winddesc, 
	     const char *bodydesc, ...);

SDF *SDFreadwind(char *name, void **btabp, int *gnobjp, int *nobjp, 
		 int stride,
		 /* char *name, offset_t offset, int *confirm */...);
