#ifndef MARK_FDF0D51E_567F_4C9B_834E_2C115498E41D
#define MARK_FDF0D51E_567F_4C9B_834E_2C115498E41D

typedef struct mHdf5File_s mHdf5File;

int init_hdf5_io(lua_State *L);
int fini_hdf5_io(lua_State *L);

mHdf5File *qlua_checkHDF5Reader(lua_State *L, int idx);
mHdf5File *qlua_checkHDF5Writer(lua_State *L, int idx);

#endif /* defined(MARK_FDF0D51E_567F_4C9B_834E_2C115498E41D) */
