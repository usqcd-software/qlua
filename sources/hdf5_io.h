#ifndef MARK_FDF0D51E_567F_4C9B_834E_2C115498E41D
#define MARK_FDF0D51E_567F_4C9B_834E_2C115498E41D

typedef struct mHdf5Reader_s mHdf5Reader;
typedef struct mHdf5Writer_s mHdf5Writer;

int init_hdf5_io(lua_State *L);
int fini_hdf5_io(lua_State *L);

mHdf5Reader *qlua_checkHDF5Reader(lua_State *L, int idx);
mHdf5Writer *qlua_checkHDF5Writer(lua_State *L, int idx);

#endif /* defined(MARK_FDF0D51E_567F_4C9B_834E_2C115498E41D) */
