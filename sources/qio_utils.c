#include "qlua.h"                                                    /* DEPS */
#include "lattice.h"                                                 /* DEPS */
#include "qio_utils.h"                                               /* DEPS */
#include "qmp.h"

#include <string.h>

/* This works only because all users of qlua_qio_std do not return open files to
 * QLUA. If they ever do, changes are required!
 */

//static QDP_Lattice *q_lat = NULL; 

static int
x_node_number_a(const int coord[], void *arg)
{
    return QDP_node_number_L((QDP_Lattice *)arg, coord);
}

static int
x_index_a(const int coord[], void *arg)
{
    return QDP_index_L((QDP_Lattice *)arg, coord);
}

static void
x_get_coords_a(int x[], int node, int index, void *arg)
{
    QDP_get_coords_L((QDP_Lattice *)arg, x, node, index);
}

static int
x_numsites_a(int node, void *arg)
{
    return QDP_numsites_L((QDP_Lattice *)arg, node);
}

static void
init_layout(QIO_Layout *layout, mLattice *S)
{
//    q_lat = S->lat;
    layout->arg             = S->lat;
    layout->node_number_a   = &x_node_number_a;
    layout->node_index_a    = &x_index_a;
    layout->get_coords_a    = &x_get_coords_a;
    layout->num_sites_a     = &x_numsites_a;
    layout->latsize         = S->dim;
    layout->latdim          = S->rank;
    layout->volume          = QDP_volume_L(S->lat);
    layout->sites_on_node   = QDP_sites_on_node_L(S->lat);
    layout->this_node       = QDP_this_node;
    layout->number_of_nodes = QMP_get_number_of_nodes();
}

QIO_Reader *
qlua_qio_std_reader(mLattice *S, const char *fname, QIO_String *file_xml)
{
    QIO_Layout       layout;
    QIO_Iflag        iflag;
    QIO_Filesystem   fs;

    init_layout(&layout, S);
    iflag.serpar       = QIO_SERIAL;
    iflag.volfmt       = QIO_SINGLEFILE;
    fs.my_io_node_a    = NULL;
    fs.master_io_node_a= NULL;
    fs.arg             = NULL;
    QIO_string_set(file_xml, "");

    return QIO_open_read(file_xml, fname, &layout, &fs, &iflag);
}

QIO_Writer *
qlua_qio_std_writer(mLattice *S, const char *fname, QIO_String *file_xml,
                    int volfmt)
{
    QIO_Layout       layout;
    QIO_Oflag        oflag;
    QIO_Filesystem   fs;

    init_layout(&layout, S);
    oflag.serpar       = QIO_SERIAL;
    oflag.mode         = QIO_TRUNC;
    oflag.ildgstyle    = QIO_ILDGNO;
    oflag.ildgLFN      = NULL;
    fs.my_io_node_a    = NULL;
    fs.master_io_node_a= NULL;
    fs.arg             = NULL;

    return QIO_open_write(file_xml, fname, volfmt, &layout, &fs, &oflag);
}

char
qlua_qio_file_precision(lua_State *L, int idx)
{
    static const struct {
        char *name;
        int value;
    } fmts[] = {
        { "double", 'D' },
        { "D",      'D' },
        { "float",  'F' },
        { "F",      'F' },
        { NULL,     -1  }
    };
    int i;
    const char *prec = luaL_checkstring(L, idx);

    for (i = 0; fmts[i].name; i++) {
        if (strcmp(fmts[i].name, prec) == 0)
            return fmts[i].value;
    }
    return luaL_error(L, "unsupported file format specification");
}

int 
qlua_qio_volume_format(lua_State *L, int idx, int tmp_idx)
{
    static const struct {
        char *name;
        int value;
    } fmts[] = {
        { "single", QDP_SINGLEFILE },
        { "multi",  QDP_MULTIFILE  },

        /* FIXME these two formats may break compatibility with QDP I/O
           which does not inherit PARTFILE from QIO; they, however, are
           expected to work fine with QIO-based functions like DD_PAIRS */
        { "part",   QIO_PARTFILE   }, 
        { "part_dir",   QIO_PARTFILE_DIR },

        { NULL,     -1  }
    };

    if (idx <= tmp_idx) {
        switch (lua_type(L, idx)) {
        case LUA_TNONE:
        case LUA_TNIL:
            return QDP_SINGLEFILE;
        default: {
            int i;
            const char *volfmt = luaL_checkstring(L, idx);
            for (i = 0; fmts[i].name; i++) {
                if (strcmp(fmts[i].name, volfmt) == 0)
                    return fmts[i].value;
            }
            return luaL_error(L, "unsupported volume format");
        }
        }
    } else {
        return QDP_SINGLEFILE;
    }
}
