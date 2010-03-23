#include "qlua.h"                                                    /* DEPS */
#include "lattice.h"                                                 /* DEPS */
#include "qio_utils.h"                                               /* DEPS */
#include "qmp.h"

static void
init_layout(QIO_Layout *layout, mLattice *S)
{
    layout->node_number     = &QDP_node_number;
    layout->node_index      = &QDP_index;
    layout->get_coords      = &QDP_get_coords;
    layout->num_sites       = &QDP_numsites;
    layout->latsize         = S->dim;
    layout->latdim          = S->rank;
    layout->volume          = QDP_volume();
    layout->sites_on_node   = QDP_sites_on_node;
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
    fs.my_io_node      = NULL;
    fs.master_io_node  = NULL;
    QIO_string_set(file_xml, "");

    return QIO_open_read(file_xml, fname, &layout, &fs, &iflag);
}

QIO_Writer *
qlua_qio_std_writer(mLattice *S, const char *fname, QIO_String *file_xml)
{
    QIO_Layout       layout;
    QIO_Oflag        oflag;
    QIO_Filesystem   fs;

    init_layout(&layout, S);
    oflag.serpar       = QIO_SERIAL;
    oflag.mode         = QIO_TRUNC;
    oflag.ildgstyle    = QIO_ILDGNO;
    oflag.ildgLFN      = NULL;
    fs.my_io_node      = NULL;
    fs.master_io_node  = NULL;

    return QIO_open_write(file_xml, fname, QIO_SINGLEFILE,
                          &layout, &fs, &oflag);
}
