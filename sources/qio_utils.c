#include <qlua.h>                                                    /* DEPS */
#include <lattice.h>                                                 /* DEPS */
#include <qio_utils.h>                                               /* DEPS */
#include <qmp.h>

static void
init_layout(QIO_Layout *layout)
{
    layout->node_number     = &QDP_node_number;
    layout->node_index      = &QDP_index;
    layout->get_coords      = &QDP_get_coords;
    layout->num_sites       = &QDP_numsites;
    layout->latsize         = qDim;
    layout->latdim          = qRank;
    layout->volume          = QDP_volume();
    layout->sites_on_node   = QDP_sites_on_node;
    layout->this_node       = QDP_this_node;
    layout->number_of_nodes = QMP_get_number_of_nodes();
}

QIO_Reader *
qlua_qio_std_reader(const char *fname, QIO_String *file_xml)
{
    QIO_Layout       layout;
    QIO_Iflag        iflag;
    QIO_Filesystem   fs;

    init_layout(&layout);
    iflag.serpar = QIO_SERIAL;
    iflag.volfmt = QIO_SINGLEFILE;
    fs.my_io_node = NULL;
    fs.master_io_node = NULL;
    QIO_string_set(file_xml, "");

    return QIO_open_read(file_xml, fname, &layout, &fs, &iflag);
}
