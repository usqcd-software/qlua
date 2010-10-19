#ifndef MARK_4D25B115_959E_4173_9A18_C3135910F1F6
#define MARK_4D25B115_959E_4173_9A18_C3135910F1F6

QIO_Reader *qlua_qio_std_reader(mLattice *S,
                                const char *fname,
                                QIO_String *file_xml);
QIO_Writer *qlua_qio_std_writer(mLattice *S,
                                const char *fname,
                                QIO_String *file_xml,
                                int volfmt);
char qlua_qio_file_precision(lua_State *L, int idx);
int qlua_qio_volume_format(lua_State *L, int fmt_idx, int tmp_idx);

#endif /* !defined(MARK_4D25B115_959E_4173_9A18_C3135910F1F6) */
