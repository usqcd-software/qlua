/* create array and create fields; return NULL on fail */
QT(QDP_D_) **QAx(qlm_vector_alloc_,_L)(QDP_Lattice *L, int n)
{
    int status = 0;
    QT(QDP_D_) **res = malloc(sizeof(QT(QDP_D_) *) * n);
    if (NULL == res) return res;
    for (int i = 0 ; i < n ; i++) 
        res[i] = NULL;
    for (int i = 0 ; i < n ; i++) {
        res[i] = QAx(QDP_D_create_,_L)(L);
        if (NULL == res[i])
            status = 1;
    }
    if (status) {
        QA(qlm_vector_free_)(res, n);   /* sic! free checks for NULLs */
        return NULL;
    } else
        return res;
}
/* free array elements and free the array */
void QA(qlm_vector_free_)(QT(QDP_D_) **res, int n) 
{
    if (NULL == res) return;
    for (int i = 0 ; i < n ; i++) 
        if (NULL != res[i])
            QA(QDP_D_destroy_)(res[i]);
    free(res);
}
/* create QLA* array and expose fields; return NULL on fail */
QT(QLA_D_) **QA(qlm_vector_expose_)(QT(QDP_D_) **qdp_x, int n)
{
    int status = 0;
    QT(QLA_D_) **qla_x = malloc(sizeof(QT(QDP_D_) *) * n);
    if (NULL == qla_x) return NULL;
    for (int i = 0 ; i < n ; i++) 
        qla_x[i] = NULL;
    for (int i = 0 ; i < n ; i++) {
        qla_x[i] = QA(QDP_D_expose_)(qdp_x[i]);
        if (NULL == qla_x[i])
            status = 1;
    }
    if (status) {
        QA(qlm_vector_reset_)(qdp_x, n);
        free(qla_x);
        return NULL;
    } else
        return qla_x;
}
void QA(qlm_vector_reset_)(QT(QDP_D_) **qdp_x, int n) 
{
    if (NULL == qdp_x) return ;
    for (int i = 0 ; i < n ; i++) {
        if (NULL != qdp_x[i])
            QA(QDP_D_reset_)(qdp_x[i]);
    }    
}
#undef QT
#undef QA
#undef QAx
