#ifndef MARK_D62788EF_1DAC_49A9_AC44_674146CFF599
#define MARK_D62788EF_1DAC_49A9_AC44_674146CFF599

#define QLUA_MAX_LATTICE_RANK 64 /* wired to make layout routine alloc-free */

void qlua_sublattice(int lo[], int hi[], int node, void *env);
int qlua_comm_map(const int *net_coord, void *env);
extern QDP_Layout qlua_layout;

#endif /* !defined(MARK_D62788EF_1DAC_49A9_AC44_674146CFF599) */
