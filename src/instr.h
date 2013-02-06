#ifdef INSTR_WITH_SCT

#include "sct.h"
#define INSTR_INIT(N,CS) sct_init(N,CS)
#define INSTR_DEF(T,S) static struct { int it; const char *s; } T = { -1, S }
#define INSTR_SET(T,V) if(!T) T = sct_new_timer(T,V)
#define INSTR_START(T) if(T.it<0) T.it = sct_new_timer(T.s); sct_start(T.it)
#define INSTR_STOP(T) sct_stop(T.it)
#define INSTR_REPORT(A,B,C) sct_report(A,B,C)
#define INSTR_FINALIZE(...)

#else

#define INSTR_INIT(...)
#define INSTR_DEF(...)
#define INSTR_SET(...)
#define INSTR_START(...)
#define INSTR_STOP(...)
#define INSTR_FINALIZE(...)

#endif
