/*! $Id: .vimrc,v 1.9 2007/08/29 18:08:37 syritsyn Exp $
 **************************************************************************
 * \file debug.h
 *  debugging primitives
 *  Define three levels of debug output
 *      0) always verbose
 *      1) function calls, enter/leave
        2) 
 *      3) detailed code-string based output on (almost) each separate action, 
 *         speculative problem anticipating output
 *      4) debug print from heavily called funargs; desperate measure
 *  TODO think of some log file separate from stderr
 *       meaning of debug level 2?
 *
 * \author Sergey N. Syritsyn
 * 
 * \date Created: 11/11/2007
 *
 ***************************************************************************/
#ifndef QDPW_DEBUG_H__ 
#define QDPW_DEBUG_H__ 

#include <stdio.h>

#ifndef DEBUG_LEVEL
# define DEBUG_LEVEL 0
#endif

#ifndef DEBUG_ECHO  // for external redefinition
# ifdef QDP_Precision // qdp present
#  define DEBUG_ECHO(...) {\
    fprintf(stderr,"[%05d]%s: ", QDP_this_node, __func__);\
    fprintf(stderr,__VA_ARGS__);\
    fflush(stderr);\
}
# else
#  define DEBUG_ECHO(...) {\
    fprintf(stderr,"%s: ",__func__);\
    fprintf(stderr,__VA_ARGS__);\
    fflush(stderr);\
}
# endif
#endif

#define DEBUG_MUTE

#if DEBUG_LEVEL>4
# warning "DEBUG level is above the highest supported"
#endif

#if DEBUG_LEVEL>3
//# warning "DEBUG level 4"
# define DEBUG_L1(...) DEBUG_ECHO(__VA_ARGS__)
# define DEBUG_L2(...) DEBUG_ECHO(__VA_ARGS__)
# define DEBUG_L3(...) DEBUG_ECHO(__VA_ARGS__)
# define DEBUG_L4(...) DEBUG_ECHO(__VA_ARGS__)
#elif DEBUG_LEVEL>2
//# warning "DEBUG level 3"
# define DEBUG_L1(...) DEBUG_ECHO(__VA_ARGS__)
# define DEBUG_L2(...) DEBUG_ECHO(__VA_ARGS__)
# define DEBUG_L3(...) DEBUG_ECHO(__VA_ARGS__)
# define DEBUG_L4(...) DEBUG_MUTE
#elif DEBUG_LEVEL>1
//# warning "DEBUG level 2"
# define DEBUG_L1(...) DEBUG_ECHO(__VA_ARGS__)
# define DEBUG_L2(...) DEBUG_ECHO(__VA_ARGS__)
# define DEBUG_L3(...) DEBUG_MUTE
# define DEBUG_L4(...) DEBUG_MUTE
#elif DEBUG_LEVEL>0
//# warning "DEBUG level 1"
# define DEBUG_L1(...) DEBUG_ECHO(__VA_ARGS__)
# define DEBUG_L2(...) DEBUG_MUTE
# define DEBUG_L3(...) DEBUG_MUTE
# define DEBUG_L4(...) DEBUG_MUTE
#else
//# warning "DEBUG level 0"
# define DEBUG_L1(...) DEBUG_MUTE
# define DEBUG_L2(...) DEBUG_MUTE
# define DEBUG_L3(...) DEBUG_MUTE
# define DEBUG_L4(...) DEBUG_MUTE
#endif

#define DEBUG_L0(...) DEBUG_ECHO(__VA_ARGS__)

#define DEBUG_ENTER DEBUG_L1("enter\n")
#define DEBUG_LEAVE DEBUG_L1("leave\n")
#define DEBUG_IN    DEBUG_ENTER
#define DEBUG_OUT   DEBUG_LEAVE

#define LOG_ERROR   DEBUG_ECHO

#define LOG_PRINT(...)      { if (!QDP_is_initialized() || 0 == QDP_this_node) printf(__VA_ARGS__); }
#define LOG_FPRINT(fp, ...) { if (!QDP_is_initialized() || 0 == QDP_this_node) fprintf(fp, __VA_ARGS__); }
#define LOG_ECHO(fp,...)    { printf("%s[%d]: ", __func__, QDP_this_node); printf(__VA_ARGS__); }
#define LOG_FECHO(fp,...)   { fprintf(fp, "%s[%d]: ", __func__, QDP_this_node); fprintf(fp, __VA_ARGS__); }

#endif/*QDPW_DEBUG_H__ */
