/* OpenMP controls */

/*    ICV                     ENV                         C access
 *
 *    dyn                     OMP_DYNAMIC                 bool         r/w *_dynamic()
 *    nested   	              OMP_NESTED                  bool         r/w *_nested()
 *    num_threads	      OMP_NUM_THREADS             int          r/w *_{max,num}_threads()
 *    schedule		      OMP_SCHEDULE                string,int   r/w *_schedule()
 *    bind		      OMP_PROC_BIND               string       r/  *_proc_bind()
 *    thread+limit	      OMP_THREAD_LIMIT            int          r/  *_thread_limit()
 *    max_active_levels	      OMP_MAX_ACTIVE_LEVELS       int          r/w *_max_active_levels()
 *    default-device	      OMP_DEFAULT_DEVICE          int          r/w *_default_device()
 *    num_devices             -                           int          r/  *_num_devices()
 *    version                 -                           string       r/  _OPENMP
 *    
 */
