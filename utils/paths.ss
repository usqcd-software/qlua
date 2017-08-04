;; This Scheme code generates VM code for various gauge paths.
;; The code generated includes action, gauge forces, and clover terms.
;; See the bottom of the file for usage.

;; change axes in the path
(define (change-path mp pl)
  (cond
   [(null? pl) '()]
   [(assq (caar pl) mp)
    => (lambda (o&n) (cons (cons (cadr o&n) (cdar pl)) (change-path mp (cdr pl))))]
   [else (error 'change-path "Unknown direction")]))

;; split list into two pieces of about equal length
(define (split-list v* k) ;; => (k a b)
  (let loop ([a '()] [b '()] [v* v*])
    (cond
     [(null? v*) (k a b)]
     [(null? (cdr v*)) (k (cons (car v*) a) b)]
     [else (loop (cons (car v*) a) (cons (cadr v*) b) (cddr v*))])))

;; double the steps
(define (double-path pl)
  (cond
   [(null? pl) pl]
   [else (append (list (car pl) (car pl)) (double-path (cdr pl)))]))

;; triple the steps
(define (triple-path pl)
  (cond
   [(null? pl) pl]
   [else (append (list (car pl) (car pl) (car pl)) (triple-path (cdr pl)))]))

;; shift path by one step
(define (shift-path pl)
  (append (cdr pl) (list (car pl))))

;; rotate the path to a given starting step
(define (set-start pl m)
  (if (equal? (car pl) m) pl
      (set-start (shift-path pl) m)))
   
;; reverse path
(define (reverse-path pl)
  (define (invert-link lk)
    (list (car lk) (if (eq? (cadr lk) '+) '-  '+)))
  (let loop ([r '()] [p pl])
    (cond
     [(null? p) r]
     [else (loop (cons (invert-link (car p)) r) (cdr p))])))

;; check that the path is closed
(define (closed-path? pl)
  (define (d->o d) (if (eq? d '+) +1 -1))
  (define (track-shift p d*)
    (cond
     [(assq (car p) d*)
      => (lambda (x&s) (set-cdr! x&s (+ (cdr x&s) (d->o (cadr p)))) d*)]
     [else (cons (cons (car p) (d->o (cadr p))) d*)]))
  (define (check-zeros d*)
    (cond
     [(null? d*) #t]
     [(zero? (cdar d*)) (check-zeros (cdr d*))]
     [else #f]))
  (let loop ([d '()] [p pl])
    (cond
     [(null? p) (check-zeros d)]
     [else (loop (track-shift (car p) d) (cdr p))])))

;; count steps m in a path
(define (count-shifts pl m)
  (cond
   [(null? pl) 0]
   [(equal? (car pl) m) (+ 1 (count-shifts (cdr pl) m))]
   [else (count-shifts (cdr pl) m)]))

;; the head of a closed path
(define (path-head pl)
  (let loop ([r '()] [p pl] [n (length pl)])
    (cond
     [(zero? n) (reverse r)]
     [(< n 0) (error 'path-head "path of odd length")]
     [else (loop (cons (car p) r) (cdr p) (- n 2))])))

;; reversed tail of a closed path
(define (path-tail pl)
  (path-head (reverse-path pl)))

;; ops
(define (mk-ref n)                    `(ref ,n))
(define (mk-adj x)                    `(adjoin ,x))
(define (mk-mul x y)                  `(mul ,x ,y))
(define (mk-add x y)                  `(add ,x ,y))
(define (mk-shift-from-forward x d)   `(shift-from-forward ,x ,d))
(define (mk-shift-from-backward x d)  `(shift-from-backward ,x ,d))
(define (mk-dot x y )                 `(dot ,x ,y))
(define (mk-real x)                   `(real ,x))

(define (mk-sum v*)
  (cond
   [(null? v*) (error 'mk-sum "Empty list of values")]
   [(null? (cdr v*)) (car v*)]
   [else (split-list v* (lambda (a b) (mk-add (mk-sum a) (mk-sum b))))]))

(define (op-ref k)       (lambda (r) `(ref (,k) (,r) ())))
(define (op-adj)         (lambda (r s) `(adjoin () (,r) (,s))))
(define (op-real)        (lambda (r s) `(real () (,r) (,s))))
(define (op-add)         (lambda (r a b) `(add () (,r) (,a ,b))))
(define (op-mul)         (lambda (r a b) `(mul () (,r) (,a ,b))))
(define (op-dot)         (lambda (r a b) `(dot () (,r) (,a ,b))))
(define (op-shift k d)   (lambda (r s) `(shift (,k ,d) (,r) (,s))))

(define (op-die r*)      (lambda () `(die ,r* () ())))
(define (op-kill r)      (lambda () `(kill (,r) () ())))
(define (op-live r*)     (lambda () `(live ,r* () ())))

;; lua does not allow more than 200 locals, hence reusing registers
(define (emit-op f mk-idx op mr* fr* nx k) ;; => (k x mr* fr* nx)
  (define (reg n) (format "~a" n))
  (define (get-in r)
    (cond
     [(assoc r mr*) => cdr]
     [else (erorr 'emit-op "Unmapped in reg")]))
  (define (get-out r mr* fr* nx kk) ;; => (kk x mr* fr* nx)
    (cond
     [(assoc r mr*) => (lambda (r&x) (kk (cdr r&x) mr* fr* nx))]
     [(null? fr*) (kk nx (cons (cons r nx) mr*) fr* (+ nx 1))]
     [else (kk (car fr*) (cons (cons r (car fr*)) mr*) (cdr fr*) nx)]))
  (define (skip v al*)
    (let loop ([r* '()] [al* al*])
      (cond
       [(null? al*) r*]
       [(equal? v (caar al*)) (append (cdr al*) r*)]
       [else (loop (cons (car al*) r*) (cdr al*))])))
  (define (reuse r* mr* fr* nx kk) ;; => (kk mr* fr* nx)
    (let loop ([r* r*] [mr* mr*] [fr* fr*])
      (cond
       [(null? r*) (kk #f mr* fr* nx)]
       [(assoc (car r*) mr*) => (lambda (r&x) (loop (cdr r*) (skip (car r*) mr*) (cons (cdr r&x) fr*)))]
       [else (error 'emit-op "Reusing unused reg")])))
  (define (finish-op ix opt* in* out* mr* fr* nx)
    (let* ([in* (map reg in*)]
           [out* (map reg out*)]
           [instr (case ix
                    [(ref)    (format "{\"g\", ~a, ~a}," (car out*) (mk-idx (car opt*)))]
                    [(adjoin) (format "{\"j\", ~a, ~a}," (car out*) (car in*))]
                    [(real)   (format "{\"r\", ~a, ~a}," (car out*) (car in*))]
                    [(add)    (format "{\"a\", ~a, ~a, ~a}," (car out*) (car in*) (cadr in*))]
                    [(mul)    (format "{\"m\", ~a, ~a, ~a}," (car out*) (car in*) (cadr in*))]
                    [(dot)    (format "{\"d\", ~a, ~a, ~a}," (car out*) (car in*) (cadr in*))]
                    [(shift)  (format "{\"~a\", ~a, ~a, ~a}," (cadr opt*) (car out*) (car in*)
                                      (mk-idx (car opt*)))]
                    ;; die, kill and live instructions
                    ;;[(die)    (and (not (null? opt*)) (format "-- die ~a" (map reg opt*)))]
                    ;;[(kill)   (format "~a = nil;" (reg (car opt*)))]
                    ;;[(live)   (format "-- live [~a] ~a" (length opt*) (map reg opt*))]
                    [(live kill die)   #f]
                    [else  (error 'emit-op "Uknown operation")])])
    (k instr mr* fr* nx)))
  (let* ([opt* (cadr op)]
         [out* (caddr op)]
         [in*  (cadddr op)])
    (if (eq? (car op) 'kill)
        (reuse opt* mr* fr* nx k)
        (let* ([in* (map get-in in*)])
          (let loop ([out* out*] [ro* '()] [mr* mr*] [fr* fr*] [nx nx])
            (cond
             [(null? out*) (finish-op (car op) opt* in* (reverse ro*) mr* fr* nx)]
             [else (get-out (car out*) mr* fr* nx
                            (lambda (x mr* fr* nx)
                              (loop (cdr out*) (cons x ro*) mr* fr* nx)))]))))))

(define (live-analysis rem* res)
  (define (remove-reg* g* live)
    (let loop ([r '()] [l* live])
      (cond
       [(null? l*) r]
       [(member (car l*) g*) (loop r (cdr l*))]
       [else (loop (cons (car l*) r) (cdr l*))])))
  (define (add-reg* g* live)
    (let loop ([r live] [g* g*])
      (cond
       [(null? g*) r]
       [(member (car g*) r) (loop r (cdr g*))]
       [else (loop (cons (car g*) r) (cdr g*))])))
  (define (kill-reg* g* live em*)
    (cond
     [(null? g*) em*]
     [(member (car g*) live) (kill-reg* (cdr g*) live em*)]
     [else (kill-reg* (cdr g*) live (cons ((op-kill (car g*))) em*))]))
  (define (diff-reg* t s)
    (let loop ([r '()] [s s])
      (cond
       [(null? s) r]
       [(member (car s) t) (loop r (cdr s))]
       [(member (car s) (cdr s)) (loop r (cdr s))]
       [else (loop (cons (car s) r) (cdr s))])))
  (let loop ([em* (list ((op-live (list res))))] [rem* rem*] [live (list res)])
    (cond
     [(null? rem*) em*]
     [else (let* ([op (car rem*)]
                  [out* (caddr op)]
                  [in* (cadddr op)]
                  [live-x (remove-reg* out* live)]
                  [em* (kill-reg* in* live-x em*)]
                  [live-y (add-reg* in* live-x)]
                  [em* (cons op em*)]
                  [em* (cons ((op-die (diff-reg* live live-y))) em*)]
                  [em* (cons ((op-live live-y)) em*)])
             (loop em* (cdr rem*) live-y))])))

(define (ops->lua skip f idx op)
  (define (mk-idx x) (cadr (assq x idx)))
  (define (build-nonop op n tbl em* k xm) (build-op op n tbl em* (xm n) k))
  (define (build-unop op n tbl em* k xm)
    (construct-op n tbl (cadr op) em*
      (lambda (r-a n tbl em*) (build-op op n tbl em* (xm n r-a) k))))
  (define (build-binop op n tbl em* k xm)
    (construct-op n tbl (cadr op) em*
      (lambda (r-a n tbl em*)
        (construct-op n tbl (caddr op) em*
          (lambda (r-b n tbl em*) (build-op op n tbl em* (xm n r-a r-b) k))))))
  (define (build-op op n tbl em* xm k)
    (k n (+ n 1) (cons (cons op n) tbl) (cons xm em*)))
  (define (construct-op n tbl op em* k) ;; => (k r n tbl em*)
    (cond
     ;; lookup already computed expressions
     [(assoc op tbl) => (lambda (op&r) (k (cdr op&r) n tbl em*))]
     [else (case (car op)
             [(ref)                  (build-nonop op n tbl em* k (op-ref (cadr op)))]
             [(adjoin)               (build-unop  op n tbl em* k (op-adj))]
             [(real)                 (build-unop  op n tbl em* k (op-real))]
             [(add)                  (build-binop op n tbl em* k (op-add))]
             [(mul)                  (build-binop op n tbl em* k (op-mul))]
             [(dot)                  (build-binop op n tbl em* k (op-dot))]
             [(shift-from-forward)   (build-unop  op n tbl em* k (op-shift (caddr op) "f"))]
             [(shift-from-backward)  (build-unop  op n tbl em* k (op-shift (caddr op) "b"))]
             [else (error 'ops->lua "Uknown operation")])]))
  (define (build-tree op)
    (define (emit-code r n tbl rem*)
      (let loop ([em* (live-analysis rem* r)] [mr* '()] [fr* '()] [nx 1])
        (cond
         [(null? em*) (cdr (assoc r mr*))]
         [else (emit-op f mk-idx (car em*) mr* fr* nx
                        (lambda (x mr* fr* nx)
                          (if x (printf "~a   ~a~%" skip x))
                          (loop (cdr em*) mr* fr* nx)))])))
    (construct-op 1 '() op '() emit-code))
  (let ([r (build-tree op)])
    (printf "~a   {\"v\", ~a}~%" skip r)))

;; print a list of strings
(define (print-list . px)
  (for-each (lambda (x) (printf "~a~%" x)) px))

;; convert path to ops
(define (path->ops pl)
  (cond
   [(null? pl) (error 'path->ops "empty path")]
   [(null? (cdr pl))
    (case (cadar pl)
      [(+) (mk-ref (caar pl))]
      [(-) (mk-adj (mk-shift-from-backward (mk-ref (caar pl)) (caar pl)))])]
   [else (let ([t (path->ops (cdr pl))])
           (case (cadar pl)
             [(+) (mk-mul (mk-ref (caar pl)) (mk-shift-from-forward t (caar pl)))]
             [(-) (mk-shift-from-backward (mk-mul (mk-adj (mk-ref (caar pl))) t) (caar pl))]))]))

;; action piece for a path
(define (path->action pl)
  (let ([a (path->ops (path-head pl))]
        [b (path->ops (path-tail pl))])
    (mk-real (mk-dot a b))))

;; action for a collection of paths
(define (compute-action pl*) (mk-sum (map path->action pl*)))

;; transport matrix of a path
(define (path->matrix pl)
  (let ([a (path->ops (path-head pl))]
        [b (path->ops (path-tail pl))])
    (mk-mul a (mk-adj b))))

;; compute the sum of transport matrices
(define (compute-matrix pl*) (mk-sum (map path->matrix pl*)))

;; force for a single closed path
(define (path->force pl)
  (let ([a (path->ops (path-head pl))]
        [b (path->ops (path-tail pl))])
    (mk-mul a (mk-adj b))))

;; force for a collection of paths
(define (compute-force pl* m)
  (define (add-loop r p)
    (let ([n (count-shifts p m)])
      (let loop ([r r] [p (set-start p m)] [i 0])
        (cond
         [(= i n) r]
         [else (loop (cons p r) (set-start (shift-path p) m) (+ i 1))]))))
  (define (build-paths r p*)
    (if (null? p*) r
        (build-paths (add-loop (add-loop r (car p*))
                               (reverse-path (car p*)))
                     (cdr p*))))
  (mk-sum (map path->force (build-paths '() pl*))))

;; compare two paths possibly with different starting points
(define (same-path? a b)
  (define (equal-path? a)
    (let ([len (length a)])
      (let loop ([i 0] [a a])
        (cond
         [(equal? a b) #t]
         [(= i len) #f]
         [else (loop (+ i 1) (shift-path a))]))))
  (or (equal-path? a)
      (equal-path? (reverse-path a))))

;; symmetrize paths over axes changes
(define (symmetrize-paths m* pl*)
  (let loop-m ([r '()] [m* m*])
    (cond
     [(null? m*) r]
     [else (let loop-p ([r r] [pl* pl*])
             (cond
              [(null? pl*) (loop-m r (cdr m*))]
              [else (loop-p (cons (change-path (car m*) (car pl*)) r)
                            (cdr pl*))]))])))

(define (invariant-set? pl* m*)
  (define (in-set? x pl*)
    (cond
     [(null? pl*) #f]
     [(same-path? x (car pl*)) #t]
     [else (in-set? x (cdr pl*))]))
  (define (inv? m)
    (let loop ([x* (map (lambda (p) (change-path m p)) pl*)])
      (cond
       [(null? x*) #t]
       [(in-set? (car x*) pl*) (loop (cdr x*))]
       [else #f])))
  (cond
   [(null? m*) #t]
   [(inv? (car m*)) (invariant-set? pl* (cdr m*))]
   [else #f]))

(define (unique-paths? pl*)
  (define (in-set? x pl*)
    (cond
     [(null? pl*) #f]
     [(same-path? x (car pl*)) #t]
     [else (in-set? x (cdr pl*))]))
  (cond
   [(null? pl*) #t]
   [(in-set? (car pl*) (cdr pl*)) #f]
   [else (unique-paths? (cdr pl*))]))

;; indices of directions up to 6-d
(define dir-index '([x 0] [y 1] [z 2] [t 3] [s 4] [w 5]))

;; interesting loops
(define plaq1x1    '([x +] [y +] [x -] [y -]))
(define rect1x2    '([x +] [x +] [y +] [x -] [x -] [y -]))
(define chair-a    '([x +] [y +] [z +] [y -] [x -] [z -]))
(define chair-b    '([y +] [x +] [z +] [x -] [y -] [z -]))
(define chair-c    '([x +] [z +] [x -] [y +] [z -] [y -]))
(define chair-d    '([y +] [x -] [z +] [x +] [y -] [z -]))
(define twist-a    '([x +] [z +] [y +] [x -] [z -] [y -]))
(define twist-b    '([y +] [x +] [z +] [y -] [x -] [z -]))
(define twist-c    '([y +] [x -] [z +] [y -] [x +] [z -]))
(define twist-d    '([x +] [y +] [z +] [x -] [y -] [z -]))

(define a1-plaq     (list plaq1x1))
(define a1-rect     (list rect1x2))
(define a1-chair    (list chair-a chair-b chair-c chair-d))
(define a1-twist    (list twist-a twist-b twist-c twist-d))

(define a2-plaq     (map double-path a1-plaq))
(define a2-rect     (map double-path a1-rect))
(define a2-chair    (map double-path a1-chair))
(define a2-twist    (map double-path a1-twist))

(define a1-clover     '(([x +] [y +] [x -] [y -])
                        ([y +] [x -] [y -] [x +])
                        ([x -] [y -] [x +] [y +])
                        ([y -] [x +] [y +] [x -])))
(define a2-clover    (map double-path a1-clover))
(define a3-clover    (map triple-path a1-clover))

(define (generate-loops lname fname gname)
  (define tname (format "~a.~a" lname fname))
  (define (gen-loop rx)
    (let ([name (car rx)]
          [loop (cadr rx)]
          [coords (caddr rx)]
          [derivs (cadddr rx)]
          [generator (car (cddddr rx))])
      (define (gen-force dx)
        (printf "     d~a = {~%" dx)
        (ops->lua "      " "U" dir-index (compute-force loop (list dx '+)))
        (printf "     },~%"))
      (printf "  ~a = ~a.~a.~a(\"~a\", {~%" name lname gname generator name)
      (printf "   action = {~%")
      (ops->lua "    " "U" dir-index (compute-action loop))
      (printf "   },~%")
      (printf "   force = {~%")
      (for-each gen-force derivs)
      (printf "   },~%")
      (printf "   count = ~a~%" (length loop))
      (printf "   }),~%")))
  (define (gen-matrix rx)
    (let ([name (car rx)]
          [loop (cadr rx)]
          [coords (caddr rx)]
          [generator (cadddr rx)])
      (printf "  ~a = ~a.~a.~a(\"~a\", {~%" name lname gname generator name)
      (printf "   matrix = {~%")
      (ops->lua "    " "U" dir-index (compute-matrix loop))
      (printf "   }}),~%")))
  (printf "-- This is magic. You have been warned.~%")
  (printf "~a = ~a or {};~%" lname lname)
  (printf "~a = {~%" tname)
  (for-each gen-loop
            (list (list "plaq1"   a1-plaq "x,y"     '[x]     "make_sym_xy")
                  (list "rect1"   a1-rect "x,y"     '[x y]   "make_gen_xy")
                  (list "chair1"  a1-chair "x,y,z"  '[x z]   "make_mix_xyz")
                  (list "twist1"  a1-twist "x,y,z"  '[x]     "make_sym_xyz")
                  (list "plaq2"   a2-plaq "x,y"     '[x]     "make_sym_xy")
                  (list "rect2"   a2-rect "x,y"     '[x y]   "make_gen_xy")
                  (list "chair2"  a2-chair "x,y,z"  '[x z]   "make_mix_xyz")
                  (list "twist2"  a2-twist "x,y,z"  '[x]     "make_sym_xyz")))
  (for-each gen-matrix
            (list (list "clover1x1"  a1-clover "x,y"  "make_matrix_xy")
                  (list "clover2x2"  a2-clover "x,y"  "make_matrix_xy")
                  (list "clover3x3"  a3-clover "x,y"  "make_matrix_xy")))
  (printf "};~%"))


(define (generate-collection lname tname gname l*)
  (define (save-element root gen record)
    (let ([fname (format "~a.qlua" root)])
      (printf ";;; saving to ~a ...~%" fname)
      (with-output-to-file fname (lambda () (gen record)))))
  (define (gen-loop rx)
    (let ([name (car rx)]
          [loop (cadr rx)]
          [coords (caddr rx)]
          [derivs (cadddr rx)]
          [generator (car (cddddr rx))])
      (define (gen-force dx)
        (printf "     d~a = {~%" dx)
        (ops->lua "      " "U" dir-index (compute-force loop (list dx '+)))
        (printf "     },~%"))
      (printf "-- This is magic. You have been warned.~%")
      (printf "require \"qcdlib/gauge-loops\";~%")
      (printf "~a.~a.~a = ~a.~a.~a(\"~a\", {~%" lname tname name lname gname generator name)
      (printf "   action = {~%")
      (ops->lua "    " "U" dir-index (compute-action loop))
      (printf "   },~%")
      (printf "   force = {~%")
      (for-each gen-force derivs)
      (printf "   },~%")
      (printf "   count = ~a~%" (length loop))
      (printf "   });~%")))
  (define (gen-matrix rx)
    (let ([name (car rx)]
          [loop (cadr rx)]
          [coords (caddr rx)]
          [generator (cadddr rx)])
      (printf "-- This is magic. You have been warned.~%")
      (printf "require \"qcdlib/gauge-loops\";~%")
      (printf "~a.~a.~a = ~a.~a.~a(\"~a\", {~%" lname tname name lname gname generator name)
      (printf "   matrix = {~%")
      (ops->lua "    " "U" dir-index (compute-matrix loop))
      (printf "   }});~%")))
  (define (gen-element record)
    (case (car record)
      [(gen-loop) (save-element (cadr record) gen-loop (cdr record))]
      [(gen-matrix) (save-element (cadr record) gen-matrix (cdr record))]
      [else (error 'gen-element (format "unknown kind of record ~a "record))]))
  (for-each gen-element l*))

(define collection
  (list (list 'gen-loop "plaq1"   a1-plaq "x,y"     '[x]     "make_sym_xy")
        (list 'gen-loop "rect1"   a1-rect "x,y"     '[x y]   "make_gen_xy")
        (list 'gen-loop "chair1"  a1-chair "x,y,z"  '[x z]   "make_mix_xyz")
        (list 'gen-loop "twist1"  a1-twist "x,y,z"  '[x]     "make_sym_xyz")
        (list 'gen-loop "plaq2"   a2-plaq "x,y"     '[x]     "make_sym_xy")
        (list 'gen-loop "rect2"   a2-rect "x,y"     '[x y]   "make_gen_xy")
        (list 'gen-loop "chair2"  a2-chair "x,y,z"  '[x z]   "make_mix_xyz")
        (list 'gen-loop "twist2"  a2-twist "x,y,z"  '[x]     "make_sym_xyz")
        (list 'gen-matrix "clover1x1"  a1-clover "x,y"  "make_matrix_xy")
        (list 'gen-matrix "clover2x2"  a2-clover "x,y"  "make_matrix_xy")
        (list 'gen-matrix "clover3x3"  a3-clover "x,y"  "make_matrix_xy")))


;; examples
;; (load "paths.ss")

;;;;;; misc examples
;; (ops->lua "" "U" "t" (compute-action a1-plaq))
;; (ops->lua "" "U" "t" (compute-action a1-rect))
;; (ops->lua "" "U" "t" (compute-force a1-chair '(x +)))
;; (ops->lua "" "U" "f" (compute-force a2-twist '(x +)))
;; (with-output-to-file "qcdlib/loop-tables.qlua" (lambda () (generate-loops "qcdLib" "LoopTables" "GaugeLoops")))

;;;;;;
;; The code below generates a set of files for various gauge loops.
;; Most of the files in $TOP/qlib/qcdlib/gauge-loops/ are produced by the
;; expression below.
;; (generate-collection "qcdLib" "LoopTables" "GaugeLoops" collection)
